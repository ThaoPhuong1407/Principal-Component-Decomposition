(ns pcd.core
  (:require
   [clojure.pprint :as pp]
   [clojure.core.matrix :as matrix]
   [clojure.core.matrix.stats :as matrix-stats]
   [clojure.core.matrix.dataset :as matrix-ds]
   [clojure.core.matrix.linear :as matrix-linear]
   [clojure.core.matrix.protocols :as matrix-proto]
   [clojure.core.matrix.operators :as matrix-operator]
   [incanter.core :as in-core]
   [incanter.stats :as in-stats]
   [clojure.java.io :as io]
   [clojure.data.csv :as csv])
  (:gen-class))

;; VARIABLES
(def infile-path "../../data/Mall_Customers.csv")

;; SUPPORTING FUNCTIONS
(defn log2 [n] (/ (Math/log n) (Math/log 2)))

(defn vec-remove
  "remove an elem from a vector. O(1)"
  [vec position]
  (into (subvec vec 0 position) (subvec vec (inc position))))

(defn process-data
  "Arguments: 
   • filename: the path to input file
   • columns (optional): columns to remove from the original data
   Return data as clojure.core.matrix.impl.dataset.DataSet"
  ([filename columns]
   (with-open [reader (io/reader filename)]
     (let [data (drop 1 (csv/read-csv reader))
           ds (matrix-ds/dataset data)
           ds-drop (matrix-ds/remove-columns ds columns)]
       ds-drop)))
  ([filename]
   (with-open [reader (io/reader filename)]
     (let [data (drop 1 (csv/read-csv reader))
           ds (matrix-ds/dataset data)]
       ds))))

(defn write-data-csv
  [data]
  (with-open [writer (io/writer "out-file.csv" :append true)]
    (csv/write-csv writer data)))

(defn matrix-float
  "Arguments: 
   • len: number of columns
   • ds: dataset
   Return: a matrix containing only float datatype"
  [len ds]
  (loop
   [index 0
    result []]
    (if (< index len)
      (recur (inc index) (into result [(vec (map #(Float/parseFloat %) (get ds index)))]))
      result)))

(defn get-dist-from-plane-to-origin
  "Note: The hyperplane starts at the mean point. It goes through the mean, and is perpendicular to the eigenvector"
  [mean eigenvector]
  (matrix-proto/vector-dot mean eigenvector))

(defn get-dist-from-plane-to-points
  [dist-from-plane-to-origin data-set eigenvector]
  (let [len-eigenvector (count (vec eigenvector))
        dist-from-points-to-origin (matrix-proto/vector-dot data-set (matrix/reshape eigenvector [len-eigenvector 1]))
        dist-from-plane-to-points (matrix-operator/- dist-from-plane-to-origin dist-from-points-to-origin)
        len-dist (count dist-from-plane-to-points)]
    (matrix/reshape dist-from-plane-to-points [len-dist])))

(defn description-length
  [distance-array]
  (let [num-elements (count distance-array)
        variance (cond
                   (> num-elements 1) (matrix-stats/variance distance-array)
                   (= num-elements 1) 0
                   :else Double/NaN)]
    (cond
      (= variance 0) 20
      (Double/isNaN variance) 0
      :else (+ 20 (* num-elements (log2 (Math/sqrt variance)))))))

;; MAIN 
;; 1. Process data
(def ds (process-data infile-path (list 0 1 2))) ;; import data and remove columns 0 1 2
(def ds-float (matrix-float (count ds) ds)) ;; convert matrix of string to float

;; 2. covariance matrix, eigenvalues and eigenvectors(
;; (def cov-matrix (in-stats/covariance ds-float))
(def svd (matrix-linear/svd ds-float))
(def eigen-values (vec (get svd :S)))
(def eigen-vectors-T (vec (get svd :V*)))

;; 3. Sort eigenvectors in the descending order of eigenvalues.
(def with-idx (map-indexed vector eigen-values))
(def sorted-eigen-val (sort-by last > with-idx))
(def sorted-idx (map #(first %) sorted-eigen-val))
(def sorted-eigen-vectors (map #(get eigen-vectors-T %) sorted-idx))

;; 4. Chopping process
(defn divide-cluster-using-DL
  [data-and-dist-vector]
  "Find subclusters if any from a dataset and its distance array"
  (let [result-cut-point (loop
                          [temp1 data-and-dist-vector
                           temp2 []
                           minDL (description-length (into [] (map last data-and-dist-vector)))
                           cut-point 0
                           dl-arr [minDL]
                           index 0]

                           (if (< index (count data-and-dist-vector))

                             (let [curr-item (nth data-and-dist-vector index)
                                   data-and-dist-1 (remove #{curr-item} temp1)
                                   data-and-dist-2 (conj temp2 curr-item)
                                   dist-1 (into [] (map last data-and-dist-1))
                                   dist-2 (into [] (map last data-and-dist-2))
                                   new-dl (+ (description-length dist-1) (description-length dist-2))
                                   cut-point-temp (if (< new-dl minDL) (+ 1 index) cut-point)]

                               (recur data-and-dist-1
                                      data-and-dist-2
                                      (if (< new-dl minDL) new-dl minDL)
                                      cut-point-temp
                                      (conj dl-arr new-dl)
                                      (inc index)))
                             [cut-point, dl-arr]))]

    (if (= (nth result-cut-point 0) 0)
      [(vec (map first data-and-dist-vector))] ;; return an vector of clusters 
      (let [sub-cluster1 (subvec data-and-dist-vector 0 (nth result-cut-point 0))
            sub-cluster2 (subvec data-and-dist-vector (nth result-cut-point 0))]
        (concat (divide-cluster-using-DL sub-cluster1) (divide-cluster-using-DL sub-cluster2))))))

(defn find-subclusters-from-data
  "find subclusters, given a eigenvector and a dataset"
  [eigenvector data]
  (let
   [;; a. Find the distance of the hyperplane, which goes through data mean, to origin
    mean (matrix-stats/mean data)
    dist-from-plane-to-origin (get-dist-from-plane-to-origin mean eigenvector)

    ;; b. Compute the list of distances from all points to the cutting plane
    dist-from-plane-to-points (get-dist-from-plane-to-points dist-from-plane-to-origin data eigenvector)
    dist-from-plane-to-points-dict (zipmap data dist-from-plane-to-points)

    ;; c. Sort the points in order of distance from the cutting hyper plane (positive to negative)
    sorted-dist-from-plane-to-points (vec (sort-by val > dist-from-plane-to-points-dict))

    ;; d. Find subclusters if any
    sub-clusters (divide-cluster-using-DL sorted-dist-from-plane-to-points)]
    ;; (println sub-clusters, (count sub-clusters))

    sub-clusters))

(defn chop-cluster
  [sorted-vectors, data]
  "find subclusters, given an array of sorted eigenvectors and a dataset"
  (loop
   [index-vec 0
    clusters []]
    (if (< index-vec (count sorted-vectors))
      (let
       [current-eigenvector (nth sorted-vectors index-vec)
        sub-clusters (if (= index-vec 0)
                       ;; 1st vector
                       (find-subclusters-from-data current-eigenvector data)

                       ;; remaining vectors 
                       (loop
                        [index-cluster 0
                         result []]
                         (if (< index-cluster (count clusters))
                           (let [curr-cluster (nth clusters index-cluster)
                                 sub-curr-clusters (find-subclusters-from-data current-eigenvector curr-cluster)]

                             (recur (inc index-cluster)
                                    (into result sub-curr-clusters)))
                           result)))]

        (println current-eigenvector, (count sub-clusters))
        (recur (inc index-vec)
               sub-clusters))
      clusters)))

(def clusters (chop-cluster sorted-eigen-vectors ds-float))

;; 5. Merging process
(defn mahalanobis-distance
  [point mean invcovar]
  ;; sqrt[(point - mean).T * invcovar * (point - mean)]

  (matrix/mget (matrix/inner-product
                (matrix/sub point mean)
                invcovar
                (matrix/transpose (matrix/sub point mean)))))

(defn euclidean-distance
  [point1 point2]
  (Math/sqrt
   (reduce + (map #(Math/pow (- %1 %2) 2) point1 point2))))

(defn cluster-assignment-cost
  [length-1 length-2]
  (if (or (= length-1 0)
          (= length-2 0))
    0
    (let [total-length (+ length-1 length-2) ;; 19 
          probability-1 (/ length-1 total-length) ;; 10/19  -0.926
          probability-2 (/ length-2 total-length)] ;; 9/19  -1.078  
      ;; - [(total-length) * [(prob-1 * log2(prob-1)) + (prob-2 * log2(prob-2))]]
      (- (* total-length
           (+
            (* probability-1 (log2 probability-1))
            (* probability-2 (log2 probability-2))))))))

(defn merge-2-clusters-using-DL
  "Return true if we should merge 2 clusters. False otherwise."
  [cluster1 cluster2]
  (let [length-c1 (count cluster1)
        length-c2 (count cluster2)
        cac (cluster-assignment-cost length-c1 length-c2) ;; cost of asignment cost
        merged-cluster (vec (concat cluster1 cluster2)) ;; merge 2 clusters

        ;; get the sorted eigenvectors
        svd (matrix-linear/svd merged-cluster)
        eigen-values (vec (get svd :S))
        eigen-vectors-T (vec (get svd :V*))
        with-idx  (map-indexed vector eigen-values)
        sorted-eigen-val (sort-by last > with-idx)
        sorted-idx (map #(first %) sorted-eigen-val)
        sorted-eigen-vectors (map #(get eigen-vectors-T %) sorted-idx)
        mean (matrix-stats/mean merged-cluster)]
    (loop
     [eigen-index 0
      return-val false]
      (if (< eigen-index (count eigen-vectors-T))
        (let [curr-eigenvector (nth sorted-eigen-vectors eigen-index)
              dist-from-plane-to-origin (get-dist-from-plane-to-origin mean curr-eigenvector)
              dist-from-plane-to-points-c1 (get-dist-from-plane-to-points dist-from-plane-to-origin cluster1 curr-eigenvector)
              dist-from-plane-to-points-c2 (get-dist-from-plane-to-points dist-from-plane-to-origin cluster2 curr-eigenvector)
              dist-from-plane-to-points-merged (get-dist-from-plane-to-points dist-from-plane-to-origin merged-cluster curr-eigenvector)
              dl-1 (description-length dist-from-plane-to-points-c1)
              dl-2 (description-length dist-from-plane-to-points-c2)
              dl-1plus2 (+ dl-1 dl-2 cac) ;; No cac here -- still testing
              dl-merged (description-length dist-from-plane-to-points-merged)]
          
          (if (< dl-merged dl-1plus2)   
           (do 
             ;; (println "cluster1 + cluster2: " dl-1plus2 "cac: " cac)
             ;; (println "cluster1 AND cluster2 : " dl-merged)
             true)
            (recur (inc eigen-index) false)))
        return-val))))

(defn sort-by-distance
  [cluster clusters]
  (let [cluster-mean (matrix-stats/mean cluster)]
    (sort-by #(if (<= (count %) 2)
                (euclidean-distance cluster-mean (matrix-stats/mean %))
                (mahalanobis-distance cluster-mean
                                      (matrix-stats/mean %)
                                      (matrix/inverse (in-stats/covariance %))))
             <
             clusters)))

#_(reduce works like this
        (reduce (fn [first rest]
                 do something
                 if (break condition) (reduced return value) ;; return immediately, works like "break" or "return" in JS
                 else first = whatever value passed in here) array) ;; The value of "first" will be replaced here
      
          Each iteration
          - "first" will be replaced by the last value in the bracket
          - "rest" is the next element in the array)

(defn merge-clusters
  [clusters]

  (let [;; 1. Sort clusters based on their lengths
        sorted-clusters (vec (sort-by count > clusters))
        length (count sorted-clusters)]
     
    ;; 2. Loop through the each cluster vs rest:
    (loop 
       [index 0
        result []]
        (if (< index length)
          (let [current-cluster (vec (nth sorted-clusters index)) ;; the current-cluster in this iteration
                rest-clusters (subvec sorted-clusters (+ index 1)) ;; all clusters after the current-cluster in this iteration (no repetition)

                ;; 2a. Sort rest-clusters based on either the euclidean
                ;; or mahalanobis distance in relative to the current cluster -- To speed up the merging process
                sorted-rest-clusters (if (>= length 2)
                                       (vec (sort-by-distance current-cluster rest-clusters)) ;; Only sort if we have more than 2 clusters
                                       rest-clusters) ;; No need to sort if we only have 1 cluster.

                ;; 2b. Try to merge the current-cluster with clusters in the sorted-rest-clusters 
                joined-clusters (vec (concat [current-cluster] sorted-rest-clusters)) 
                sorted-cluster-idx (atom 0) ;; If a merge happens, we can remove the cluster from sorted-rest-clusters using this index
                output (reduce (fn [first-cluster sorted-cluster] 
                                   (if (= true (merge-2-clusters-using-DL first-cluster sorted-cluster))
                                     ;; If a merge happens, return the new vector of clusters, where the merged clusters is added, and the 2 clusters are removed
                                     (do 
                                       (println "Merge with: " sorted-cluster)
                                       ;; (println (subvec sorted-clusters 0 (+ index 1)))
                                       ;; (println [(vec (concat first-cluster sorted-cluster))])
                                       ;; (println (vec-remove sorted-rest-clusters @sorted-cluster-idx))

                                       (reduced
                                        (concat
                                         (vec-remove sorted-rest-clusters @sorted-cluster-idx) ;; remove the sorted-cluster from sorted-rest-clusters    
                                         [(vec (concat first-cluster sorted-cluster))]         ;; the merged cluster
                                         (subvec sorted-clusters 0 index))))             ;; all clusters iterated (not included in sorted-rest-clusters)
                                     ;; Else, do nothing, the last iteration will return first-cluster, which should be "current-cluster"
                                     (do
                                       (swap! sorted-cluster-idx inc)
                                       first-cluster))) joined-clusters)]

            (if (= output current-cluster)
              ;; No merge
              (do
                (println "NO MERGED" (count sorted-clusters))
                (recur (inc index) ;; move onto the next cluster
                       clusters))  ;; result = original clusters since no merge
              ;; A merge happened
              (do
                (println "MERGED" (count sorted-clusters) (count output))
                ;; (println "Vector before Merge")
                ;; (pp/pprint sorted-clusters)
                ;; (println "Vector after Merge")
                ;; (pp/pprint output)
                (merge-clusters (vec output))))) ;; recursively call merge-clusters with the new vector
          result))))

;; TESTING
(def merged-clusters (merge-clusters clusters))

(write-data-csv (map #(into ["x", "y"] %) merged-clusters))
(map #(write-data-csv (into [["empty", "x", "y"]] %)) merged-clusters)


(defn -main
  "I don't do a whole lot ... yet."
  [& args])



