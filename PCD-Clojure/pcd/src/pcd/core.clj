(ns pcd.core
  (:require
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

;; MAIN 
;; 1. Process data
(def ds (process-data infile-path (list 0 1 2))) ;; import data and remove columns 0 1 2
(def ds-float (matrix-float (count ds) ds)) ;; convert matrix of string to float

;; 2. covariance matrix, eigenvalues and eigenvectors(
(def cov-matrix (in-stats/covariance ds-float))
(def svd (matrix-linear/svd ds-float))
(def eigen-values (vec (get svd :S)))
(def eigen-vectors-T (vec (get svd :V*)))

;; 3. Sort eigenvectors in the descending order of eigenvalues.
(def with-idx (map-indexed vector eigen-values))
(def sorted-eigen-val (sort-by last > with-idx))
(def sorted-idx (map #(first %) sorted-eigen-val))
(def sorted-eigen-vectors (map #(get eigen-vectors-T %) sorted-idx))

;; 4. Chopping process
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
        variance (if (> num-elements 0) (matrix-stats/variance distance-array) Double/NaN)]
    (cond
      (= variance 0) 20
      (Double/isNaN variance) 0
      :else (+ 20 (* num-elements (log2 (Math/sqrt variance)))))))


;;  RESULT IS NOT CORRECT, DO THIS!
(defn divide-clusters-using-DL
  [data-and-dist-vector]
  (let [result-cut-point (loop
                          [index 0
                           temp1 data-and-dist-vector
                           temp2 []
                           minDL (description-length (into [] (map last data-and-dist-vector)))
                           cut-point 0
                           dl-arr [minDL]]

                           (if (< index (count data-and-dist-vector))

                             (let [curr-item (nth data-and-dist-vector index)
                                   data-and-dist-1 (remove #{curr-item} temp1)
                                   data-and-dist-2 (conj temp2 curr-item)
                                   dist-1 (into [] (map last data-and-dist-1))
                                   dist-2 (into [] (map last data-and-dist-2))
                                   new-dl (+ (description-length dist-1) (description-length dist-2))
                                   cut-point-temp (if (< new-dl minDL) (+ 1 index) cut-point)]

                               (recur (inc index)
                                      data-and-dist-1
                                      data-and-dist-2
                                      (if (< new-dl minDL) new-dl minDL)
                                      cut-point-temp
                                      (conj dl-arr new-dl)))
                             [cut-point, dl-arr]))]
    ;; (println "dl-arr" (nth result-cut-point 1))

    (if (= (nth result-cut-point 0) 0)
      (do
        (println data-and-dist-vector)
        data-and-dist-vector)

      (let [sub-cluster1 (subvec data-and-dist-vector 0 (nth result-cut-point 0))
            sub-cluster2 (subvec data-and-dist-vector (nth result-cut-point 0))]
        (conj (divide-clusters-using-DL sub-cluster1) (divide-clusters-using-DL sub-cluster2))))))


(defn find-subcluster-from-data
  "find subclusters, given a eigenvector and a dataset"
  [eigenvector data]
  (let
   [;; a. Find the distance of the hyperplane, which goes through data mean, to origin
    mean (matrix-stats/mean data)
    dist-from-plane-to-origin (get-dist-from-plane-to-origin mean eigenvector)

    ;; b. Compute the list of distances from all points to the cutting plane
    dist-from-plane-to-points (get-dist-from-plane-to-points dist-from-plane-to-origin data eigenvector)
    dist-from-plane-to-points-dict (zipmap ds-float dist-from-plane-to-points)

    ;; c. Sort the points in order of distance from the cutting hyper plane (positive to negative)
    sorted-dist-from-plane-to-points (vec (sort-by val > dist-from-plane-to-points-dict))
    sub-clusters (divide-clusters-using-DL sorted-dist-from-plane-to-points)]
    (println "eigenvector " eigenvector)
    (println "sub-clusters " sub-clusters)
    (println "sub-clusters length " (count sub-clusters))
    sub-clusters))

(defn chop
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
                       (find-subcluster-from-data current-eigenvector data)
                       ;; remaining vectors 
                       (loop
                        [index-cluster 0
                         result []]
                         (if (< index-cluster (count clusters))
                           (let [curr-cluster (nth clusters index-cluster)
                                 curr-data (if (= 2 (count (matrix/shape curr-cluster)))
                                             (conj [] (first curr-cluster))
                                             (into [] (vec (map first curr-cluster))))]
                            ;;  (println curr-cluster)
                             (recur (inc index-cluster)
                                    (into result (find-subcluster-from-data current-eigenvector curr-data))))
                           result)))]

        (recur (inc index-vec)
               sub-clusters))
      clusters)))


(def sub-clusters (find-subcluster-from-data (nth sorted-eigen-vectors 0) ds-float))

;; (def clusters (chop sorted-eigen-vectors ds-float))
;; (def sub-clusters (divide-clusters-using-DL clusters))



;;    b. Get a list of distances of all points to the hyperplane

;;       (def dist_hyperplane_point (matrix-operator/- dist_hyperplane_origin (matrix-proto/vector-dot ds-float e-vec)))
;;    3. Sort distance_mean_point

;;     # c) Sort distance_mean_point
;;     sorted_index = np.argsort(dist_hyperplane_point)[::-1]
;;     sorted_dist_hyperplane_point = dist_hyperplane_point[sorted_index]
;;     sorted_data = data[sorted_index]
;;     return [sorted_dist_hyperplane_point, sorted_data]

(defn -main
  "I don't do a whole lot ... yet."
  [& args])



