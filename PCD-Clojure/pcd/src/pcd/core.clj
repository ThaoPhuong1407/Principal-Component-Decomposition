(ns pcd.core
  (:require
   [clojure.pprint :as pp]
   [clojure.java.io :as io]
   [clojure.data.csv :as csv]

   ;; stats
   [clojure.core.matrix :as matrix]
   [clojure.core.matrix.stats :as matrix-stats]
   [clojure.core.matrix.dataset :as matrix-ds]
   [clojure.core.matrix.linear :as matrix-linear]
   [clojure.core.matrix.protocols :as matrix-proto]
   [clojure.core.matrix.operators :as matrix-operator]
   [incanter.core :as in-core]
   [incanter.stats :as in-stats]

   ;; Python
   [libpython-clj2.require :refer [require-python]]
   [libpython-clj2.python :refer [py. py.. py.-] :as py]

   ;; display the plot
   [clojure.java.shell :as sh])
  (:gen-class))

;; Matplotlib 
(require-python '[matplotlib.pyplot :as plot])
(require-python '[matplotlib.cm :as color])

#_(defmacro with-show
    "Takes forms with mathplotlib.pyplot to then show locally"
    [save-file & body]
    `(let [_# (plot/clf) ;; clear the current figure
           fig# (plot/figure)]
       ~(cons 'do body) ;; what we want to plot
       (plot/savefig ~save-file)))

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
  [data out-file-name]
  (with-open [writer (io/writer (str "./src/pcd/data-output/" out-file-name) :append true)]
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

(defn euclidean-distance
  [point1 point2]
  (Math/sqrt
   (reduce + (map #(Math/pow (- %1 %2) 2) point1 point2))))

(defn get-dist-from-plane-to-origin
  "Note: The hyperplane starts at the mean point. It goes through the mean, and is perpendicular to the eigenvector"
  [mean eigenvector]
  ;; (println "here" mean eigenvector)
  ;; (println (matrix-proto/vector-dot mean eigenvector))
  (matrix-proto/vector-dot mean eigenvector))

(defn get-dist-from-plane-to-points
  [dist-from-plane-to-origin data-set eigenvector]
  (let [len-eigenvector (count (vec eigenvector))
        dist-from-points-to-origin (matrix-proto/vector-dot data-set (matrix/reshape eigenvector [len-eigenvector 1]))
        dist-from-plane-to-points (matrix-operator/- dist-from-plane-to-origin dist-from-points-to-origin)
        len-dist (count dist-from-plane-to-points)]
    (matrix/reshape dist-from-plane-to-points [len-dist])))

(defn get-dist-from-plane-to-points-V2

  "References: https://www.khanacademy.org/math/linear-algebra/vectors-and-spaces/dot-cross-products/v/point-distance-to-plane
   - plane equation: ax + by + cz = d, input is in form: [a, b, c, -d]
   - normal vec = [a b c]
   - Distance from a point to a plane = (dot(point, normalVec) - d) / normalVecLength "

  [planeEquation points]
  (let [normVec (butlast planeEquation) ;; a, b, c
        d (last planeEquation) ;; -d
        arrayOfZero (make-array Double/TYPE (count normVec))
        normalVecLength (euclidean-distance normVec arrayOfZero)
        distanceList (loop [index 0
                            distances []]

                       ;; Break-out condition 
                       (if (>= index (count points))
                         ;; return the distance list
                         distances

                         ;; calculate the distance from point to plane
                         (let [curr-point (nth points index)
                               dotVal (+ (matrix-proto/vector-dot curr-point normVec) d)
                               distance (/ dotVal normalVecLength)]
                           (recur (inc index)
                                  (conj distances distance)))))]
    distanceList))

(defn get-hyperplane-equation
  "References: https://math.stackexchange.com/questions/82151/find-the-equation-of-the-plane-passing-through-a-point-and-a-vector-orthogonal
   The return value is Ax + By = dot(mean, eigenvector) = [A, B, -dot(mean,eigenvector)]"
  [mean eigenvector]
  (let [dotProd  (- (matrix-proto/vector-dot mean eigenvector))
        equation (conj (vec eigenvector) dotProd)]
    equation))

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

(defn cluster-assignment-cost
  [length-1 length-2]
  (if (or (= length-1 0)
          (= length-2 0))
    0
    (let [total-length (+ length-1 length-2)
          probability-1 (/ length-1 total-length)
          probability-2 (/ length-2 total-length)]
      ;; Formula: -[(total-length) * [(prob-1 * log2(prob-1)) + (prob-2 * log2(prob-2))]]
      (- (* total-length
            (+
             (* probability-1 (log2 probability-1))
             (* probability-2 (log2 probability-2))))))))

(defn get-sorted-eigen-vectors
  [data]
  (let [;; 1. Find eigenvalues and eigenvectors using SVD
        cov-matrix (in-stats/covariance data)
        svd (matrix-linear/svd cov-matrix) ;; 2x2
        eigen-values (vec (get svd :S))
        eigen-vectors-T (vec (get svd :V*))

        ;; 2. Sort eigenvectors in the descending order of eigenvalues.
        assign-idx (map-indexed vector eigen-values)
        sorted-eigen-val (sort-by last > assign-idx)
        sorted-idx (map #(first %) sorted-eigen-val)
        sorted-eigen-vectors (map #(get eigen-vectors-T %) sorted-idx)]
    sorted-eigen-vectors))

(defn get-descending-sort-dist
  "Return a descending sorted list of distances between points to hyperplane"
  [eigenvector data is-plot sorted-eigenvecs]
  (let
   [;; a. Construct the hyperplane equation, that goes through the mean
    mean (matrix-stats/mean data)
    hyperplane (get-hyperplane-equation mean eigenvector)

    ;; b. Compute the list of distances from all points to the cutting plane
    dist-from-plane-to-points (get-dist-from-plane-to-points-V2 hyperplane data)
    dist-from-plane-to-points-dict (zipmap data dist-from-plane-to-points)

    ;; c. Sort the points in order of distance from the cutting hyper plane (positive to negative)
    sorted-dist-from-plane-to-points (vec (sort-by val > dist-from-plane-to-points-dict))

    ;; For drawing: 2 eigenvectors
    first-eigen (first sorted-eigenvecs)
    second-eigen (second sorted-eigenvecs)]

    ;;;;;-------- VISUALIZATION 
    ;; a. Initialized figure with 3 subplots
    (plot/subplots :figsize [15 5])

    (when is-plot
    ;; b. Drawing the 1st subplot 
      (plot/subplot 1 3 1 :adjustable "datalim" :aspect 1)
    ;; Draw dataset
      (plot/scatter  (into [] (map first data)) (into [] (map last data)) :c "black")
    ;; Draw 2 eigenvectors going through the mean point
      (plot/quiver (first mean) (second mean) (first first-eigen) (second first-eigen) :color "pink" :scale 2)
      (plot/quiver (first mean) (second mean) (first second-eigen) (second second-eigen) :color "blue" :scale 2)
    ;; Set legend and subplot title
      (plot/legend ["data" "1st eigenvec" "2nd eigenvec"])
      (plot/title "Data with Eigenvecs"))

    sorted-dist-from-plane-to-points))

(defn draw-final
  [clusters]
  (let [len (count clusters)
        colors-chopped (color/get_cmap "tab20" len)]
    (loop
     [i 0]
      (when (< i len)
            ;; loop 
        (let [x-vector (into [] (map first (nth clusters i)))
              y-vector (into [] (map last (nth clusters i)))
              color [(colors-chopped i)]]
          (plot/scatter x-vector y-vector :c color)
          (recur (inc i)))))))

;; 4. Chopping process
(defn maybe-divide-cluster-using-DL
  "Divide a cluster into 2 sub-clusters if makes sense
   The input is in form of [[data1, distance1], [data2, distance2], etc.]. 
      Reason: Although we return the actual data set, we need the distances for calculation"
  [data-and-dist-vector eigenvector is-plot pic-name]

  ;; DIVIDING PROCESS
  ;; result-cut-point = [cut-point, dl-arr]
  (let [data (map first data-and-dist-vector)
        mean (matrix-stats/mean data)
        ;; slope (/ (nth eigenvector 1) (nth eigenvector 2))
        result-cut-point (loop
                          [temp1 data-and-dist-vector
                           temp2 []
                           minDL (description-length (into [] (map last data-and-dist-vector)))
                           cut-point 0
                           dl-arr [minDL] ;; For visualization, which is not yet implemented for Clojure
                           index 0]

                           (if (< index (count data-and-dist-vector))

                             (let [curr-item (nth data-and-dist-vector index)
                                   data-and-dist-1 (remove #{curr-item} temp1)
                                   data-and-dist-2 (conj temp2 curr-item)
                                   dist-1 (into [] (map last data-and-dist-1))
                                   dist-2 (into [] (map last data-and-dist-2))
                                   new-dl (+ #_(cluster-assignment-cost (count dist-1) (count dist-2)) (description-length dist-1) (description-length dist-2))
                                   cut-point-temp (if (< new-dl minDL) (+ 1 index) cut-point)]
                              ;;  (println minDL new-dl)
                               (recur data-and-dist-1
                                      data-and-dist-2
                                      (if (< new-dl minDL) new-dl minDL)
                                      cut-point-temp
                                      (conj dl-arr new-dl)
                                      (inc index)))

                             [cut-point, dl-arr]))
        cut-point (nth result-cut-point 0)
        dl-arr (nth result-cut-point 1)
        cluster-1 (subvec data-and-dist-vector 0 cut-point)
        cluster-2 (subvec data-and-dist-vector cut-point)
        cluster-data-1 (vec (map first cluster-1))
        cluster-data-2 (vec (map first cluster-2))]

    ;;;;; SAVE DL DATA
    (write-data-csv [dl-arr] "description-len-chop.csv")

    ;;;;;-------- VISUALIZATION 
    (when is-plot
      ;; Drawing the 2nd subplot: Visualize the description length 
      (let [data-y (nth result-cut-point 1)
            data-x (range (count data-y))]
        (plot/subplot 1 3 2 :adjustable "datalim" :aspect 1 )
        (plot/plot data-x data-y)
        (plot/title "DL"))


      ;; Drawing the 3rd subplot: hyperplane and 2 clusters if cut
      (let [cluster1X (into [] (map first cluster-data-1)) ;; x
            cluster1Y (into [] (map last cluster-data-1))  ;; y
            cluster2X (into [] (map first cluster-data-2))
            cluster2Y (into [] (map last cluster-data-2))]

        (plot/subplot 1 3 3 :adjustable "datalim" :aspect 1)

        ;; Draw points in 2 sub-clusters (2 different colors)
        (plot/scatter cluster1X cluster1Y :c "red")
        (plot/scatter cluster2X cluster2Y :c "blue")

        ;; Draw the cutpoint, and current eigenvector 
        (plot/quiver (first mean)  (second mean)  (first eigenvector) (second eigenvector) :color "orange" :scale 2)  ;; current eigenvector
        ;; (plot/scatter (first cluster2X) (first cluster2Y) :c "yellow") ;; cut-point

        ;; Draw the hyperplane going thru the cutpoint
        (plot/quiver (first cluster2X)  (first cluster2Y) (second eigenvector) (- (first eigenvector)) :color "red" :scale 2)

        ;; Set legend and title
        (plot/legend ["first cluster" "second cluster" "eigenvec" "hyperplane"])

        (plot/title "Cut Result")
        (plot/tight_layout)
        (plot/savefig (str "./src/pcd/picture-output/Chop/" pic-name) :bbox_inches "tight")))

    ;;;;; RETRUN DATA
    (if (= cut-point 0)
      ;; if cut-point is still = 0, simply return the input data
      [(vec (map first data-and-dist-vector))] ;; return an [] array of data, no distances 

      ;; else, there is a cut at result-cut-point, return 2 clusters of data
      [cluster-data-1 cluster-data-2])))

(defn chop-cluster
  "find subclusters, given an array of sorted eigenvectors and a dataset"
  [data is-plot]
    ;; 2. Try to chop each sub-cluster inside the unchopped-clusters
  (loop [unchopped-clusters [data] ;; All clusters in this array still need to be processed  
         final-chopped-clusters []
         cut-time 0
         pic-name (str "iteration-" cut-time)];; All clusters in this array are in their best shape and cannot be chopped down anymore

      ;; Debugging println
    (println "")
    (println "-----------")
    (println "final-chopped-clusters:" (count final-chopped-clusters))

      ;; if no more clusters in the unchopped-clusters
    (if (== (count unchopped-clusters) 0)
        ;; if true, return the final-chopped-list
      final-chopped-clusters

        ;; else
      (let [curr-cluster (first unchopped-clusters)
              ;; Get the sorted eigenvectors for the current cluster
            sorted-eigen-vectors (get-sorted-eigen-vectors curr-cluster)

              ;; for each eigen-vector, try to chop the curr-cluster. if chopped, break
            chopped-clusters   (loop
                                [idx-eigvec   0
                                 chopped      false
                                 after-chopped-clusters []]

                                   ;; Loop break-out condition
                                 (if (or (>= idx-eigvec (count sorted-eigen-vectors)))
                                         (= chopped true))

                                     ;; if done, return the after-chopped-clusters
                                   after-chopped-clusters

                                     ;; else, try to chop
                                   (let [curr-eivector (nth sorted-eigen-vectors idx-eigvec)
                                         sorted-dist (get-descending-sort-dist curr-eivector curr-cluster is-plot sorted-eigen-vectors)
                                          ;; sub-curr-clusters should be an array of either 1 or 2 clusters. No more or less        
                                         sub-curr-clusters (maybe-divide-cluster-using-DL sorted-dist curr-eivector is-plot pic-name)]

                                     (println "eigenvectors:" sorted-eigen-vectors)
                                     (println "trying to chop using: " curr-eivector)
                                     (recur (inc idx-eigvec)
                                            (if (== (count sub-curr-clusters) 1) false true)
                                            sub-curr-clusters)))) ;; chopped

            is-chopped (if (== (count chopped-clusters) 2) true false)]
        (println "chopped: " is-chopped)

        (if is-chopped
            ;; there is a chopped: (1) remove the processed cluster, (2) add the 2 sub clusters to the unchopped-clusters
          (recur (into (drop 1 unchopped-clusters) chopped-clusters)
                 final-chopped-clusters
                 (inc cut-time)
                 (str "iteration-" (inc cut-time)))
            ;; all set, no more chop: (1) remove the processed cluster, (2) add that removed cluster to final-chopped-clusters
          (recur (drop 1 unchopped-clusters)
                 (into final-chopped-clusters chopped-clusters)
                 (inc cut-time)
                 (str "iteration-" (inc cut-time))))))))

;; 5. Merging process
(defn mahalanobis-distance
  [point mean invcovar]
  ;; sqrt[(point - mean).T * invcovar * (point - mean)]

  (matrix/mget (matrix/inner-product
                (matrix/sub point mean)
                invcovar
                (matrix/transpose (matrix/sub point mean)))))

(defn merge-2-clusters-using-DL
  "Return true if we should merge 2 clusters. False otherwise."
  [cluster1 cluster2 is-plot picname]
  (let [;; Merged
        merged-cluster (vec (concat cluster1 cluster2)) ;; merge 2 clusters
        merged-sorted-vectors (get-sorted-eigen-vectors merged-cluster)
        merged-mean (matrix-stats/mean merged-cluster)

        ;; 2 clusters
        c1-sorted-vectors (get-sorted-eigen-vectors cluster1)
        c2-sorted-vectors (get-sorted-eigen-vectors cluster2)
        c1-mean (matrix-stats/mean cluster1)
        c2-mean (matrix-stats/mean cluster2)

        ;; assignment cost
        c1-length (count cluster1)
        c2-length (count cluster2)
        cac 0 #_(cluster-assignment-cost c1-length c2-length)
        
        max-iter (- (count merged-sorted-vectors) 1)] ;; cost of asignment cost
        
    (loop
     [eigen-index 0
      return-val false]
    
      ;; Loop Break-out condition
      (if (or
           (>= eigen-index max-iter)
           (= return-val true))

        ;; Return value
        return-val

        ;; Loop
        (let [;; get distances from points in merged cluster to hyperplane of merged cluster
              curr-evector-m (nth merged-sorted-vectors eigen-index)
              hyperplane-merged (get-hyperplane-equation merged-mean curr-evector-m)
              dist-from-plane-to-points-merged (get-dist-from-plane-to-points-V2 hyperplane-merged merged-cluster)

              ;; get distances from points in cluster 1 to hyperplane of cluster 1
              curr-evector-c1 (nth c1-sorted-vectors eigen-index)
              hyperplane-c1 (get-hyperplane-equation c1-mean curr-evector-c1)
              dist-from-plane-to-points-c1 (get-dist-from-plane-to-points-V2 hyperplane-c1 cluster1)

              ;; get distances from points in cluster 2 to hyperplane of cluster 2
              curr-evector-c2 (nth c2-sorted-vectors eigen-index)
              hyperplane-c2 (get-hyperplane-equation c2-mean curr-evector-c2)
              dist-from-plane-to-points-c2 (get-dist-from-plane-to-points-V2 hyperplane-c2 cluster2)

              dl-1 (description-length dist-from-plane-to-points-c1)
              dl-2 (description-length dist-from-plane-to-points-c2)
              dl-1plus2 (+ dl-1 dl-2 cac)
              dl-merged (description-length dist-from-plane-to-points-merged)
              
              merged? (< dl-merged dl-1plus2)]

          ;; VISUALIZATION
          ;; Draw the 2 input clusters with plot title: MERGED? TRUE or FALSE
          ;; Show dl-1plus2 vs dl-merged
          (when is-plot
            (plot/subplots :figsize [15 5])

            ;; Draw merged cluster, with sorted eigen-vec
            (plot/subplot 1 3 1 :adjustable "datalim" :aspect 1)
            (plot/scatter  (into [] (map first cluster1)) (into [] (map last cluster1)) :c "pink")
            (plot/scatter  (into [] (map first cluster2)) (into [] (map last cluster2)) :c "blue")
            (plot/quiver (first merged-mean)  (second merged-mean)  (first (first merged-sorted-vectors)) (second  (first merged-sorted-vectors)) :color "orange" :scale 2)
            (plot/quiver (first merged-mean)  (second merged-mean) (first (second merged-sorted-vectors)) (second  (second merged-sorted-vectors)) :color "red" :scale 2)
            (plot/legend ["cluster1" "cluster2" "1st evec" "2nd evec"])
            (plot/title "Merged Cluster")

            ;; Draw merged cluster, with current-eigen and hyperplane
            (plot/subplot 1 3 2 :adjustable "datalim" :aspect 1)
            (plot/scatter  (into [] (map first merged-cluster)) (into [] (map last merged-cluster)) :c "black")
            (plot/quiver (first merged-mean)  (second merged-mean)  (first curr-evector-m) (second curr-evector-m) :color "orange" :scale 2)
            (plot/quiver (first merged-mean)  (second merged-mean) (second curr-evector-m) (- (first curr-evector-m)) :color "red" :scale 2)
            (plot/legend ["data" "cur evec" "plane"])
            (plot/title "Merged Cluster")

             ;; Text Info 
            (plot/subplot 1 3 3 :adjustable "datalim" :aspect 1)
            (plot/text 0.5 0.1 (str "Merged: " merged?) :ha "center" :va "center" :size 10)
            (plot/text 0.5 0.5 (str "dl-1 + dl2 = " dl-1plus2) :ha "center" :va "center" :size 10)
            (plot/text 0.5 0.7 (str "dl-merged = " dl-merged) :ha "center" :va "center" :size 10)

            ;; Save figure
            (plot/savefig (str "./src/pcd/picture-output/Merge/" picname)))
          
          (if merged?
            (recur (inc eigen-index) true)
            (recur (inc eigen-index) false)))))))

(defn sort-by-distance
  [cluster clusters]
  (let [cluster-mean (matrix-stats/mean cluster)]
    (sort-by #(if (<= (count %) 2)
                (euclidean-distance cluster-mean (matrix-stats/mean %))
                (mahalanobis-distance cluster-mean
                                      (matrix-stats/mean %)
                                      (matrix/inverse (in-stats/covariance %))))
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
  [clusters is-plot]

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
                               (if (= true (merge-2-clusters-using-DL first-cluster sorted-cluster is-plot (str "Merge-" length)))   
                                 ;; If a merge happens, return the new vector of clusters, where the merged clusters is added, and the 2 clusters are removed  
                                 (reduced
                                    (concat
                                     (vec-remove sorted-rest-clusters @sorted-cluster-idx) ;; remove the sorted-cluster from sorted-rest-clusters    
                                     [(vec (concat first-cluster sorted-cluster))]         ;; the merged cluster
                                     (subvec sorted-clusters 0 index)))                  ;; all clusters iterated (not included in sorted-rest-clusters)   
                                 ;; Else, do nothing, the last iteration will return first-cluster, which should be "current-cluster"
                                 (do
                                   (swap! sorted-cluster-idx inc)
                                   first-cluster))) joined-clusters)]

          (println "-------------")
          (if (= output current-cluster)
              ;; No merge
            (do
              (println "NO MERGED" (count sorted-clusters))
              (recur (inc index) ;; move onto the next cluster
                     clusters))  ;; result = original clusters since no merge
              ;; A merge happened
            (do
              (println "MERGED")
              (merge-clusters (vec output) is-plot)))) ;; recursively call merge-clusters with the new vector
        result))))

(defn -main
  "I don't do a whole lot ... yet."
  [& args]

  ;; Need this to handle "NSWindow drag regions should only be invalidated on the Main Thread!"
  (let [mplot (py/import-module "matplotlib")]
    (py/call-attr mplot "use" "WebAgg")
    (py/call-attr mplot "use" "Agg"))

  ;; The PCD process
  (let [;; 1. Process data
        ds (process-data infile-path (list 0 1 2)) ;; import data and remove columns 0 1 2
        ds-float (matrix-float (count ds) ds)      ;; convert matrix of string to float
        is-plot true

        ;; 3. Chop the dataset
        chopped-clusters (chop-cluster ds-float false)
        
        ;; 4. Merge the clusters from (3)
        merged-clusters (merge-clusters chopped-clusters is-plot)]

    ;; 5. OPTIONAL: VISUALIZATION
    (when is-plot
      (plot/subplots :figsize [10 5])
      ;; 1. Draw chopped-clusters
      (plot/subplot 1 2 1 :adjustable "box" :aspect "equal")
      (draw-final chopped-clusters)


       ;; 2. Draw chopped-clusters
      (plot/subplot 1 2 2 :adjustable "box" :aspect "equal")
      (draw-final merged-clusters)

      ;; Save image
      (plot/tight_layout)
      (plot/savefig "./src/pcd/picture-output/Final-Result"))))