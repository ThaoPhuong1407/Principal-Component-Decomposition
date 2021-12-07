;;; -*- Syntax: Yolambda; -*-

;;; Matrix manipulation package

(define zero (sffloat 0 (makeIeeeSingleFloatResult)))
(define one (sffloat 1 (makeIeeeSingleFloatResult)))
(define inverse (lambda (x) (/ one x)))

;(define sum
;  (lambda (n fn)
;    (define (sumHelper i partialSum)
;      (if (zero? i)
;	  (+ (fn 0) partialSum)
;	  (sumHelper (- i 1) (+ (fn i) partialSum))))
;    (sumHelper (- n 1) zero)))

;(define sum
;  (lambda (n fn)
;    (let ((partialSum zero))
;      (do ((i (- n 1) (- i 1)))
;	  ((zero? i) (+ partialsum (fn 0)))
;	  (set! partialSum (+ partialSum (fn i)))))))

(defineClass Matrix ()
  ((dims)
   (data initform nil)))

(define (print Matrix|object port slash? . options)
  (withSlots (dims data) object
    (let ((dim (car dims)))
      (format port "{Matrix[")
      (dotimes (i dim)
        (format port "(")
        (let ((row (vectorRef data i)))
          (dotimes (j (cadr dims))
            (format port " ~a " (vectorRef row j)))
          (format port ")")))
      (format port "]}"))))

(define makeZeroMatrix
  (lambda (dim . dim2)
    (if (null? dim2) (set! dim2 dim) (set! dim2 (car dim2)))
    (let ((mat (make Matrix dims: (list dim dim2) data: (makeVector dim))))
      (withSlots (data) mat
        (dotimes (i dim)
          (vectorSet! data i (makeVector dim2 zero)))
        mat))))

(define makeMatrix
  (lambda (dim . dim2)
    (if (null? dim2) (set! dim2 dim) (set! dim2 (car dim2)))
    (makeZeroMatrix dim dim2)))

(define (matrixRef Matrix|mat i j)
  (withSlots (data) mat
    (vectorRef (vectorRef data i) j)))

(define (matrixSet! Matrix|mat i j value)
  (withSlots (data) mat
    (vectorSet! (vectorRef data i) j value)
    mat))

(define (matrixInc! Matrix|mat i j inc)
  (matrixSet! mat i j (+ (matrixRef mat i j) inc)))

(define (matrixDec! Matrix|mat i j inc)
  (matrixSet! mat i j (- (matrixRef mat i j) inc)))

(define (matrixMult! Matrix|mat i j inc)
  (matrixSet! mat i j (* (matrixRef mat i j) inc)))

(define (matrixDiv! Matrix|mat i j inc)
  (matrixSet! mat i j (/ (matrixRef mat i j)  inc)))

(define (matrixDimension Matrix|mat . n)
  (if (null? n) (set! n 0) (set! n (car n)))
  (withSlots (dims data) mat
    (if (= n 0) (car dims) (cadr dims))))

(define makeUnitMatrix
  (lambda (dim)
    (let ((mat (makeZeroMatrix dim)))
      (dotimes (i dim)
        (matrixSet! mat i i one))
      mat)))

;(define (matrixSwapElts Martix|mat i1 j1 i2 j2)
;  (withSlots (data) mat
;    (let ((temp (matrixRef data i1 j1)))
;      (matrixSet! data i1 j1 (matrixRef data i2 j2))
;      (matrixSet! data i2 j2 temp))))

(define (matrixSwapRows! Matrix|mat row1 row2)
  (withSlots (data) mat
    (let ((temp (vectorRef data row1)))
      (vectorSet! data row1 (vectorRef data row2))
      (vectorSet! data row2 temp))))

(define (copyMatrix Matrix|mat)
  (let* ((n (matrixDimension mat))
	 (m (matrixDimension mat 1))
         (copy (makeMatrix n m)))
    (dotimes (i n)
      (dotimes (j m)
        (matrixSet! copy i j (matrixRef mat i j))))
    copy))

(define (matrixTranspose Matrix|mat)
  (withSlots (dims) mat
    (let* ((n (matrixDimension mat))
	   (m (matrixDimension mat 1))
	   (copy (makeMatrix m n)))
      (dotimes (i n)
	(dotimes (j m)
	  (matrixSet! copy j i (matrixRef mat i j))))
      copy)))

(define (rowVectorMatrixMultiply row Matrix|mat)
  (let* ((dim (matrixDimension mat))
         (result (makeVector dim nil)))
    (if (not (= dim (vectorLength row)))
        (error "Can't multiply a row (~a) by a ~a by ~a matrix." (vectorLength row) dim dim))
    (dotimes (i dim result)
      (let ((sum zero))
        (dotimes (j dim)
          (set! sum (+ sum (* (vectorRef row j) (matrixRef mat i j)))))
        (vectorSet! result i sum)))))

(define vectorMult
  (lambda (row col)
    (let* ((dim (vectorLength row))
           (result zero))
      (if (not (= dim (vectorLength col)))
          (error "Can't multiply a row (~a) by a column (~a) vector." dim (vectorLength col)))
      (dotimes (i dim result)
        (set! result (+ result (* (vectorRef row i) (vectorRef col i))))))))

;;; matrix multiplication

;(define (matrixMultiply Matrix|mat1 mat2)
;  (let* ((n (matrixDimension mat1))
;         (ans (makeZeroMatrix n)))
;    (dotimes (i n)
;      (dotimes (j n)
;        (matrixSet!
;          ans i j
;          (sum n (lambda (k)
;                   (* (matrixRef mat1 i k)
;                      (matrixRef mat2 k j)))))))
;    ans))

(define (matrixMultiply Matrix|mat1 mat2)
  (let* ((n (matrixDimension mat1))
	 (m (matrixDimension mat1 1))
	 (x (matrixDimension mat2))
	 (y (matrixDimension mat2 1))
	 (ans (makeZeroMatrix n y)))
    (if (not (= m x))
	(error "Can't multiply a ~a by ~a matrix and a ~a by ~a matrix." n m x y))
    (dotimes (i n)
      (dotimes (j y)
        (let ((sum zero))
          (dotimes (k m)
            (set! sum (+ sum (* (matrixRef mat1 i k)
                                (matrixRef mat2 k j)))))
          (matrixSet! ans i j sum))))
    ans))

(define (matrixSubtract Matrix|mat1 mat2)
  (let* ((n (matrixDimension mat1))
	 (m (matrixDimension mat1 1))
         (ans (makeMatrix n m)))
    (dotimes (i n)
      (dotimes (j m)
        (matrixSet! ans i j (- (matrixRef mat1 i j)
                               (matrixRef mat2 i j)))))
    ans))

(define (matrixAdd Matrix|mat1 mat2)
  (let* ((n (matrixDimension mat1))
	 (m (matrixDimension mat1 1))
         (ans (makeMatrix n m)))
    (dotimes (i n)
      (dotimes (j m)
        (matrixSet! ans i j (+ (matrixRef mat1 i j)
                               (matrixRef mat2 i j)))))
    ans))

;;; gauss-jordan does matrix inversion.  Note that this destructively changes 
;;; the argument matrix mat

(define (gaussJordan Matrix|mat)
  (let* ((n (matrixDimension mat))
         (inv (makeUnitMatrix n)))
    (dotimes (i n)
      (do ((pi i (+ 1 pi)))
          ((or (= pi n)
               (not (zero? (matrixRef mat pi i))))
           (cond ((= pi n) (error "GAUSS-JORDAN: matrix to invert is singular"))
                 ((not (= i pi)) (matrixSwapRows! mat i pi)
                  (matrixSwapRows! inv i pi)))))
      (let ((pivot (matrixRef mat i i)))
        (dotimes (j n)
          (if (not (= j i))
              (let ((factor (/ (matrixRef mat j i) pivot)))
                (matrixSet! mat j i zero)
                (do ((k (+ 1 i) (+ 1 k)))
                    ((= k n))
                    (matrixDec! mat j k (* (matrixRef mat i k) factor)))
                (dotimes (k n)
                  (matrixDec! inv j k (* (matrixRef inv i k) factor))))))
        (matrixSet! mat i i one)
        (do ((j (+ 1 i) (+ 1 j)))
            ((= j n))
            (matrixDiv! mat i j pivot))
        (dotimes (j n)
          (matrixDiv! inv i j pivot))))
    inv))

(define (matrixTrace Matrix|mat)
  (withSlots (data) mat
    (let ((rows (matrixDimension mat))
          (cols (vectorLength (vectorRef data 0)))
          (result zero))
      (if (= rows cols)
          (dotimes (i rows result)
            (set! result (+ result (matrixRef mat i i))))
          (error "Can't compute trace of a non-square matrix [~a,~a]." cols rows)))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                               MATLAB Interface                   ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define mxREAL 0)
(define mxCOMPLEX 1)

(define *matlabhandle* #f)
(define *matlabserver* "")                ; local server
(define *lastmatlabcommand* #f)
(define *tempvar* (symbol->string (gensym "ML")))

(defineClass MatlabMatrix (Matrix)
  ((data initform nil)
   (rows initform 1)
   (cols initform 1)
   (name initform #f)))

(define openMatlab
  (lambda ()
    (if (not *matlabhandle*)
        (set! *matlabhandle* (engOpen *matlabserver*)))))

(define closeMatlab
  (lambda ()
    (if *matlabhandle* (engClose *matlabhandle*))
    (set! *matlabhandle* #f)))

(define restartMatlab
  (lambda ()
    (closeMatlab)
    (openMatlab)))

(define (print MatlabMatrix|object port slash? . options)
  (withSlots (data rows cols name) object
    (if (null? data)
        (format "{Matrix destroyed}")
        (begin
          (format port "{Matrix ~a [" name)
          (dotimes (i rows)
            (format port "(")
            (dotimes (j cols)
              (format port " ~a " (matlabGetArrayValue data (+ (* rows j) i) (makeIeeeSingleFloatResult))))
            (format port ")"))
          (format port "]}")))))

(define makeZeroMatlabMatrix
  (lambda (dim . rest)
    (let ((mat (make MatlabMatrix
                     data: (mxCreateDoubleMatrix dim (if (pair? rest) (car rest) dim) mxREAL)
                     rows: dim
                     cols: (if (pair? rest) (car rest) dim)
                     name: (symbol->string (gensym "ML")))))
      (withSlots (data name) mat
        (openMatlab)                        ; Make sure that Matlab has been started
        (mxSetName data name)
        (engPutArray *matlabhandle* data))
      mat)))

(define makeMatlabMatrix
  (lambda (name dim . rest)
    (apply makeZeroMatlabMatrix dim rest)))

(define makeNamedMatlabMatrix
  (lambda (name dim . rest)
    (let ((mat (make MatlabMatrix
                     data: (mxCreateDoubleMatrix dim (if (pair? rest) (car rest) dim) mxREAL)
                     rows: dim
                     cols: (if (pair? rest) (car rest) dim)
                     name: name)))
      (withSlots (data name) mat
        (openMatlab)                        ; Make sure that Matlab has been started
        (mxSetName data name)
        (engPutArray *matlabhandle* data))
      mat)))

(define (matrixRef MatlabMatrix|mat i j)
  (withSlots (data rows) mat
    (matlabGetArrayValue data (+ (* rows j) i) (makeIeeeSingleFloatResult))))

(define (matrixSet! MatlabMatrix|mat i j value)
  (withSlots (data rows) mat
    (matlabSetArrayValue data (+ (* rows j) i) (+ zero value))
    (engPutArray *matlabhandle* data)
    mat))

(define (matrixDimension MatlabMatrix|mat)
  (withSlots (rows) mat
    rows))

(define makeUnitMatlabMatrix
  (lambda (dim)
    (let ((mat (makeZeroMatlabMatrix dim)))
      (dotimes (i dim)
        (matrixSet! mat i i one))
      mat)))

(define (matrixGetM MatlabMatrix|mat)
  (withSlots (data) mat
    (mxGetM data)))

(define (matrixGetN MatlabMatrix|mat)
  (withSlots (data) mat
    (mxGetN data)))

(define (copyMatrix MatlabMatrix|mat)
  (let* ((n (matrixDimension mat))
         (copy (makeMatlabMatrix n)))
    (dotimes (i n)
      (dotimes (j n)
        (matrixSet! copy i j (matrixRef mat i j))))
    copy))

(define (matrixMultiply MatlabMatrix|mat1 mat2)
  (let* ((n (matrixDimension mat1))
         (ans (makeZeroMatlabMatrix n)))
    (dotimes (i n)
      (dotimes (j n)
        (let ((sum zero))
          (dotimes (k n)
            (set! sum (+ sum (* (matrixRef mat1 i k)
                                (matrixRef mat2 k j)))))
          (matrixSet! ans i j sum))))
    ans))

(define (matrixSubtract MatlabMatrix|mat1 mat2)
  (let* ((n (matrixDimension mat1))
         (ans (makeMatlabMatrix n)))
    (dotimes (i n)
      (dotimes (j n)
        (matrixSet! ans i j (- (matrixRef mat1 i j)
                               (matrixRef mat2 i j)))))
    ans))

(define (matrixTrace MatlabMatrix|mat)
  (withSlots (data rows cols) mat
    (let ((result zero))
      (if (= rows cols)
          (dotimes (i rows result)
            (set! result (+ result (matrixRef mat i i))))
          (error "Can't compute trace of a non-square matrix [~a,~a]." cols rows)))))

(define (matrixTranspose MatlabMatrix|mat)
  (withSlots (name) mat
    (let ((res (symbol->string (gensym "ML"))))
      (mlEval (format #f "~a=transpose(~a)" res name))
      (getMatlabMatrix res))))

(define getMatlabMatrix
  (lambda (name)
    (let* ((mat (make MatlabMatrix name: name)))
      (withSlots (data rows cols name) mat
        (set! data (engGetArray *matlabhandle* name))
        (set! rows (mxGetM data))
        (set! cols (mxGetN data))
      mat))))

(define (putMatlabMatrix MatlabMatrix|mat)
  (withSlots (data) mat
    (engPutArray *matlabhandle* data)))

(define getMatlabDouble
  (lambda (name)
    (let* ((data (engGetArray *matlabhandle* name))
           (res (matlabGetArrayValue data 0 (makeIeeeSingleFloatResult))))
      (mxDestroyArray data)
      res)))

(define (matrixRefresh MatlabMatrix|mat)
  (withSlots (data rows cols name) mat
    (mxDestroyArray data)
    (set! data (engGetArray *matlabhandle* name))
    (set! rows (mxGetM data))
    (set! cols (mxGetN data)))
  mat)

(define (matrixDestroy MatlabMatrix|mat)
  (withSlots (data rows cols name) mat
    (mxDestroyArray data)
    (set! data nil)
    (set! rows nil)
    (set! cols nil)
    (set! name nil))
  mat)    

(define (matrixName MatlabMatrix|m)
  (withSlots (name) m name))

(define (mahalanobisDistance MatlabMatrix|x mean invcovar)
  (mlEval (format #f "~a-~a;sqrt((ans*~a)*transpose(ans));"
                  (matrixName x) (matrixName mean) (matrixName invcovar)))
  (getMatlabDouble "ans"))

(define (euclideanDistance  MatlabMatrix|x mean invcovar)
  (let ((sumsq zero)
        (dims (matrixGetN x)))
    (dotimes (i dims)
      (let ((diff (- (matrixRef x 0 i) (matrixRef mean 0 i))))
        (set! sumsq (+ sumsq (* diff diff)))))
    (sqrt sumsq)))

(define (bigDistance MatlabMatrix|x mean invcovar)
  (let ((sumsq zero)
        (dims (matrixGetN x)))
    (dotimes (i dims)
      (let ((diff (- (matrixRef x 0 i) (matrixRef mean 0 i))))
        (set! sumsq (+ sumsq (* diff diff)))))
    (* (sqrt sumsq) 1000))) ;+++ arbitrary number
    
(define mlEval
  (lambda (command)
    (set! *lastmatlabcommand* command)
    (openMatlab)                        ; Make sure that Matlab has been started
    (engEvalString *matlabhandle* command)))
    
(define (matrixInverse MatlabMatrix|mat)
  (withSlots (name) mat
    (let ((resname (symbol->string (gensym "ML"))))
      (mlEval (format #f "~a=inv(~a)" resname name))
      (getMatlabMatrix resname))))

(define (matrixEigenvalues MatlabMatrix|mat)
  (withSlots (name) mat
    (let ((resname (symbol->string (gensym "ML"))))
      (mlEval (format #f "~a=eig(~a)" resname name))
      (getMatlabMatrix resname))))

(define (matrixEigenvectorsAndEigenvalues MatlabMatrix|mat)
  (withSlots (name) mat
    (let ((valname (symbol->string (gensym "ML")))
          (vecname (symbol->string (gensym "ML"))))
      (mlEval (format #f "[~a,~a]=eig(~a)" vecname valname name))
      (list (getMatlabMatrix vecname) (getMatlabMatrix valname)))))

(define (matrixDeterminant MatlabMatrix|mat)
  (withSlots (name) mat
    (mlEval (format #f "~a=det(~a)" *tempvar* name))
    (matlabGetArrayValue
      (engGetArray *matlabhandle* *tempvar*)
      0
      (makeIeeeSingleFloatResult))))

;;; converts a yolambda vector into a Matlab column vector
(define vector->matlabMatrix
  (lambda (vec)
    (let ((mlv (makeZeroMatlabMatrix (size vec) 1)))
      (dotimes (i (size vec))
        (matrixSet! mlv i 0 (vectorRef vec i)))
      mlv)))

;;; Converts a Yolambda list into a Matlab column vector
(define list->matlabMatrix
  (lambda (lst)
    (let* ((ll (length lst))
           (mlv (makeZeroMatlabMatrix ll 1)))
      (dotimes (i ll)
        (matrixSet! mlv i 0 (car lst))
        (set! lst (cdr lst)))
      mlv)))

(define (asMatlabMatrix Vector|vec)
  (vector->matlabMatrix vec))

(define (asMatlabMatrix Pair|lst)
  (list->matlabMatrix lst))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; PCA decomposition and categorizing


;;; Computes the eigenvalues and eigenvectors of the provided covarince matrix
;;; and sorts them into decreasing order of eigenvalue.

(define computeSortedEigenvectorsFromCovarianceMatrix
  (lambda (covarmat)
    (destructuringBind (evecs evals) (matrixEigenvectorsAndEigenvalues covarmat)
       ;;; Now sort them
       (withSlots (rows cols) evecs
         (dotimes (i cols (list evecs evals))
           (do ((j (+ i 1) (+ j 1)))
               ((= j cols) #f)
               (let ((ival (matrixRef evals i i))
                     (jval (matrixRef evals j j)))
                  (when (> jval ival)
                     (dotimes (s rows)
                       (let ((tmp (matrixRef evecs s j)))
                         (matrixSet! evecs s j (matrixRef evecs s i))
                         (matrixSet! evecs s i tmp)))
                     (matrixSet! evals j j ival)
                     (matrixSet! evals i i jval)))))))))

(define meanAndCovariance
  (lambda (data) ; Data is a list of vectors one for each dimension.
    (let ((cdv (symbol->string (gensym "ML")))
          (means (symbol->string (gensym "ML")))
          (cvars (symbol->string (gensym "ML")))
          (points (matrixGetM (car data)))
          (dims (length data))
          (dmat "["))
      (cond ((= points 1)
             ;; We have only a single data point so the mean withh be that point
             ;; and the covariance matrix will be zero.
             (let ((means (makeZeroMatlabMatrix 1 dims))
                   (covars (makeZeroMatlabMatrix dims dims)))
              (do ((dl data (cdr dl))
                   (ind 0 (+ ind 1)))
                  ((null? dl) (list means covars))
                  (withSlots (name) (car dl)
                    (matrixSet! means 0 ind (matrixRef (car dl) 0 0))))))
            (else
              (do ((dl data (cdr dl)))
                  ((null? dl) (set! dmat (stringAppend dmat "]")))
                  (withSlots (name) (car dl)
                    (set! dmat (stringAppend dmat name))
                    (if (not (null? (cdr dl))) (set! dmat (stringAppend dmat ",")))))
              (mlEval (format #f "~a=~a;" cdv dmat))
              (mlEval (format #f "~a=mean(~a);" means cdv))
              (mlEval (format #f "~a=cov(~a);" cvars cdv))
              (list (getMatlabMatrix means) (getMatlabMatrix cvars)))))))

;;; Compute the equation of the hyperplane that is the selected eigenvector
;;; as a normal form passing through the (supplied) mean.

(define hyperplaneEquationPerpendicularForm
  (lambda (means eigs eignum)
    (let* ((evecs (car eigs))
           (normname (symbol->string (gensym "ML"))))
      (withSlots (name) evecs
        (mlEval (format #f "~a=~a(:,~a).';" normname name (+ eignum 1))))
      (let ((norm (getMatlabMatrix normname))
            (n 0))
        (dotimes (i (matrixGetN means))
          (set! n (+ n (* (matrixRef means i 0) (matrixRef norm i 0)))))
        (if (< n 0)
            (dotimes (i (matrixGetN means) (set! n (- n)))
              (matrixSet! norm i 0 (- (matrixRef norm i 0)))))
        (list norm n)))))

;;; Compute the distance of a point from a hyperplane.  The hyperplane is in the
;;; form of a list (norm n) being the equation of the hyperplane "r.norm=n".
;;; The point is a row matrix.

;;; (define pointDistanceFromHyperplane
;;;  (lambda (point hyperplane)
;;;    (destructuringBind (norm n) hyperplane
;;;      (let ((dist n))
;;;        (dotimes (i (matrixGetN norm) dist)
;;;          (set! dist (- dist (* (matrixRef point i 0) (matrixRef norm i 0)))))))))

;;; Compute the distance of a point from a hyperplane.  The hyperplane is in the
;;; form of a list (norm n) being the equation of the hyperplane "r.norm=n".
;;; The point is a list if vectors of data (data) and an index (point).

(define pointDistanceFromHyperplaneList
  (lambda (data point hyperplane)
    (destructuringBind (norm n) hyperplane
      (let ((dist n))
        (dotimes (i (matrixGetN norm) dist)
          (set! dist (- dist (* (matrixRef (at data i) point 0) (matrixRef norm i 0)))))))))

(define log2 (lambda (n) (/ (log n) (log 2))))

(define computeDistanceList
  (lambda (data n hyperplane)
    (let ((dists (makeZeroMatlabMatrix n 1)))
      (dotimes (i n)
        (matrixSet! dists i 0 (pointDistanceFromHyperplaneList data i hyperplane)))
      dists)))

;;; This one is the winner!  It's fast and accurate.
(define clusterDescriptionLength
  (lambda (dists from to mean variance)
    (let* ((points (- to from)))
      (+ 20
        (* points (if (= variance 0) 0 (log2 (sqrt variance))))))))

;;; This computes the same entropy 'point at a time' as the above version.
;;; It is much too slow to use in practice but it illustrates how entropy
;;; is defined.  Compared to the above version there is a constant addition
;;; governed by the arbitrary resolution (0.01).  Since we are only interested
;;; in relative changes in entropy they are both identical (except in performance).

(define clusterDescriptionLength-slow
  (lambda (dists from to mean variance)
    (if (= variance 0)
        0
        (let* ((points (- to from))
               (work (makeZeroMatlabMatrix points 1))
               (sd (sqrt variance))
               (wname (withSlots (name) work name)))
          (dotimes (i points)
            (matrixSet! work i 0 (- (matrixRef dists (+ i from) 0) mean)))
          (mlEval (format #f "ans=erf((~a+0.01)/~a)-erf((~a-0.01)/~a)" wname sd wname sd))
          (mlEval (format #f "ones(1,~a)*(-log2(ans))" points))
          (matrixDestroy work)
          (getMatlabDouble "ans")))))

(define clusterAssignmentCost
  (lambda (pc1 pc2)
    (if (or (zero? pc1) (zero? pc2))
        zero
        (let* ((pc0 (+ pc1 pc2 zero))
               (pr1 (/ pc1 pc0))
               (pr2 (/ pc2 pc0)))
          (- (* pc0 (+ (* pr1 (log2 pr1)) (* pr2 (log2 pr2)))))))))

(define clusterAssignmentCostGlobal
  (lambda (pc1 pc2 tot)
    (set! tot (+ tot zero))
    (if (or (zero? pc1) (zero? pc2))
        zero
        (let* ((pc0 (+ pc1 pc2 zero))
               (pr0 (/ pc0 tot))
               (pr1 (/ pc1 tot))
               (pr2 (/ pc2 tot)))
          (- (- (* pc0 (+ (* pr1 (log2 pr1)) (* pr2 (log2 pr2)))))
             (- (* pc0 (+ (* pr0 (log2 pr0))))))))))

;;; Look for an opportunity to divide the population in the direction
;;; that is orthogonal to the hyperplane represented by one of the eigenvectors.
;;; Eigen vectors are checked in (decreasing) order of related eigenvalue (by virtue
;;; of the eigenvectors being sorted beforehand). If the population is split
;;; a list of the two datasets is returned as a result.  If the population
;;; was not split #f is returned.
;;; variance is computed incrementally as follows:
;;;
;;; n*mean^2 + Sigma<i=1,n> x(i)^2 -2*mean*Sigma<i=1,n> x(i)
;;; --------------------------------------------------------
;;;                            n-1

(define maybeDividePopulation                          ; +++ make cleanup of matrices vork
  (lambda (data n objs mv eigs numpoints plot?)
    (let* ((md (makeVector n 0))
           (vd (makeVector n 0))
           (ed (makeVector n 0))
           (mu (makeVector n 0))
           (vu (makeVector n 0))
           (eu (makeVector n 0))
           (etot (makeZeroMatlabMatrix (+ n 1) 1))
           (best 0)
           (bestindex 0)
           (result #f)
           (sorted #f)
           (order #f)
           (sortednm (symbol->string (gensym "ML")))
           (ordernm (symbol->string (gensym "ML"))))
      (do ((eignum 0 (+ 1 eignum)))
          ((or (>= eignum (matrixGetN (car eigs))) result)
           (when plot?
             (mlEval "subplot(2,1,2);")
             (mlEval (format #f "plot(~a);" (withSlots (name) *etot* name)))
             (mlEval "title('Change in description length vs cluster division point');")
             (mlEval "xlabel('Points in order along normal to the hyperplane','fontsize',8);")
             (mlEval "ylabel('Relative Entropy','fontsize',8);")
             (mlEval "hold on;")
             (mlEval "legend('Change in entropy',-1);")
             (when result ; if we split
               (mlEval (format #f "plot([~a],[~a],'hr')" (+ bestindex 2) best))
               (mlEval "legend('Change in entropy','Minimum entropy cut point',-1);")))
           result)
         ;(format "Trying to divide along eigenvector #~a~%" eignum)
         (let* ((hyperplane (hyperplaneEquationPerpendicularForm (car mv) eigs eignum))
                (dists (computeDistanceList data n hyperplane)))
           (mlEval (format #f "[~a,~a]=sort(~a);" sortednm ordernm (withSlots (name) dists name)))
           (set! sorted (getMatlabMatrix sortednm))
           (set! order (getMatlabMatrix ordernm))
           ;; Compute mean, variances and entropy down
           (let ((sum 0) (sumsq 0))
             (dotimes (i n)
               (let ((xi (matrixRef sorted i 0)))
                 (set! sum (+ sum xi))
                 (let ((mi (/ sum (+ i 1))))
                   (vectorSet! md i mi)        ; Set the mean
                   (set! sumsq (+ sumsq (* xi xi)))
                   (let* ((var (if (> i 0) ; Set the variance
                                   (/
                                     (+ (* (+ i 1) (* mi mi)) ;n*mean^2
                                        sumsq                 ;sigma(xi^2)
                                        (- (* 2 mi sum)))     ;-2*mean*sigma(xi)
                                     i)                       ;n-1
                                   0))
                          (ent (clusterDescriptionLength sorted 0 (+ i 1) mi var)))
                     (vectorSet! vd i var)
                     (vectorSet! ed i ent)))))) ; entropy down
           ;; Compute mean, variances and entropy up
           (let ((sum 0) (sumsq 0))
             (dotimes (i n)
               (let ((xi (matrixRef sorted (- n (+ i 1)) 0)))
                 (set! sum (+ sum xi))
                 (let ((mi (/ sum (+ i 1))))
                   (vectorSet! mu (- n (+ i 1)) mi)        ; Set the mean
                   (set! sumsq (+ sumsq (* xi xi)))
                   (let* ((var (if (> i 0) ; Set the variance
                                   (/
                                     (+ (* (+ i 1) (* mi mi)) ;n*mean^2
                                        sumsq                 ;sigma(xi^2)
                                        (- (* 2 mi sum)))     ;-2*mean*sigma(xi)
                                     i)                       ;n-1
                                   0))
                          (ent (clusterDescriptionLength sorted (- n (+ i 1)) n mi var)))
                     (vectorSet! vu (- n (+ i 1)) var) ; variance up
                     (vectorSet! eu (- n (+ i 1)) ent)))))) ; entropy up
           (set! best one)
           (dotimes (i (+ n 1))
             (let* ((etotal (+        
                              (clusterAssignmentCost i (- n i) numpoints)
                              (if (< i n) (vectorRef eu i) zero)
                              (if (> i 0) (vectorRef ed (- i 1)) zero)
                              (- (vectorRef eu 0)))))
               (matrixSet! etot i 0 etotal)
               (when (< etotal best)
                 (set! bestindex (- i 1))
                 (set! best etotal))))
           (set! *md* md)
           (set! *vd* vd)
           (set! *ed* ed)
           (set! *mu* mu)
           (set! *vu* vu)
           (set! *eu* eu)
           (set! *etot* etot)
           (set! *order* order)
           (set! *sorted* sorted)
           ;; Cleanup work matrices
           ;(mlEval (format #f "clear ~a ~a;"
                            ;              (withSlots (name) dists name)
                            ;              (withSlots (name) sorted name)
                            ;              (withSlots (name) order name)))
           ;(destroyMatrix dists)
           ;(destroyMatrix sorted)
           ;(destroyMatrix order)
           ;(format "Best=~s Bestindex=~a~%" best bestindex)
           (if (and (>= bestindex 0) (< bestindex (- n 1)))
               (set! result (splitAt data n objs order bestindex eigs eignum plot?)))
           result))))) [[data, data], [data, data]]

;;; Check whether merging two populations will result in reduced entropy.
;;; If merging reduces entropy the resulting merged cluster is returned,
;;; otherwise #f is returned.

(define mergeClusters
  (lambda (pop1 pop2)
    (let* ((npop1 (car pop1))
           (npop2 (car pop2)))
      (cond ((zero? npop1) pop2)
            ((zero? npop2) pop1)
            (else
              (let* ((cpop (+ npop1 npop2))
                     (cobjs (makeVector cpop #f))
                     (cdata (map (lambda (x) (makeZeroMatlabMatrix cpop 1)) (cadr pop1)))
                     (pos 0))
                (dotimes (i npop1)
                  (vectorSet! cobjs i (vectorRef (caddr pop1) i)))
                (dotimes (i npop2)
                  (vectorSet! cobjs (+ npop1 i) (vectorRef (caddr pop2) i)))
                (do ((dims cdata (cdr dims))
                     (p1dims (cadr pop1) (cdr p1dims))
                     (p2dims (cadr pop2) (cdr p2dims)))
                    ((null? dims) (list cpop cdata cobjs))
                   ;; Populate from p1
                   (dotimes (i npop1)
                     (matrixSet! (car dims) i 0 (matrixRef (car p1dims) i 0)))
                   ;; Populate from p2
                   (dotimes (i npop2)
                     (matrixSet! (car dims) (+ i npop1) 0 (matrixRef (car p2dims) i 0))))))))))

(define varGivenMean
  (lambda (mn data)
    (let* ((n (matrixGetM data))
           (sumsq zero))
      (dotimes (i n)
        (let ((val (- (matrixRef data i 0) mn)))
          (set! sumsq (+ sumsq (* val val)))))
      (/ sumsq (- n 1)))))

(define variance
  (lambda (data)
    (let* ((n (matrixGetM data))
           (sumsq zero)
           (mn zero))
      (dotimes (i n)
        (set! mn (+ mn (matrixRef data i 0))))
      (set! mn (/ mn n))
      (dotimes (i n)
        (let ((val (- (matrixRef data i 0) mn)))
          (set! sumsq (+ sumsq (* val val)))))
      (/ sumsq (- n 1)))))

(define maybeMergePopulations
  (lambda (pop1 pop2 numpoints plot?)
    (let* ((cac (clusterAssignmentCost (car pop1) (car pop2) numpoints))
           (set (mergeClusters pop1 pop2))
           (mv (meanAndCovariance (cadr set)))
           (eigs (computeSortedEigenvectorsFromCovarianceMatrix (cadr mv))))
      (callWithCurrentContinuation
        (lambda (return)
          (do ((eignum 0 (+ 1 eignum)))
              ((>= eignum (matrixGetN (car eigs)))
               set)
              (let* ((hyperplane (hyperplaneEquationPerpendicularForm (car mv) eigs eignum))
                     (cdists (computeDistanceList (cadr set) (car set) hyperplane))
                     (p1dists (computeDistanceList (cadr pop1) (car pop1) hyperplane))
                     (p2dists (computeDistanceList (cadr pop2) (car pop2) hyperplane))
                     (cvar (variance cdists)) ;(varGivenMean zero cdists)
                     (p1var (variance p1dists))
                     (p2var (variance p2dists))
                     (ccost (clusterDescriptionLength cdists 0 (car set) zero cvar))
                     (p1cost (clusterDescriptionLength p1dists 0 (car pop1) zero p1var))
                     (p2cost (clusterDescriptionLength p2dists 0 (car pop2) zero p2var)))
                ;(format "MMP: ccost=~a p1cost=~a p2cost=~a cac=~a  merge? ~a~%" ccost p1cost p2cost cac (not (< (+ p1cost p2cost cac) ccost)))
                (if (< (+ p1cost p2cost cac) ccost) (return #f)))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                                    splitAt
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define splitAt
  (lambda (data n objs order bestindex eigs eignum plot?)
    (let* ((index0 (- (inexact->exact (matrixRef order bestindex 0)) 1))
           (index1 (- (inexact->exact (matrixRef order (+ bestindex 1) 0)) 1))
           (numDims (length data))
           (x0 (matrixRef (car data) index0 0))
           (y0 (matrixRef (cadr data) index0 0))
           (x1 (matrixRef (car data) index1 0))
           (y1 (matrixRef (cadr data) index1 0))
           (mx (/ (+ x0 x1) 2))
           (my (/ (+ y0 y1) 2))
           (evecs (car eigs))                ; display only
           (evals (cadr eigs))                ; display only
           (oen (if (zero? eignum) 1 0))
           (sd1 (sqrt (matrixRef evals oen oen))) ; display only
           (ev1x (matrixRef evecs 0 oen))        ; display only
           (ev1y (matrixRef evecs 1 oen))) ; display only
      ;; display only
      (when (and (= numDims 2) plot?)
        ;(mlEval (format #f "plot([~a],[~a],'kh');" mx my))
        (mlEval (format #f "plot([~a, ~a],[~a,~a],'-k');"
                        (- mx (* 4 sd1 ev1x)) (+ mx (* 3 sd1 ev1x))
                        (- my (* 4 sd1 ev1y)) (+ my (* 3 sd1 ev1y))))
        (mlEval "legend('Data', 'Principal Eigenvector', 'Secondary Eigenvector', 'Division Line',-1)"))
      ;; Now split the data into the two sets [0 .. bestindex] [bestindex+1 .. n]
      (let* ((size0 (+ bestindex 1))
             (size1 (- n (+ bestindex 1)))
             (obj0 (makeVector size0))
             (obj1 (makeVector size1))
             (d0 nil)
             (d1 nil))
        (dotimes (i numDims)
          (set! d0 (cons (makeZeroMatlabMatrix size0 1) d0))
          (set! d1 (cons (makeZeroMatlabMatrix size1 1) d1)))
        (dotimes (i size0)
          (let ((index (- (inexact->exact (matrixRef order i 0)) 1)))
            (vectorSet! obj0 i (vectorRef objs index))))
        (dotimes (i size1)
          (let ((index (- (inexact->exact (matrixRef order (+ i size0) 0)) 1)))
            (vectorSet! obj1 i (vectorRef objs index))))
        (do ((dd0 d0 (cdr dd0))
             (dd1 d1 (cdr dd1))
             (dds data (cdr dds)))
            ((null? dd0) #f)
           (dotimes (i size0)
             (let ((index (- (inexact->exact (matrixRef order i 0)) 1)))
               (matrixSet! (car dd0) i 0 (matrixRef (car dds) index 0))))
           (dotimes (i size1)
             (let ((index (- (inexact->exact (matrixRef order (+ i size0) 0)) 1)))
               (matrixSet! (car dd1) i 0 (matrixRef (car dds) index 0)))))
        (if (> size0 size1)
            (list (list size0 d0 obj0) (list size1 d1 obj1))
            (list (list size1 d1 obj1) (list size0 d0 obj0)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;                                  Categorize
;;;
;;; Principal component decomposition
;;;  Args:
;;;    data - a list of column matrices.  Each has all points for a
;;;           particular dimension.
;;;    objects - a list of objects.
;;;    iterations - the maximum number of iterations before a forced
;;;           stop.
;;;    plot? - #f if no plot required.  Otherwise produce 2D plots.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define categorize
  (lambda (data objects iterations plot?)
    (let* ((numpoints (matrixGetM (car data))) ; number of points from first dim.
            ;; NEW: starting point, should be empty at the end
           (categories (list (list numpoints data objects))) ; data divided into categories 
           (categorized ())) 
      ;; Principal Component cluster division
      (do ((it 0 (+ it 1)))
          ((or (>= it iterations) (null? categories))
           ;; orphan merging
		   ;; 1. Sort clusters based on its length
           (let ((res (map
                        (lambda (cluster)
                          (let ((mv (meanAndCovariance (cadr cluster))))
                            (list (car cluster)     ; Number of points
                                  (car mv)          ; Mean
                                  (cadr mv)         ; Variance
                                  (cadr cluster)    ; Data Points
                                  (caddr cluster)))); The objects 
                        (bsort (append categories categorized) head>)))
                 (merged nil))
             ;; Each element of res is (points mean invcovar (dims))
             ;; Sorted in order of increasing size.
			 
			
             (do ((clusters res (cdr clusters)))
                 (callWithCurrentContinuation
                   (lambda (return)
				     ;; 2. Loop through the 1st cluster vs rest, sort the rest based on the distance to the 1st cluster
                     (let* ((cands (map (lambda (cl)
                                         (list
                                           ((if (> (car cl) 2)
                                               mahalanobisDistance
                                               bigDistance)
                                            (cadr (car clusters)) ;mean(x)
                                            (cadr cl) ; mean(y)
                                            (matrixInverse (caddr cl))) ; invcovar(y)
                                           cl))
                                       (cdr clusters)))
                            (scands (map cadr (bsort cands head>))))
						;; 3. Loop through the 1st cluster vs the sorted rest (scands), try to merge them
                       (do ((mergers scands (cdr mergers)))
                           ((null? mergers) 
                            (push! (list (car (car clusters))           ; number of points
                                         (cadddr (car clusters))        ; data
                                         (car (cddddr (car clusters)))) ; objects
                                   merged)
                            (return #f))
                           (let ((mc (maybeMergePopulations
                                       (list (car (car clusters))
                                             (cadddr (car clusters))
                                             (car (cddddr (car clusters))))
                                       (list (car (car mergers))
                                             (cadddr (car mergers))
                                             (car (cddddr (car mergers))))
                                       0
                                       plot?)))
									   
							 ;; 4 If merge is true, we replace the 2 clusters with the merged cluster, repeat the whole process.
                             (when mc
                               ;(format "Orphan merged~%")
                               (push! mc merged)
                               (setCdr! clusters (delete (car mergers) (cdr clusters)))
                               (return #f))))))))))
							   
          (let* ((biggestset (car categories))
                 (setsize (car biggestset)) # num of points
                 (set (cadr biggestset))  # points
                 (setobj (caddr biggestset)))
            (set! categories (cdr categories)) # remove 1st from current array
            (let* ((mv (meanAndCovariance set))
                   (eigs (computeSortedEigenvectorsFromCovarianceMatrix (cadr mv))))
              (let ((split (maybeDividePopulation set setsize setobj mv eigs numpoints plot?)))
                (if split
                    ;; We split the data.  Add it to the list of categories and resort.
                    (set! categories (nreverse (bsort (append split categories) head>)))
                    ;; We didn't split the data add the set to the 'done' pile
                    (set! categorized (cons biggestset categorized))))))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;Visualization aids.                                               ;;;
;;;                                                                  ;;;
;;; A collection of procedures for displaying data in Matlab for:    ;;;
;;;   a. Debugging algorithms                                        ;;;
;;;   b. Demonstrating performance                                   ;;;
;;;   c. Producing publishable results                               ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Plot 2D test data as a scatter plot

(define plot2Ddata
  (lambda (data color title)
    (destructuringBind (dx dy) data
      (mlEval (format #f "plot(~a,~a, '~a');"              ;"scatter(~a,~a);"
                      (withSlots (name) dx name)
                      (withSlots (name) dy name)
                      color))
      (if title (mlEval (format #f "title('~a');" title))))))

;;; Plots the 2D eigenvectors +/- 1 standard deviation centered on the mean.
;;; The primary axes is plotted in green and the secondary axis is plotted in blue.

(define plot2DEigenvectors
  (lambda (means eigs)
    (let* ((evecs (car eigs))
           (evals (cadr eigs))
           (meanx (matrixRef means 0 0))
           (meany (matrixRef means 1 0))
           (sd0 (sqrt (matrixRef evals 0 0)))
           (sd1 (sqrt (matrixRef evals 1 1)))
           (ev0x (matrixRef evecs 0 0))
           (ev0y (matrixRef evecs 1 0))
           (ev1x (matrixRef evecs 0 1))
           (ev1y (matrixRef evecs 1 1)))
      (mlEval "hold on;")
      (mlEval (format #f "plot([~a, ~a],[~a,~a],'-g');"
                      (- meanx (* sd0 ev0x)) (+ meanx (* sd0 ev0x))
                      (- meany (* sd0 ev0y)) (+ meany (* sd0 ev0y))))
      (mlEval (format #f "plot([~a, ~a],[~a,~a],'-b');"
                      (- meanx (* sd1 ev1x)) (+ meanx (* sd1 ev1x))
                      (- meany (* sd1 ev1y)) (+ meany (* sd1 ev1y))))
      (mlEval "legend('Data', 'Principal Eigenvector', 'Secondary Eigenvector',-1)"))))

(define *plotColorDatabase* '("+r" "og" "xb" "*k" "sm" "dc" "vr" "^g" "<b" ">k" "pm" "hc"))
(define *plotColorDatabase2* '(".r" ".g" ".b" ".k" ".m" ".c" ".r" ".g" ".b" ".k" ".m" ".c"))

(define plotCategories
  (lambda (categories)
    (let ((colors *plotColorDatabase2*))
      (dolist (c categories)
        (mlEval "hold on")
        (plot2Ddata (cadr c) (car colors) #f)
        (set! colors (cdr colors))))))        

;;; Fin




