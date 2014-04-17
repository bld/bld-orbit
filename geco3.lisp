(in-package :bld-orbit)

;; SIA chromosome

(defclass sia-chromosome (chromosome)
  ()
  (:documentation "Sun incidence angles array chromosome"))

(defmethod size ((self sia-chromosome))
  10)

(defmethod locus-arity ((self sia-chromosome) locus-index)
  90)

(defmethod loci-printable-form ((self sia-chromosome))
  (loci self))

;; Turn SIA organism to list of angles in radians

(defmethod sia-o-to-rad ((sia-o sia-organism))
  (let ((sia-c (first (genotype sia-o))))
    (loop for sia across (loci sia-c)
       collect (* (/ pi 180) (- sia (/ (locus-arity sia-c 0) 2d0))))))

;; SIA organism

(defclass sia-organism (organism)
  ()
  (:documentation "SIA array organism"))

(defmethod chromosome-classes ((self sia-organism))
  '(sia-chromosome))

;; SIA population

(defclass sia-population (generational-population minimizing-score-mixin)
  ()
  (:documentation "Population of SIA lookup array organisms with minimized scores"))

(defmethod organism-class ((self sia-population))
  'sia-organism)

;; SIA genetic plans: minimize state contraints

;; Evaluate SIA organism

;; Example problem: polar 2D eom

(defvar *sia* 0)

(defparameter *sail-polar*
  (make-hash :be 0.034 :mu 1))

(defparameter *x0* (make-hash :r 1 :th 0 :vr 0 :vt 1))

(defparameter *tspan* '(0 2 4 6 8 10 12 14 16 18 20))

(defun dx (tm x p)
  (let* ((sia (* (/ pi 180) *sia*))
	 (co (cos sia))
	 (si (sin sia)))
    (with-keys (be mu) p
      (with-keys (r th vr vt) x
	(make-hash
	 :r vr
	 :th (/ vt r)
	 :vr (+ (/ (expt vt 2) r)
		(/ (* mu (1- (* be (expt co 2) (abs co))))
		   (expt r 2)))
	 :vt (- (/ (* mu be (expt co 2) si)
		   (expt r 2))
		(/ (* vr vt)
		   r)))))))

(defun prop-sail-polar (sia-o)
  (loop for *sia* in (sia-o-to-rad sia-o)
     for i = 0 then (incf i)
     for t0 = (nth i *tspan*)
     for tf = (nth (1+ i) *tspan*)
     for x0 = *x0* then (second (car (last traj)))
     for traj = (rka #'dx t0 tf x0 :hmax (- tf t0) :param *sail-polar*)
     collect traj))

(defun export-polar-traj-data (fname sia-o)
  (let ((sias (sia-o-to-rad sia-o))
	(trajs (prop-sail-polar sia-o)))
    (with-open-file (s fname :direction :output :if-exists :supersede)
      (format s "# Time SIA R TH VR VT~%")
      (loop for sia in sias
	 for traj in trajs
	 do (loop for (tm x) in traj
	       do (with-keys (r th vr vt) x
		    (format 
		     s "~&~a"
		     (substitute 
		      #\E #\d
		      (format nil "~a ~a ~a ~a ~a ~a" tm sia r th vr vt)))))))))

