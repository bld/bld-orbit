;;; Solar sail dynamics using Lamba-Over-Let Pandoric closures

(in-package :sail-lol)

(defmethod norminfx ((x g))
  "Infinity norm of geometric algebra state variable"
  (norminf x))

;; Constants

(defparameter *c* 2.99792458e8 "Speed of light (m/s)")

;; Derived values

(defderived rm2 (norme2 (gethash :r x)))
(defderived rm (sqrt (rm2 self x)))
(defderived ru (/ (gethash :r x) (rm self x)))

;; Bodies

(let* ((mu 1.327585d20)
       (ls 1361d0)
       (au 1.4959787d11)
       (au2 (expt au 2)))
  (defpfun sun (tm) (self mu ls au au2)
    (ve2)))

;; Frame functions

(defun rv-basis-2d (x)
  (lethash (r v) x
    (let* ((h (unitg (*o r v)))
	   (x (unitg r))
	   (y (unitg (*i r h))))
      (list x y))))

(defun rv-basis-3d (x)
  (lethash (r v) x
    (let* ((h (unitg (*o r v)))
	   (x (unitg r))
	   (y (unitg (*i r h))))
      (list x y (dual h)))))

;; Forces

(defpan solar-pressure (tm x) (self ls au2)
  (lethash (r) x
    (let ((rsun (funcall self tm)))
      (/ (* ls (/ au2 (norme2 (- r rsun)))) *c*))))

(defpan gravity (tm x) (self mu)
  (lethash (r) x
    (let* ((rcb (funcall self tm))
	   (rcb-sc (- r rcb)))
      (- (/ (* mu (unitg rcb-sc)) (norme2 rcb-sc))))))

(defpan ideal-sail-normal (tm x) (cb mass area)
  (lethash (r) x
    (/ (* 2 (solar-pressure cb tm x) (unitg r) area)
       mass)))

(defpan ideal-sail-fixed (tm x) (mass area cb rs)
  (lethash (r) x
    (let* ((ru (unitg r))
	   (n (rotateg ru rs)))
      (/ (* 2 (solar-pressure cb tm x) area (expt (scalar (*i ru n)) 2) n)
	 mass))))

;; Cartesian equations of motion

(defpan cart-kepler-eom (tm x) (cb gfun)
  (lethash (r v) x
    (let ((dx (make-hash
	       :r v
	       :v (funcall gfun cb tm x))))
      dx)))

(defpan cart-eom (tm x) (self gfun afun cb)
  (lethash (r v) x
    (let ((dx (make-hash
	       :r v
	       :v (+ (funcall gfun cb tm x)
		     (funcall afun self tm x)))))
      dx)))

;; 2D cartesian examples

(defparameter *x02dcart* (make-hash :r (ve2 :e1 (get-pandoric #'sun 'au)) :v (ve2 :e2 (sqrt (/ (get-pandoric #'sun 'mu) (get-pandoric #'sun 'au))))))

(let ((cb #'sun)
      (t0 0)
      (tf (* 60 60 24 365.25))
      (x0 *x02dcart*)
      (fname "lol-2d-kepler-eg.dat")
      (gfun #'gravity)
      (derived (make-hash-table)))
  (defpfun 2d-kepler-eg (tm x &optional p) (self t0 tf x0 fname gfun cb derived)
    (cart-kepler-eom self tm x)))

(let ((cb #'sun)
      (mass 52)
      (area 1260)
      (t0 0)
      (tf (* 60 60 24 365.25))
      (x0 *x02dcart*)
      (fname "lol-2d-normal-eg.dat")
      (gfun #'gravity)
      (afun #'ideal-sail-normal)
      (bfun #'rv-basis-2d)
      (derived (make-hash-table)))
  (defpfun 2d-normal-eg (tm x &optional p) (self mass area t0 tf x0 fname gfun afun cb)
    (cart-eom self tm x)))

(let ((cb #'sun)
      (mass 52)
      (area 1260)
      (rs (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2)))))
      (t0 0)
      (tf (* 60 60 24 365.25))
      (x0 *x02dcart*)
      (fname "lol-2d-fixed-eg.dat")
      (gfun #'gravity)
      (afun #'ideal-sail-fixed)
      (bfun #'rv-basis-2d))
  (defpfun 2d-fixed-eg (tm x &optional p) (self mass area rs t0 tf x0 fname gfun afun cb)
    (cart-eom self tm x)))

;; Export

(defun write-cart-traj (file trajdata)
  "Write cartesian orbit data to format plotable by Gnuplot.
Columns:
1: Time
2,3,4: Position
5,6,7: Velocity"
  (with-open-file (s file :direction :output :if-exists :supersede)
    (format s "# Time X Y Z VX VY VZ~%")
    (loop for (tm x) in trajdata
       do (lethash (r v) x
	    (format s "~&~a"
		    (substitute #\E #\d
				(format nil "~a ~a ~a ~a ~a ~a ~a" 
					tm 
					(gref r :e1) (gref r :e2) (aif (gref r :e3) it 0)
					(gref v :e1) (gref v :e2) (aif (gref v :e3) it 0))))))))

;; Propagate a sail problem

(defpan propagate () (self t0 tf x0 fname)
  (let ((traj (rka self t0 tf x0 :hmax (- tf t0))))
    (when fname (write-cart-traj fname traj))
    traj))
