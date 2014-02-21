;;; Solar sail dynamics using Lamba-Over-Let Pandoric closures
(defpackage :sail-lol
  (:use :cl :lol :bld-ga :bld-e2 :bld-e3 :bld-ode :bld-utils :anaphora)
  (:import-from :alexandria plist-hash-table)
  (:shadowing-import-from
   :bld-gen
   + - * / expt
   sin cos tan
   atan asin acos
   sinh cosh tanh
   asinh acosh atanh
   log exp sqrt abs
   min max signum))

(in-package :sail-lol)

(defmethod norminfx ((x g))
  "Infinity norm of geometric algebra state variable"
  (norminf x))

;; Bodies

(let ((mu 1.327585d20)
      (ls 1361d0))
  (defpfun sun (tm) (self mu ls)
    (ve3))
  (defpfun sun2d (tm) (self mu ls)
    (ve2)))

(let ((mu 3.9873878d14)
      (au 1.4959787d11))
  (defpfun earth (tm) (self mu au)
    (ve3 :e1 au))
  (defpfun earth2d (tm) (self mu ls)
    (ve2 :e1 au)))

;; Forces

(defpan gravity (x) (cb)
  (with-pandoric (mu) cb
    (lethash (r) x
      (- (/ (* mu (unitg r)) (norme2 r))))))

(defpan ideal-sail-normal (x) (cb lightness)
  (with-pandoric (mu) cb
    (lethash (r v) x
      (* lightness mu (/ (norme2 r))
	 (unitg r)))))

(defpan ideal-sail-fixed (x) (lightness cb rs)
  (with-pandoric (mu) cb
    (lethash (r) x
      (let* ((ru (unitg r))
	     (n (rotateg ru rs)))
	(* lightness mu (/ (norme2 r))
	   (expt (scalar (*i ru n)) 2)
	   n)))))

;; Cartesian equations of motion

(defpan cart-kepler-eom (x) (self gfun)
  (lethash (r v) x
    (make-hash
     :r v
     :v (funcall gfun x self))))

(defpan cart-eom (x) (self gfun afun)
  (lethash (r v) x
    (make-hash
     :r v
     :v (+ (funcall gfun x self)
	   (funcall afun x self)))))

;; 2D cartesian examples

(let* ((cb #'sun2d)
       (t0 0)
       (tf (* 60 60 24 365.25))
       (x0 (make-hash :r (ve2 :e1 (get-pandoric #'earth 'au)) :v (ve2 :e2 (sqrt (/ (get-pandoric #'sun 'mu) (get-pandoric #'earth 'au))))))
       (fname "lol-2d-kepler-eg.dat")
       (gfun #'gravity))
  (defpfun 2d-kepler-eg (tm x &optional p) (self t0 tf x0 fname gfun cb)
    (cart-kepler-eom x self)))

(let* ((cb #'sun2d)
       (lightness 0.1)
       (t0 0)
       (tf (* 60 60 24 365.25))
       (x0 (make-hash :r (ve2 :e1 (get-pandoric #'earth 'au)) :v (ve2 :e2 (sqrt (/ (get-pandoric #'sun 'mu) (get-pandoric #'earth 'au))))))
       (fname "lol-2d-normal-eg.dat")
       (gfun #'gravity)
       (afun #'ideal-sail-normal))
  (defpfun 2d-normal-eg (tm x &optional p) (self lightness t0 tf x0 fname gfun afun cb)
    (cart-eom x self)))

(let ((cb #'sun2d)
      (lightness 0.1)
      (rs (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2)))))
      (t0 0)
      (tf (* 60 60 24 365.25))
      (x0 (make-hash :r (ve2 :e1 (get-pandoric #'earth 'au)) :v (ve2 :e2 (sqrt (/ (get-pandoric #'sun 'mu) (get-pandoric #'earth 'au))))))
      (fname "lol-2d-fixed-eg.dat")
      (gfun #'gravity)
      (afun #'ideal-sail-fixed))
  (defpfun 2d-fixed-eg (tm x &optional p) (self lightness rs t0 tf x0 fname gfun afun cb)
    (cart-eom x self)))

;; 3D cartesian examples

(let* ((cb #'sun2d)
       (t0 0)
       (tf (* 60 60 24 365.25))
       (x0 (make-hash :r (ve3 :e1 (get-pandoric #'earth 'au)) :v (ve3 :e2 (sqrt (/ (get-pandoric #'sun 'mu) (get-pandoric #'earth 'au))))))
       (fname "lol-3d-kepler-eg.dat")
       (gfun #'gravity))
  (defpfun 3d-kepler-eg (tm x &optional p) (self t0 tf x0 fname gfun cb)
    (cart-kepler-eom x self)))

(let* ((cb #'sun2d)
       (lightness 0.1)
       (t0 0)
       (tf (* 60 60 24 365.25))
       (x0 (make-hash :r (ve3 :e1 (get-pandoric #'earth 'au)) :v (ve3 :e2 (sqrt (/ (get-pandoric #'sun 'mu) (get-pandoric #'earth 'au))))))
       (fname "lol-3d-normal-eg.dat")
       (gfun #'gravity)
       (afun #'ideal-sail-normal))
  (defpfun 3d-normal-eg (tm x &optional p) (self lightness t0 tf x0 fname gfun afun cb)
    (cart-eom x self)))

(let ((cb #'sun2d)
      (lightness 0.1)
      (rs (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2)))))
      (t0 0)
      (tf (* 60 60 24 365.25))
      (x0 (make-hash :r (ve3 :e1 (get-pandoric #'earth 'au)) :v (ve3 :e2 (sqrt (/ (get-pandoric #'sun 'mu) (get-pandoric #'earth 'au))))))
      (fname "lol-3d-fixed-eg.dat")
      (gfun #'gravity)
      (afun #'ideal-sail-fixed))
  (defpfun 3d-fixed-eg (tm x &optional p) (self lightness rs t0 tf x0 fname gfun afun cb)
    (cart-eom x self)))

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
