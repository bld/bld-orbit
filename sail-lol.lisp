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

(defparameter *x0* (make-hash :r (ve2 :e1 1) :v (ve2 :e2 1)))

(defpan gravity (x) (mu)
  (lethash (r) x
    (- (/ (* mu (unitg r)) (norme2 r)))))

(let ((mu 1)
      g)
  (defpfun kepler (tm x &optional p) (mu g)
    (lethash (r v) x
      (make-hash
       :r v
       :v (setf g (gravity x #'kepler))))))

(defpan ideal-sail-normal (x) (mu lightness)
  (lethash (r v) x
    (* lightness mu (/ (norme2 r))
       (unitg r))))

(let ((mu 1)
      (lightness 0.1)
      g
      a)
  (defpfun normal (tm x &optional p) (mu lightness g a)
    (lethash (r v) x
      (make-hash
       :r v
       :v (+ (setf g (gravity x #'normal))
	     (setf a (ideal-sail-normal x #'normal)))))))

(defpan ideal-sail-acc (x) (lightness mu rs)
  (lethash (r) x
    (let* ((ru (unitg r))
	   (n (rotateg ru rs)))
      (* lightness mu (/ (norme2 r))
	 (expt (scalar (*i ru n)) 2)
	 n))))

(let ((mu 1)
      (lightness 0.1)
      (rs (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2))))))
  (defpfun fixed (tm x &optional p) (mu lightness rs)
    (lethash (r v) x
      (make-hash
       :r v
       :v (+ (gravity x #'fixed)
	     (ideal-sail-acc x #'fixed))))))

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
