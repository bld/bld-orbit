;;; Solar sail dynamics using Lamba-Over-Let Pandoric closures
(defpackage :sail-lol
  (:use :cl :lol :bld-ga :bld-e2 :bld-e3 :bld-ode :bld-utils)
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

(defparameter *x0* (plist-hash-table (list :r (ve2 :e1 1) :v (ve2 :e2 1))))

(defun gravity (x eom)
  (with-pandoric (mu) eom
    (lethash (r) x
      (- (/ (* mu (unitg r)) (norme2 r))))))

(let ((mu 1))
  (defpfun kepler (tm x &optional p) (mu)
    (lethash (r v) x
      (plist-hash-table
       (list :r v
	     :v (gravity x #'kepler))))))

(defun ideal-sail-normal (x eom)
  (with-pandoric (mu lightness) eom
    (lethash (r v) x
      (* lightness mu (/ (norme2 r))
	 (unitg r)))))

(let ((mu 1)
      (lightness 0.1))
  (defpfun normal (tm x &optional p) (mu lightness)
    (lethash (r v) x
      (plist-hash-table
       (list
	:r v
	:v (+ (gravity x #'normal)
	      (ideal-sail-normal x #'normal)))))))

(defun ideal-sail-acc (x eom)
  (with-pandoric (lightness mu rs) eom
    (lethash (r) x
      (let* ((ru (unitg r))
	     (n (rotateg ru rs)))
	(* lightness mu (/ (norme2 r))
	   (expt (scalar (*i ru n)) 2)
	   n)))))

(let ((mu 1)
      (lightness 0.1)
      (rs (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2))))))
  (defpfun fixed (tm x &optional p) (mu lightness rs)
    (lethash (r v) x
      (plist-hash-table
       (list
	:r v
	:v (+ (gravity x #'fixed)
	      (ideal-sail-acc x #'fixed)))))))

