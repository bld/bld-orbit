(ql:quickload :lol)
(ql:quickload :bld-utils)
(ql:quickload :unit-formulas)
(ql:quickload :bld-e3)
(ql:quickload :bld-ga)
(ql:quickload :bld-gen)

(defpackage :chebyshev
  (:use :cl :lol :bld-utils :unit-formulas :bld-ga :bld-e3)
  (:shadowing-import-from :bld-gen
    + - * / expt
    sin cos tan
    atan asin acos
    sinh cosh tanh
    asinh acosh atanh
    log exp sqrt abs
    min max signum))

(in-package :chebyshev)

(defmethod print-object ((o hash-table) stream)
  (printhash o stream))

(defclass record ()
  ((start :initarg :start)
   (duration :initarg :duration)
   (xcoeffs :initarg :xcoeffs)
   (ycoeffs :initarg :ycoeffs)
   (zcoeffs :initarg :zcoeffs)))

(defmethod rv ((record record) date)
  (with-slots (start duration xcoeffs ycoeffs zcoeffs) record
    "Position and velocity from Chebyshev polynomials"
    (let* ((tm (/ (- (* 2 (- date start))
		     duration)
		  duration))
	   (twotm (* 2 tm))
	   (pkm1 tm)
	   (pk tm)
	   (xp (aref xcoeffs 0))
	   (yp (aref ycoeffs 0))
	   (zp (aref zcoeffs 0))
	   (qkm1 0)
	   (qk 1)
	   (xv 0)
	   (yv 0)
	   (zv 0)
	   pkm2
	   qkm2)
      (loop for k from 1 below (length xcoeffs)
	 do (incf xp (* (aref xcoeffs k) pk))
	   (incf yp (* (aref ycoeffs k) pk))
	   (incf zp (* (aref zcoeffs k) pk))
	   (incf xv (* (aref xcoeffs k) qk))
	   (incf yv (* (aref ycoeffs k) qk))
	   (incf zv (* (aref zcoeffs k) qk))
	   (setf pkm2 pkm1)
	   (setf pkm1 pk)
	   (setf pk (- (* twotm pkm1) pkm2))
	   (setf qkm2 qkm1)
	   (setf qkm1 qk)
	   (setf qk (+ (* twotm qkm1) (* 2 pkm1) (- qkm2)))
	 finally (return
		   (let ((vscale (/ 2 duration)))
		     (make-hash
		      :r (ve3 :e1 xp :e2 yp :e3 zp)
		      :v (ve3 :e1 (* xv vscale) :e2 (* yv vscale) :e3 (* zv vscale)))))))))

(defmethod import-record (record file position order)
  (with-open-file (s file)
    ;; Skip until specified record is reached
    (loop until (equal (read s) record))
    ;; Skip number after record
    (read s)
    ;; Read in Julian start/stop dates
    (let ((jd0 (read s))
	  (jdf (read s)))
      ;; Skip until position of body is reached
      (loop repeat (- position 3) do (read s))
      ;; Read in coefficients, collecting into X, Y, and Z
      (loop repeat order
	 collect (read s) into x
	 collect (read s) into y
	 collect (read s) into z
	 finally (return
		   (make-instance
		    'record
		    :start jd0
		    :duration (- jdf jd0)
		    :xcoeffs (make-array order :initial-contents x)
		    :ycoeffs (make-array order :initial-contents y)
		    :zcoeffs (make-array order :initial-contents z)))))))

(defparameter *mercury1* (import-record 1 "p:/src/ephemeris/ascp1900.421" 3 14))

(defparameter *earth1* (import-record 1 "p:/src/ephemeris/ascp1900.421" 231 13))
