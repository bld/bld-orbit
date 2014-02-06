;;; Force related functions

(in-package :bld-orbit)

;;; Acceleration functions

(defun solar-pressure (rm b)
  "Inverse square solar pressure as function of radius from body and body object"
  (/ (* (slot-value b 'ls) (expt (/ *au* rm) 2)) *c*))

(defun gravity (s r mu)
  "Gravitational acceleration"
  (- (* (/ mu (norme2 r)) (unitg r))))

(defun sail-ideal-acc (s x sc)
  "Ideal solar sail acceleration"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (lightness cb pointfun) sc
      (with-slots (mu) cb
	(let ((n (funcall pointfun s x sc)))
	  (* lightness mu (/ (norme2 r))
	     (expt (scalar (*i (unitg r) n)) 2)
	     n))))))

;;; Pointing functions

(defun sail-pointing-normal (s x sc)
  "Return sail normal vector from unit sun vector"
  (unitg (slot-value (to-cartesian x s sc) 'r)))

(defun sail-pointing-fixed (s x sc)
  "Return sail normal vector from fixed RVBASIS rotor RS"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (rs basis) sc
      (let* ((rvbasis (rvbasis r v))
	     (rvr (recoverrotor rvbasis basis))
	     (rsrv (*g rvr rs)))
	(rotateg (first basis) rsrv)))))

(defun sail-frame-fixed (s x sc)
  "Return sail frame (and rotor) from fixed RVBASIS rotor RS"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (rs basis) sc
      (let* ((rvbasis (rvbasis r v))
	     (rvr (recoverrotor rvbasis basis))
	     (rsrv (*g rvr rs)))
	(values (new-frame rsrv basis) rsrv)))))
  
(defun sail-pointing-table (s x sc)
  "Return sail normal vector from lookup table of RVBASIS rotors RS"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (rs basis) sc
      (let* ((rvbasis (rvbasis r v))
	     (rvr (recoverrotor rvbasis basis))
	     (rsi (second (find s rs :test #'<= :key #'first))))
	(rotateg (first basis) (*g rvr rsi))))))

(defun sail-frame-table (s x sc)
  "Return sail frame from lookup table of RVBASIS rotors RS"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (rs basis) sc
      (let* ((rvbasis (rvbasis r v))
	     (rvr (recoverrotor rvbasis basis))
	     (rsi (second (find s rs :test #'<= :key #'first)))
	     (rsframe (*g rvr rsi)))
	(values (new-frame rsail rsframe) rsframe)))))

;;; Forces

