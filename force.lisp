;;; Force related functions

(in-package :bld-orbit)

;;; Acceleration functions

(defun solar-pressure (rm b)
  "Inverse square solar pressure as function of radius from body and body object"
  (/ (* (slot-value b 'ls) (expt (/ *au* rm) 2)) *c*))

(defun gravity (s r mu)
  "Gravitational acceleration"
  (- (* (/ mu (norme2 r)) (unitg r))))

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

(defun sail-frame-sun-normal (s x sc)
  "Point sail at the sun"
  (with-slots (r v) (to-cartesian x s sc)
    (rvbasis r v)))

(defun sail-frame-sun-fixed (s x sc)
  "Return sail frame from fixed RVBASIS rotor"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (rs basis) sc
      (let* ((rvbasis (rvbasis r v))
	     (rvr (recoverrotor rvbasis basis))
	     (rsframe (*g rvr rs)))
	(values (new-frame rsframe basis) rsframe)))))

(defun sail-frame-sun-table (s x sc)
  "Return sail frame from lookup table of RVBASIS rotors RS"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (rs basis t0 tf) sc
      (let* ((stmp (if (> s tf) tf s))
	     (rvbasis (rvbasis r v))
	     (rvr (recoverrotor rvbasis basis))
	     (rsi (second (find stmp rs :test #'<= :key #'first)))
	     (rsframe (*g rvr rsi)))
	(values (new-frame rsframe basis) rsframe)))))

;;; Forces

(defun sail-flat-optical-acc (rsail-sun sframe sail sun)
  "Flat plate sail with optical properties from 'Propulsive Reflectivity and Photoflexibility' by Derbes and Lichodziejewski (AIAA 2006-4520)"
  (with-slots (area mass optical) sail
    (with-slots (reft refp refp-0 thetap-25 ef eb bf bb) optical
      (let* ((n (first sframe))
	     (rm (norme r))
	     (nt (- (*x u n n)))
	     (p (solar-pressure rm sun))
	     (cossia (- (scalar (*i n u))))
	     (costhetap (* (/ (cos (rad thetap-25)) (cos (rad 25))) cossia))
	     (s (unitg (+ (* costhetap n)
			  (* (/ (sin (rad thetap-25))
				(sin (rad 25)))
			     nt))))
	     (fabs (* p area cossia u))
	     (femit (* -1 (- 1 reft) p area cossia
		       (- (* ef bf) (* eb bb))
		       (/ (+ ef eb))
		       n))
	     (fref (- (* p area cossia (* reft refp-0) s))))
	(+ fabs femit fref)))))

(defun sail-flat-ideal-acc (rsail-sun sframe sail sun)
  "Ideal flat sail acceleration as function of:
RSAIL-SUN: Sail position wrt sun
SFRAME: Sail frame
SAIL: Sail object
SUN: Sun object"
  (with-slots (area mass) sail
    (let* ((ru (unitg rsail-sun))
	   (n (first sframe))
	   (p (solar-pressure (norme rsail-sun) sun)))
      (/ (* 2 p area (expt (scalar (*i ru n)) 2) n) mass))))

(defun no-sail (rsail-sun sframe sail sun) (newg rsail-sun))

(defun sail-ideal-acc (s x sc)
  "Ideal solar sail acceleration"
  (with-slots (r v) (to-cartesian x s sc)
    (with-slots (lightness cb pointfun) sc
      (with-slots (mu) cb
	(let ((n (funcall pointfun s x sc)))
	  (* lightness mu (/ (norme2 r))
	     (expt (scalar (*i (unitg r) n)) 2)
	     n))))))

