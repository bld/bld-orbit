;;; Force related functions

(in-package :bld-orbit)

;;; Acceleration functions

(defun solar-pressure (tm x sc)
  "Solar pressure for given time, state, and spacecraft data"
  (with-slots (sun) sc
    (with-slots (ls) sun
      (with-derived (rm) tm x sc
	(/ (* ls (expt (/ *au* rm) 2)) *c*)))))

(defun gravity (tm x sc)
  "Gravitational acceleration from central body"
  (with-derived (rm2 ru) tm x sc
    (with-slots (cb) sc
      (with-slots (mu) cb
	(- (* (/ mu rm2) ru))))))

;;; Forces

(defun no-sail (tm x sc)
  "No sail acceleration"
  (newg (r x)))

(defun sail-ideal-acc-normal (tm x sc)
  "Ideal sail pointing directly at the sun"
  (with-slots (area mass) sc
    (with-derived (p ru-sc-sun) tm x sc
      (/ (* 2 p area ru-sc-sun) mass))))

(defun sail-ideal-acc-fixed (tm x sc)
  (with-slots (area mass) sc
    (with-derived (p ru-sc-sun n) tm x sc
      (/ (* 2 p area (expt (scalar (*i ru-sc-sun n)) 2) n) mass))))

#+null(defun sail-flat-optical-acc (rsail-sun sframe sail sun)
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
