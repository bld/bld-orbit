(in-package :bld-orbit)

(defun gravity2 (r b)
  (with-slots (mu) b
    (- (* (/ mu (norme2 r)) (unitg r)))))

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
      (/ (* 2 p area (expt (*i ru n) 2) n) mass))))

(defmethod eom2 (s x sc)
  (with-slots (cb nbodies) sc
    ;; position & velocity of SC relative to central body
    (with-slots (r v) (to-cartesian x s sc)
      (let* ((tm (time-of s x)) ; UTC
	     (r-cb (position-velocity cb tm)) ; Position of central body
	     (r-sc (+ r-cb r)) ; position of SC in system
	     (r-nbodies ; positions of other bodies in system
	      (mapcar 
	       #'(lambda (b)
		   (slot-value (position-velocity b tm) 'r))
	       nbodies))
	     (rsc-nbodies ; positions of SC relative to other bodies
	      (mapcar 
	       #'(lambda (rb) (- r-sc rb))
	       r-nbodies))
	     (g-cb (gravity2 r cb)) ; gravity of central body
	     (g-nbodies ; gravities of other bodies
	      (mapcar 
	       #'(lambda (rscb b) (gravity2 rscb b))
	       rsc-nbodies nbodies))
	     (sframe (funcall pointfun s x sc)) ; sail frame
	     (fs-cb 
	     
