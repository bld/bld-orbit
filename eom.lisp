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
      (/ (* 2 p area (expt (scalar (*i ru n)) 2) n) mass))))

(defun no-sail (rsail-sun sframe sail sun) (newg rsail-sun))

(defgeneric dxds (s x sc f g-cb))

(defmethod eom2 (s x sc)
  (with-slots (cb sun pointfun accfun) sc
    ;; position & velocity of SC relative to central body
    (with-slots (r v) (to-cartesian x s sc)
      (let* ((tm (time-of s x)) ; UTC
	     (r-cb (slot-value (position-velocity cb tm) 'r)) ; State of central body
	     (r-sc (+ r-cb r)) ; position of SC in system
	     (r-sun (slot-value (position-velocity sun tm) 'r))
	     (rsc-sun (- r-sc r-sun)) ; position of SC wrt sun
	     (g-cb (gravity2 r cb)) ; gravity of central body
	     (sframe (funcall pointfun s x sc)) ; sail frame
	     (f-sun (funcall accfun rsc-sun sframe sc sun))) ; sun sail force
	(dxds s x sc f-sun g-cb))))) ; derivative of state

(defmethod dxds (tm (x cartstate) sc f g-cb)
  (with-slots (r v) x
    (make-instance
     'cartstate
     :r v
     :v (+ f g-cb))))

(defmethod dxds (s (x ksstate) sc f g-cb)
  (with-slots (alpha beta e tm) x
    (let* ((xspin (to-spinor x s sc))
	   (xcart (to-cartesian xspin s sc))
	   (r (slot-value xcart 'r))
	   (v (slot-value xcart 'v))
	   (x0 (slot-value sc 'x0))
	   (e0 (slot-value x0 'e))
	   (hk0 (- e0))
	   (w0 (sqrt (/ hk0 2)))
	   (rm (norme r))
	   (u (slot-value xspin 'u))
	   (hk (- e))
	   (w (/ (- hk hk0) 2))
	   (ff (- (/ (*g f r u) 2) (* w u))))
      (make-instance
       'ksstate
       :alpha (- (* ff (/ (sin (* w0 s)) w0)))
       :beta (* ff (/ (cos (* w0 s)) w0))
       :e (* rm (scalar (*i v f)))
       :tm rm))))


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
	     (rsi (second (find s rs :test #'<= :key #'first)))
	     (rsframe (*g rvr rsi)))
	(values (new-frame rsframe basis) rsframe)))))
