(in-package :bld-orbit)

;;; Classical orbital element conversion

(defun rv2coe (rv vv mu basis)
  "Return classical orbital elements given position & velocity vectors, gravitational parameter, and a list of 3 basis vectors"
  (flet ((energy-rv (r v mu) ; Orbit energy from position, velocity, & gravitational parameter
	   (- (/ (* v v) 2)
	      (/ mu r)))
	 (sma-rv (r v mu) ; semi-major axis from position, velocity, & gravitational parameter
	   (/ mu (* -2 (energy-rv r v mu))))
	 (eccv-rv (rv vv mu) ; eccentricity vector from position, velocity, & gravitational parameter
	   (/ (- (* rv (- (norme2 vv) (/ mu (norme rv))))
		 (* vv (scalar (*i rv vv))))
	      mu))
	 (mombv-rv (rv vv) ; orbital momentum bivector from position & velocity
	   (*o rv vv))
	 (nodev-rv (mombv basis) ; ascending node vector from momentum and basis
	   (- (*i (third basis) mombv)))
	 (inc-rv (mombv basis) ; inclination from momentum bivector and basis
	   (acos (scalar (dual (*o (third basis) (unitg mombv))))))
	 (lan-rv (nodev basis) ; Longitude of ascending node from node vector and basis
	   (let ((tmp (atan (scalar (*i nodev (second basis)))
			    (scalar (*i nodev (first basis))))))
	     (if (< tmp 0) (+ tmp (* 2 pi)) tmp)))
	 (aop-rv (nodev eccv mombv) ; argument of perigee from node vector, eccentricity vector, and momentum bivector
	   (let* ((n (unitg nodev))
		  (e (unitg eccv))
		  (h (unitg mombv))
		  (ov (*i n h))
		  (tmp (atan (scalar (*i e ov)) (scalar (*i e n)))))
	     (if (< tmp 0) (+ tmp (* 2 pi)) tmp)))
	 (truan-rv (eccv rv mombv) ; true anomaly from eccentricity, position, and velocity vectors"
	   (let* ((e (unitg eccv))
		  (h (unitg mombv))
		  (ruv (unitg rv))
		  (ov (*i e h))
		  (tmp (atan (scalar (*i ruv ov))
			     (scalar (*i ruv e)))))
	     (if (< tmp 0) (+ tmp (* 2 pi)) tmp))))
    (let* ((mombv (mombv-rv rv vv))
	   (nodev (nodev-rv mombv basis))
	   (eccv (eccv-rv rv vv mu)))
      (make-hash
       :sma (sma-rv (norme rv) (norme vv) mu)
       :ecc (norme eccv)
       :inc (inc-rv mombv basis)
       :lan (raan-rv nodev basis)
       :aop (aop-rv nodev eccv mombv)
       :truan (truan-rv eccv rv mombv)))))

(defun coe2rv (sma ecc inc lan aop truan mu basis)
  "Convert cassical orbital elements to position & velocity vectors.
SMA semi major axis
ECC eccentricity
INC inclination
LAN longitude of ascending node
AOP argument of perigee
TRUAN true anomaly
MU gravitational parameter
BASIS list of 3 orthogonal basis vectors to express position & velocity in"
  (let* ((r-lan (rotor (*o (first basis) (second basis)) lan))
	 (basis-lan (new-frame r-lan basis))
	 (r-inc (rotor (- (dual (first basis-lan))) inc))
	 (basis-inc (new-frame r-inc basis-lan))
	 (r-aop (rotor (- (dual (third basis-inc))) aop))
	 (basis-aop (new-frame r-aop basis-inc))
	 (r-truan (rotor (- (dual (third basis-aop))) truan))
	 (basis-truan (new-frame r-truan basis-aop))
	 (p (* sma (- 1 (expt ecc 2))))
	 (r (/ p (+ 1 (* ecc (cos truan)))))
	 (rv (* r (first basis-truan)))
	 (vv (* (sqrt (/ mu p))
		(- (* (+ ecc (cos truan)) (second basis-aop))
		   (* (sin truan) (first basis-aop))))))
    (values rv vv)))

;; COE equations of motion

(defclass coestate ()
  ((a :initarg :a :documentation "Semi-major axis")
   (e :initarg :e :documentation "Eccentricity")
   (i :initarg :i :documentation "Inclination")
   (lan :initarg :lan :documentation "Longitude of the ascending node")
   (aop :initarg :aop :documentation "Argument of periapsis")
   (tm :initarg :tm :documentation "Time"))
  (:documentation "Equations of motion for classical orbital element state with true anomaly as independent variable"))

(defstatearithmetic coestate (a e i lan aop tm))

(let ((slots '(a e i lan aop tm)))
  (defmethod print-object ((x coestate) stream)
    (format stream "#<COESTATE ~{~a~^ ~}>"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value x slot))))))

(defmethod to-cartesian ((x coestate) f sc)
  (with-slots (a e i lan aop tm) x
    (with-slots (cb iframe) sc
      (with-slots (mu) cb
	(multiple-value-bind (r v) (coe2rv a e i lan aop f mu iframe)
	  (values
	   (make-instance
	    'cartstate
	    :r r
	    :v v)
	   tm))))))

(defmethod to-coe ((x cartstate) tm sc)
  "Convert cartesian to classical orbital element state"
  (with-slots (r v) x
    (with-slots (cb iframe) sc
      (with-slots (mu) cb
	(let* ((mombv (mombv-rv r v))
	       (nodev (nodev-rv mombv iframe))
	       (eccv (eccv-rv r v mu)))
	  (values
	   (make-instance
	    'coestate
	    :a (sma-rv (norme r) (norme v) mu)
	    :e (norme eccv)
	    :i (inc-rv mombv iframe)
	    :lan (lan-rv nodev iframe)
	    :aop (aop-rv nodev eccv mombv)
	    :tm tm)
	   (truan-rv eccv r mombv)))))))

(defmethod energy ((x coestate) sc)
  (with-slots (a) x
    (with-slots (cb) sc
      (with-slots (mu) cb
	(- (/ mu 2 a))))))

(defmethod eom (f (x coestate) sc)
  "Classical orbital element equations of motion. Watch for singularities."
  (with-slots (a e i lan aop tm) x
    (with-slots (cb accfun iframe) sc
      (with-slots (mu) cb
	(multiple-value-bind (r v) (coe2rv a e i lan aop f mu iframe) ; position & velocity
	  (let* ((rm (norme r)) ; radius
		 (p (* a (- 1 (expt e 2)))) ; semi latus rectum
		 (fv (funcall accfun f x sc)) ; force vector
		 (obasis (rv-frame f x sc)) ; orbit basis from position/velocity
		 (sf (scalar (*i fv (first obasis)))) ; radial force
		 (tf (scalar (*i fv (second obasis)))) ; transverse force
		 (wf (scalar (*i fv (third obasis)))) ; orbit normal force
		 (dlan (* (/ (expt rm 3) mu p (sin i)) ; derivative of longitude of ascending node
			  (sin (+ f aop))
			  wf)))
	    (make-instance
	     'coestate
	     :a (* (/ (* 2 p (expt rm 2))
		      mu (expt (- 1 (expt e 2)) 2))
		   (+ (* sf e (sin f)) (/ (* tf p) rm)))
	     :e (* (/ (expt rm 2) mu)
		   (+ (* sf (sin f))
		      (* tf (+ 1 (/ rm p)) (cos f))
		      (* tf (/ rm p) e)))
	     :i (* (/ (expt rm 3) mu p)
		   (cos (+ f aop))
		   wf)
	     :lan dlan
	     :aop (- (* (/ (expt rm 2) mu e)
			(- (* tf (+ 1 (/ rm p)) (sin f))
			   (* sf (cos f))))
		     (* dlan (cos i)))
	     :tm (* (/ (expt rm 2) (sqrt (* mu p)))
		    (- 1
		       (* (/ (expt rm 2) mu e)
			  (- (* sf (cos f))
			     (* tf (+ 1 (/ rm p)) (sin f)))))))))))))
