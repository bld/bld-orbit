#|
Orbital Mechanics Library
=========================

Includes:
* Cartesian coordinate equations of motion
* Kustaanheimo-Stiefel equations of motion cast into geometric algebra by Hestenes
* Simple Keplerian trajectories
* Solar sail trajectories
|#

(in-package :bld-orbit)

(defparameter *au* 1.49597870691d8 "Astronomical unit (km)")

(defparameter *c* 299792.458d0 "Speed of light (km/s)")

(defparameter *J2000* (list (ve3 :e1 1d0)
			    (ve3 :e2 1d0)
			    (ve3 :e3 1d0)) "J2000 inertial frame")

;;; Utility functions

(defun new-frame (rotor frame)
  "New reference frame from rotor & original frame"
  (mapcar #'(lambda (basisvector) (rotateg basisvector rotor)) frame))

(defun make-vector (coefs basis)
  (apply #'+
	 (loop for coef in coefs
	    for bv in basis
	    collect (* coef bv))))

(defmethod norminfx ((x g))
  "Infinity norm of geometric algebra state variable"
  (norminf x))

(defgeneric energy-cb (x cb) (:documentation "Specific orbital energy from state & central body"))

(defgeneric energy (x sc) (:documentation "Keplerian orbit energy"))

(defgeneric time-of (s x) (:documentation "Universal time"))

(defgeneric eom (s x &optional param) (:documentation "Equations of motion given independent variable S & state X"))

(defgeneric to-cartesian (x s sc) (:documentation "Convert a state to cartesian state"))

(defgeneric to-spinor (x s sc) (:documentation "Convert a state to spinor state"))

(defgeneric to-initial-spinor (x s sc) (:documentation "Convert to an initial spinor state"))

(defgeneric to-ks (x s sc) (:documentation "Convert to a Kustaanheimo-Stiefel orbit element state"))

(defgeneric to-initial-ks (x s sc) (:documentation "Generate initial Kustaanheimo-Stiefel orbit element state from another state & corresponding central body"))

;;; Rotor, spinor and basis functions

(defun recoverrotor3d (fs es)
  "Recover a basis given new and original basis vectors"
  (let* ((esr (apply #'recipbvs es))
	 (psi (apply #'+ 1 (mapcar #'(lambda (f er) (*g f er)) fs esr))))
    (unitg psi)))

(defun recoverrotor2d (fs es)
  "Recover a basis given new and original basis vectors"
  (let* ((esr (apply #'recipbvs es))
	 (psi (apply #'+ 2 (mapcar #'(lambda (f er) (*g f er)) fs esr))))
    (unitg psi)))

(defun recoverrotor (fs es)
  (cond
    ((apply #'= 3 (mapcar #'dimension (append fs es))) (recoverrotor3d fs es))
    ((apply #'= 2 (mapcar #'dimension (append fs es))) (recoverrotor2d fs es))
    (t (error "FS and ES must all be dimension 2 or 3"))))

(defun recoverspinor (r fs es)
  "Recover a spinor given orbit radius, new basis vectors, and original basis vectors"
  (* (recoverrotor fs es) (sqrt r)))

(defun rvbasis (r v)
  "Form a basis from the position and velocity vectors. First is the unit position vector. Second is the complement of the first in the orbit plane. Third (if it exists) is the orbit plane normal vector."
  (let* ((h (unitg (*o r v)))
	 (x (unitg r))
	 (y (unitg (*i r h))))
    (cond
      ((= (dimension r) (dimension v) 2) (list x y))
      ((= (dimension r) (dimension v) 3)
       (list x y (dual h)))
      (t (error "R and V must be 2 or 3 dimensional")))))

(defun orbit-frame (rsun basis)
  "Produce orbit frame from vector to sun"
  (let* ((o3 (unitg rsun))
	 (o1 (unitg (*i (third basis) (dual o3))))
	 (o2 (*i o3 (dual o1))))
    (list o1 o2 o3)))

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

;;; Sail classes

(defclass sail ()
  ((eom :initarg :eom :initform #'eom :documentation "Equations of motion")
   (cb :initarg :cb :initform (make-instance 'body :mu 1d0) :documentation "Central body")
   (sun :initarg :sun :documentation "Sun (if it's different from CB)")
   (accfun :initarg :accfun :initform #'sail-ideal-acc :documentation "Sail acceleration function")
   (pointfun :initarg :pointfun :initform #'sail-pointing-fixed :documentation "Sail pointing function")
   (lightness :initarg :lightness :initform 0)
   (mass :initarg :mass)
   (area :initarg :area)
   (optical :initarg :optical)
   (basis :initarg :basis :initform (list (ve2 :e1 1) (ve2 :e2 1)))
   (t0 :initarg :t0 :initform 0d0 :documentation "Initial time")
   (tf :initarg :tf :initform (* 2 pi) :documentation "Final time")
   (x0 :initarg :x0 :initform (make-instance 'cartstate :r (ve2 :e1 1) :v (ve2 :e2 1)))
   (rs :initarg :rs :initform (re2 :s 1d0) :documentation "Sail orientation rotor wrt orbital position frame")
   (outfile :initarg :outfile :initform "sail-2d-cart-kepler-eg.dat" :documentation "Output filename"))
  (:documentation "Solar sail orbit problem. Default is 2D Kepler with circular orbit, and units of AU and TU."))

(let ((slots '(eom cb accfun pointfun lightness basis t0 tf x0 rs)))
  (defmethod print-object ((sail sail) s)
    (format s "#<SAIL ~{~a~^~%  ~}>"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value sail slot))))))

;;; Propagate a trajectory

(defgeneric propagate (sc &key))

(defmethod propagate ((sc sail) &key (outfile (slot-value sc 'outfile)) (hmax-factor 1))
  "Propagate sailcraft trajectory. Default maximum stepsize is the difference between T0 and TF."
  (with-slots (eom t0 tf x0) sc
    (let ((results (rka eom t0 tf x0 :param sc :hmax (/ (- tf t0) hmax-factor))))
      (when outfile (write-cart-traj outfile (to-cartesian-traj results sc)))
      results)))

(defmethod propagate ((h hash-table) &key)
  "Propagate a hash table of objects & write data file"
  (maphash2 #'(lambda (k v) (propagate v)) h))

;;; Cartesian state

(defclass cartstate ()
  ((r :initarg :r :documentation "Position vector")
   (v :initarg :v :documentation "Velocity vector"))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a>" (slot-value x 'r) (slot-value x 'v)))

(defstatearithmetic cartstate (r v))

(defmethod energy-cb ((x cartstate) cb)
  (with-slots (mu) cb
    (with-slots (r v) x
      (- (/ (scalar (exptg v 2)) 2) (/ mu (norme r))))))

(defmethod energy ((x cartstate) sc)
  (with-slots (cb) sc
    (with-slots (mu) cb
      (with-slots (r v) x
	(- (/ (scalar (exptg v 2)) 2) (/ mu (norme r)))))))

(defmethod time-of (s (x cartstate))
  s)

(defmethod eom (tm (x cartstate) &optional sc)
  "Solar sail cartesian equations of motion"
  (with-slots (r v) x
    (with-slots (cb accfun) sc
      (with-slots (mu) cb
	(make-instance
	 'cartstate
	 :r v
	 :v (+ (funcall accfun tm x sc)
	       (gravity tm r mu)))))))

(defmethod to-cartesian ((x cartstate) tm sc)
  (values x tm))

;;; Spinor states

(defclass spinorstate ()
  ((u :initarg :u :documentation "Spinor of position: r = u e1 (revg u)")
   (duds :initarg :duds :documentation "Spinor derivative wrt s")
   (tm :initarg :tm :documentation "Time"))
  (:documentation "Spinor state"))

(defmethod print-object ((x spinorstate) stream)
  (format stream "#<SPINORSTATE :U ~a :DUDS ~a :TM ~a>" (slot-value x 'u) (slot-value x 'duds) (slot-value x 'tm)))

(defstatearithmetic spinorstate (u duds tm))

(defmethod energy-cb ((x spinorstate) cb)
  (with-slots (u duds) x
    (with-slots (mu) cb
      (/ (- (* 2 (norme2 duds)) mu) (norme2 u)))))

(defmethod energy ((x spinorstate) sc)
  "Keplerian orbit energy"
  (with-slots (u duds) x
    (with-slots (cb) sc
      (with-slots (mu) cb
	(/ (- (* 2 (norme2 duds)) mu) (norme2 u))))))

(defmethod time-of (s (x spinorstate))
  (slot-value x 'tm))

(defmethod eom (s (x spinorstate) &optional sc)
  "Spinor equations of motion: 2 dU/ds - E U = f r U, dt/ds = |r|"
  (with-slots (u duds tm) x
    (with-slots (cb basis accfun) sc
      (with-slots (mu) cb
	(let* ((rm (norme2 u))
	       (e (/ (- (* 2 (norme2 duds)) mu) rm))
	       (r (spin (first basis) u))
	       (v (* 2 (*g3 duds (first basis) (revg u))))
	       (f (funcall accfun s x sc)))
	  (make-instance
	   'spinorstate
	   :u duds
	   :duds (/ (+ (*g3 f r u) (* e u)) 2)
	   :tm (norme2 u)))))))

(defmethod to-cartesian ((x spinorstate) s sc)
  "Convert spinor state to cartesian coordinates given S and X"
  (with-slots (u duds tm) x
    (with-slots (basis) sc
      (let* ((r (norme2 u))
	     (dudt (/ duds r)))
	(values
	 (make-instance
	  'cartstate
	  :r (spin (first basis) u)
	  :v (graden (* 2 (*g3 dudt (first basis) (revg u))) 1))
	 tm)))))

(defmethod to-spinor ((x cartstate) tm sc)
  "Convert cartesian state to spinor given time and X"
  (with-slots (r v) x
    (with-slots (basis) sc
      (let ((u (recoverspinor (norme r) (rvbasis r v) basis)))
	(values
	 (make-instance
	  'spinorstate
	  :u u
	  :duds (* 0.5d0 (*g3 v u (first basis)))
	  :tm tm)
	 0)))))

(defmethod to-initial-spinor ((x cartstate) tm sc)
  "Convert to an initial spinor state from cartesian state & time"
  (with-slots (r v) x
    (with-slots (basis) sc
      (let ((u (recoverspinor (norme r) (rvbasis r v) basis)))
	(values
	 (make-instance
	  'spinorstate
	  :u u
	  :duds (/ (*g3 v u (first basis)) 2)
	  :tm tm)
	 0)))))

(defmethod to-spinor (x s sc)
  "Default: try converting to cartesian first"
  (multiple-value-bind (xc tm) (to-cartesian x s sc)
    (to-spinor xc tm sc)))

;;; Kustaanheimo-Stiefel Orbit Element equations of motion

(defclass ksstate ()
  ((alpha :initarg :alpha :documentation "U at s=0")
   (beta :initarg :beta :documentation "dU/ds / w0 at s=0")
   (e :initarg :e :documentation "Keplerian specific orbital energy")
   (tm :initarg :tm :documentation "Time"))
  (:documentation "Kustaanheimo-Stiefel orbital element state"))

(let ((slots '(alpha beta e tm)))
  (defmethod print-object ((x ksstate) s)
    (format s "#<KSSTATE ~{~a~^ ~}>"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value x slot))))))

(defstatearithmetic ksstate (alpha beta e tm))

(defmethod eom (s (x ksstate) &optional sc)
  "Kustaanheimo-Stiefel orbit element equations of motion from Arakida and Fukushima"
  (with-slots (alpha beta e tm) x
    (with-slots (basis cb accfun x0) sc
      (with-slots (mu) cb
	(with-slots ((e0 e)) x0
	  (let* ((hk0 (- e0))
		 (w0 (sqrt (/ hk0 2)))
		 (hk (- e))
		 (w (/ (- hk hk0) 2))
		 (u (+ (* alpha (cos (* w0 s)))
		       (* beta (sin (* w0 s)))))
		 (duds (* w0
			  (- (* beta (cos (* w0 s)))
			     (* alpha (sin (* w0 s))))))
		 (rv (spin (first basis) u))
		 (rm (norme rv))
		 (vv (* 2 (/ rm) (*g3 duds (first basis) (revg u))))
		 (fv (funcall accfun s x sc))
		 (ff (- (/ (*g fv rv u) 2) (* w u))))
	    (make-instance
	     'ksstate
	     :alpha (- (* ff (/ (sin (* w0 s)) w0)))
	     :beta (* ff (/ (cos (* w0 s)) w0))
	     :e (* rm (scalar (*i vv fv)))
	     :tm rm)))))))

(defmethod to-spinor ((x ksstate) s sc)
  "Convert KS state to spinor state"
  (with-slots (alpha beta e tm) x
    (with-slots (x0) sc
      (let* ((e0 (energy x0 sc))
	     (hk0 (- e0))
	     (w0 (sqrt (/ hk0 2)))
	     (hk (- e)))
	(values
	 (make-instance
	  'spinorstate
	  :u (+ (* alpha (cos (* w0 s)))
		(* beta (sin (* w0 s))))
	  :duds (* w0
		   (- (* beta (cos (* w0 s)))
		      (* alpha (sin (* w0 s)))))
	  :tm tm)
	 s)))))

(defmethod to-cartesian ((x ksstate) s sc)
  "Convert KS state to cartesian"
  (to-cartesian (to-spinor x s sc) s sc))

(defmethod energy ((x ksstate) sc)
  (slot-value x 'e))

(defmethod time-of (s (x ksstate))
  (slot-value x 'tm))

(defmethod to-ks ((x spinorstate) s sc)
  "Convert spinor state to KS"
  (with-slots (u duds tm) x
    (with-slots (x0) sc
      (let ((e0 (energy x0 sc)))
	(let* ((hk0 (- e0))
	       (w0 (sqrt (/ hk0 2)))
	       (e (energy x sc))
	       (hk (- e))
	       (alpha (- (* u (cos (* w0 s)))
			 (* (/ duds w0) (sin (* w0 s)))))
	       (beta (+ (* u (sin (* w0 s)))
			(* (/ duds w0) (cos (* w0 s))))))
	  (make-instance
	   'ksstate
	   :alpha alpha
	   :beta beta
	   :e e
	   :tm tm))))))

(defmethod to-ks ((x cartstate) tm sc)
  "Convert cartesian state to KS"
  (to-ks (to-spinor x 0 sc) tm sc))

(defmethod to-initial-ks ((x spinorstate) s sc)
  "Generate initial KS state from spinor state, s (assumed 0), and SC"
  (with-slots (u duds tm) x
    (let* ((e0 (energy x sc))
	   (hk0 (- e0))
	   (w0 (sqrt (/ hk0 2)))
	   (alpha0 u)
	   (beta0 (/ duds w0)))
      (values
       (make-instance
	'ksstate
	:alpha alpha0
	:beta beta0
	:e e0
	:tm tm)
       0))))

(defmethod to-initial-ks ((x cartstate) tm sc)
  (to-initial-ks (to-initial-spinor x tm sc) 0 sc))  

(defmethod to-ks (x s sc)
  "Default: try cartesian first"
  (multiple-value-bind (xc tm) (to-cartesian x s sc)
    (to-ks xc tm sc)))

;;; Export functions

(defun write-cart-traj (file trajdata)
  "Write cartesian orbit data to format plotable by Gnuplot.
Columns:
1: Time
2,3,4: Position
5,6,7: Velocity"
  (with-open-file (s file :direction :output :if-exists :supersede)
    (format s "# Time X Y Z VX VY VZ~%")
    (loop for (tm x) in trajdata
       do (with-slots (r v) x
	    (format s "~&~a"
		    (substitute #\E #\d
				(format nil "~a ~a ~a ~a ~a ~a ~a" 
					tm 
					(gref r :e1) (gref r :e2) (aif (gref r :e3) it 0)
					(gref v :e1) (gref v :e2) (aif (gref v :e3) it 0))))))))

(defun to-cartesian-traj (trajdata sc)
  "Convert trajectory data to cartesian"
  (loop for (s x) in trajdata
     collect (reverse (multiple-value-list (to-cartesian x s sc)))))

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
    (with-slots (cb basis) sc
      (with-slots (mu) cb
	(multiple-value-bind (r v) (coe2rv a e i lan aop f mu basis)
	  (values
	   (make-instance
	    'cartstate
	    :r r
	    :v v)
	   tm))))))

(defmethod to-coe ((x cartstate) tm sc)
  "Convert cartesian to classical orbital element state"
  (with-slots (r v) x
    (with-slots (cb basis) sc
      (with-slots (mu) cb
	(let* ((mombv (mombv-rv r v))
	       (nodev (nodev-rv mombv basis))
	       (eccv (eccv-rv r v mu)))
	  (values
	   (make-instance
	    'coestate
	    :a (sma-rv (norme r) (norme v) mu)
	    :e (norme eccv)
	    :i (inc-rv mombv basis)
	    :lan (lan-rv nodev basis)
	    :aop (aop-rv nodev eccv mombv)
	    :tm tm)
	   (truan-rv eccv r mombv)))))))

(defmethod energy ((x coestate) sc)
  (with-slots (a) x
    (with-slots (cb) sc
      (with-slots (mu) cb
	(- (/ mu 2 a))))))

(defmethod eom (f (x coestate) &optional sc)
  "Classical orbital element equations of motion. Watch for singularities."
  (with-slots (a e i lan aop tm) x
    (with-slots (cb accfun basis) sc
      (with-slots (mu) cb
	(multiple-value-bind (r v) (coe2rv a e i lan aop f mu basis) ; position & velocity
	  (let* ((rm (norme r)) ; radius
		 (p (* a (- 1 (expt e 2)))) ; semi latus rectum
		 (fv (funcall accfun f x sc)) ; force vector
		 (obasis (rvbasis r v)) ; orbit basis from position/velocity
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

;;; Modified equinoctial orbital elements

(defclass meestate ()
  ((p :initarg :p :documentation "Semi latus rectum")
   (f :initarg :f :documentation "e cos (AOP + LAN)")
   (g :initarg :g :documentation "e sin (AOP + LAN)")
   (h :initarg :h :documentation "tan (i/2) cos (LAN)")
   (k :initarg :k :documentation "tan (i/2) sin (LAN)")
   (l :initarg :l :documentation "True longitude: AOP + LAN + true anomaly"))
  (:documentation "Equinoctial orbital element state"))

(defstatearithmetic meestate (p f g h k l))

(let ((slots '(p f g h k l)))
  (defmethod print-object ((x meestate) s)
    (format s "#<MEESTATE ~{~a~^ ~}"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value x slot))))))

(defgeneric to-mee (x s sc))

(defmethod to-mee ((x coestate) truan sc)
  (with-slots (a e i lan aop lan tm) x
    (values
     (make-instance
      'meestate
      :p (* a (- 1 (expt e 2)))
      :f (* e (cos (+ aop lan)))
      :g (* e (sin (+ aop lan)))
      :h (* (tan (/ i 2)) (cos lan))
      :k (* (tan (/ i 2)) (sin lan))
      :l (+ lan aop truan))
     tm)))

(defmethod to-cartesian ((x meestate) tm sc)
  "Modified equinoctial element state to cartesian"
  (with-slots (p f g h k l) x
    (with-slots (cb basis) sc
      (with-slots (mu) cb
	(let* ((alpha2 (- (expt h 2) (expt k 2)))
	       (s2 (+ 1 (expt h 2) (expt k 2)))
	       (w (+ 1 (* f (cos l)) (* g (sin l))))
	       (r (/ p w)))
	  (values
	   (make-instance
	    'cartstate
	    :r (make-vector
		(list
		 (* (/ r s2)
		    (+ (cos l) (* alpha2 (cos l)) (* 2 h k (sin l))))
		 (* (/ r s2)
		    (+ (sin l) (- (* alpha2 (sin l))) (* 2 h k (cos l))))
		 (* (/ (* 2 r) s2)
		    (- (* h (sin l)) (* k (cos l)))))
		basis)
	    :v (make-vector
		(list
		 (* (/ -1 s2)
		    (sqrt (/ mu p))
		    (+ (sin l)
		       (* alpha2 (sin l))
		       (* -2 h k (cos l))
		       g
		       (* -2 f h k)
		       (* alpha2 g)))
		 (* (/ -1 s2)
		    (sqrt (/ mu p))
		    (+ (- (cos l))
		       (* alpha2 (cos l))
		       (* 2 h k (sin l))
		       (- f)
		       (* 2 g h k)
		       (* alpha2 f)))
		 (* (/ 2 s2)
		    (sqrt (/ mu p))
		    (+ (* h (cos l))
		       (* k (sin l))
		       (* f h)
		       (* g k))))
		basis))
	   tm))))))

(defmethod to-coe ((x meestate) tm sc)
  "Classical orbital element state from modified equinoctial"
  (with-slots (p f g h k l) x
    (values
     (make-instance
      'coestate
      :a (/ p (- 1 (expt f 2) (expt g 2)))
      :e (sqrt (+ (expt f 2) (expt g 2)))
      :i (atan (* 2 (sqrt (+ (expt h 2) (expt k 2))))
	       (- 1 (expt h 2) (expt k 2)))
      :aop (atan (- (* g h) (* f k))
		 (+ (* f h) (* g k)))
      :lan (atan k h))
     (- l lan aop))))


(defmethod eom (tm (x meestate) &optional sc)
  "Modified equinoctial orbital elements equations of motion"
  (with-slots (p f g h k l) x
    (with-slots (cb basis accfun) sc
      (with-slots (mu) cb
	(with-slots (r v) (to-cartesian x tm sc)
	  (let* ((rvbasis (rvbasis r v))
		 (fv (funcall accfun tm x sc))
		 (fr (scalar (*i fv (first rvbasis))))
		 (ft (scalar (*i fv (second basis))))
		 (fn (scalar (*i fv (third basis))))
		 (w (+ 1 (* f (cos l)) (* g (sin l))))
		 (s2 (+ 1 (expt h 2) (expt k 2))))
	    (make-instance
	     'meestate
	     :p (* 2 p (/ w) (sqrt (/ p mu)) ft)
	     :f (* (sqrt (/ p mu))
		   (+ (* fr (sin l))
		      (* (+ (* (+ w 1) (cos l)) f) (/ ft w))
		      (* (- (* k (cos l)) (* h (sin l))) g fn (/ w))))
	     :g (* (sqrt (/ p mu))
		   (+ (- (* fr (cos l)))
		      (* (+ (* (+ w 1) (sin l)) g) (/ ft w))
		      (- (* (- (* h (sin l)) (* k (cos l))) (/ (* g fn) w)))))
	     :h (* (sqrt (/ p mu))
		   s2 fn (cos l) (/ (* 2 w)))
	     :k (* (sqrt (/ p mu))
		   s2 fn (sin l) (/ (* 2 w)))
	     :l (+ (* (sqrt (* mu p)) (expt (/ w p) 2))
		   (* (/ w) (sqrt (/ p mu)) (- (* h (sin l)) (* k (cos l))) fn))
	     )))))))

