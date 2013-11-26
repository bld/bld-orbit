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

(defparameter *au* 149597870691d0 "Astronomical unit (m)")

(defparameter *c* 299792458d0 "Speed of light (m/s)")

(defvar *lsun* 1361.7913d0 "Solar luminosity at 1 AU (W/m^2)")

(defparameter *J2000* (list (ve3 :e1 1d0)
			    (ve3 :e2 1d0)
			    (ve3 :e3 1d0)) "J2000 inertial frame")

(defparameter *musun* 1.32712440018d20 "Solar gravitational parameter, m^3/s^2")

(defun new-frame (rotor frame)
  "New reference frame from rotor & original frame"
  (mapcar #'(lambda (basisvector) (rotateg basisvector rotor)) frame))

;; Define infinity norm method for BLD-ODE Runge-Kutta method on geometric algebra objects
(defmethod norminfx ((x g))
  (norminf x))

(defgeneric energy (x sc) (:documentation "Keplerian orbit energy"))

(defgeneric eom (s x &optional param) (:documentation "Equations of motion given independent variable S & state X"))

(defgeneric to-cartesian (s x &optional sc) (:documentation "Convert a state to cartesian state"))

(defgeneric to-spinor (s x &optional sc) (:documentation "Convert a state to spinor state"))

(defgeneric to-ks (s x &optional sc) (:documentation "Convert to a Kustaanheimo-Stiefel orbit element state"))

;; Rotor, spinor and basis functions

(defun recoverrotor3d (fs es)
  "Recover a basis given new and original basis vectors"
  (let* ((esr (apply #'recipbvs es))
	 (psi (mapcar #'(lambda (f er) (+ (*g f er) 1)) fs esr)))
    (unitg (first psi))))

(defun recoverspinor3d (r fs es)
  "Recover a spinor given orbit radius, new basis vectors, and original basis vectors"
  (* (recoverrotor3d fs es) (sqrt r)))

(defun rvbasis (rv vv)
  "Return a set of basis vectors derived from position and velocity"
  (let* ((mombv (*o rv vv))
	 (x (unitg rv))
	 (y (unitg (*i rv mombv)))
	 (z (when (= 3 (dimension rv) (dimension vv))
	      (*x2 x y))))
    (if (= 2 (dimension rv) (dimension vv)) ; 2D or 3D?
	(list x y)
	(list x y z))))

;; Acceleration functions

(defun gravity (s r mu)
  "Gravitational acceleration"
  (- (* (/ mu (norme2 r)) (unitg r))))

(defun sailidealacc (s r v sc)
  "Sail acceleration with fixed orientation relative to sail position frame"
  (with-slots (lightness mu pointfun) sc
    (let ((n (funcall pointfun s r v sc)))
      (* lightness mu (/ (norme2 r))
	 (expt (scalar (*i (unitg r) n)) 2)
	 n))))

;; Pointing functions

(defun sailpointingfixed (s r v sc)
  "Sail pointed at fixed orientation with respect to position/velocity frame defined by rotor RS"
  (with-slots (rs basis) sc
    (let* ((rvbasis (rvbasis r v))
	   (rrv (recoverrotor3d rvbasis basis)))
      (rotateg (first basis) (*g rrv rs)))))

(defun sailpointingnormal (s r v sc)
  "Sail normal to sunlight"
  (unitg r))

(defun sailpointingtable (s r v sc)
  "Return sail normal vector"
  (with-slots (rs basis) sc
    (let* ((rvbasis (rvbasis r v))
	   (rrv (recoverrotor3d rvbasis basis))
	   (rsi (second (find s rs :test #'<= :key #'first))))
      (rotateg (first basis) (*g rrv rsi)))))

;; Sail class

(defclass sail ()
  ((eom :initarg :eom :initform #'eom :documentation "Equations of motion")
   (mu :initarg :mu :initform 1d0 :documentation "Gravitational parameter of central body at position (0 0 0) wrt spacecraft")
   (accfun :initarg :accfun :initform #'sailidealacc :documentation "Sail acceleration function")
   (pointfun :initarg :pointfun :initform #'sailpointingfixed :documentation "Sail pointing function")
   (lightness :initarg :lightness :initform 0d0 :documentation "Ratio of solar to gravitational acceleration")
   (basis :initarg :basis :initform (list (ve2 :e1 1) (ve2 :e2 1)))
   (t0 :initarg :t0 :initform 0d0 :documentation "Initial time")
   (tf :initarg :tf :initform (* 2 pi) :documentation "Final time")
   (x0 :initarg :x0 :initform (make-instance 'cartstate :r (ve2 :e1 1) :v (ve2 :e2 1)))
   (rs :initarg :rs :initform (re2 :s 1d0) :documentation "Sail orientation rotor wrt orbital position frame"))
  (:documentation "Solar sail orbit problem"))

(let ((slots '(eom mu accfun pointfun lightness basis t0 tf x0 rs)))
  (defmethod print-object ((sail sail) s)
    (format s "#<SAIL ~{~a~^~%  ~}>"
	    (loop for slot in slots
	       collect (format nil "~a ~a" slot (slot-value sail slot))))))

;; Propagate a trajectory
(defmethod propagate ((sc sail))
  "Propagate sailcraft trajectory"
  (with-slots (eom t0 tf x0) sc
    (rka eom t0 tf x0 :param sc :hmax (- tf t0))))

;; Cartesian state

(defclass cartstate ()
  ((r :initarg :r :documentation "Position vector")
   (v :initarg :v :documentation "Velocity vector"))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a>" (slot-value x 'r) (slot-value x 'v)))

(defstatearithmetic cartstate (r v))

(defmethod energy ((x cartstate) sc)
  (with-slots (mu) sc
    (with-slots (r v) x
      (- (/ (scalar (exptg v 2)) 2) (/ mu (norme r))))))

(defmethod eom (tm (x cartstate) &optional sc)
  "Solar sail cartesian equations of motion"
  (with-slots (r v) x
    (with-slots (lightness mu accfun) sc
      (make-instance
       'cartstate
       :r v
       :v (+ (funcall accfun tm r v sc)
	     (gravity tm r mu))))))

(defmethod to-cartesian (tm (x cartstate) &optional sc)
  (values tm x))

(defparameter *sail-2d-cart-kepler-eg*
  (make-instance 
   'sail
   :tf (* 4 pi)
   :x0 (make-instance
	'cartstate
	:r (ve2 :e1 1)
	:v (ve2 :e2 1.1))))

(defparameter *sail-2d-cart-normal-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :pointfun #'sailpointingnormal
   :lightness 0.1d0))

(defparameter *sail-2d-cart-fixed-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :lightness 0.1
   :rs (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2))))))

(defparameter *sail-2d-cart-table-eg*
  (make-instance
   'sail
   :lightness 0.1
   :tf (* 4 pi)
   :pointfun #'sailpointingtable
   :rs (list (list (* 2 pi) (rotor (bve2 :e1e2 1) (/ pi 2)))
	     (list (* 4 pi) (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2))))))))

;; Spinor states

(defclass spinorstate ()
  ((u :initarg :u :documentation "Spinor of position: r = u e1 (revg u)")
   (duds :initarg :duds :documentation "Spinor derivative wrt s")
   (tm :initarg :tm :documentation "Time"))
  (:documentation "Spinor state"))

(defmethod print-object ((x spinorstate) stream)
  (format stream "#<SPINORSTATE :U ~a :DUDS ~a :TM ~a>" (slot-value x 'u) (slot-value x 'duds) (slot-value x 'tm)))

(defstatearithmetic spinorstate (u duds tm))

(defmethod energy ((x spinorstate) sc)
  "Keplerian orbit energy"
  (with-slots (u duds) x
    (with-slots (mu) sc
      (/ (- (* 2 (norme2 duds)) mu) (norme2 u)))))

(defmethod eom (s (x spinorstate) &optional sc)
  "Spinor equations of motion: 2 dU/ds - E U = f r U, dt/ds = |r|"
  (with-slots (u duds tm) x
    (with-slots (mu basis accfun) sc
      (let* ((rm (norme2 u))
	     (e (/ (- (* 2 (norme2 duds)) mu) rm))
	     (r (spin (first basis) u))
	     (v (* 2 (*g3 duds (first basis) (revg u))))
	     (f (funcall accfun s r v sc)))
      (make-instance
       'spinorstate
       :u duds
       :duds (/ (+ (*g3 f r u) (* e u)) 2)
       :tm (norme2 u))))))

(defmethod to-cartesian (s (x spinorstate) &optional sc)
  "Convert spinor state to cartesian coordinates given S and X"
  (with-slots (u duds tm) x
    (with-slots (basis) sc
      (let* ((r (norme2 u))
	     (dudt (/ duds r)))
	(values
	 tm
	 (make-instance
	  'cartstate
	  :r (spin (first basis) u)
	  :v (* 2 (*g3 dudt (first basis) (revg u)))))))))

(defmethod to-spinor (tm (x cartstate) &optional sc)
  "Convert cartesian state to spinor given time and X"
  (with-slots (r v) x
    (with-slots (basis) sc
      (let ((u (recoverspinor3d (norme r) (rvbasis r v) basis)))
	(make-instance
	 'spinorstate
	 :u u
	 :duds (* 0.5d0 (*g3 v u (first basis)))
	 :tm tm)))))

(defparameter *sail-2d-spin-kepler-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :accfun #'(lambda (s r v sc) (ve2))
   :lightness 0.0d0
   :x0 (to-spinor 0 (slot-value *sail-2d-cart-kepler-eg* 'x0) *sail-2d-cart-kepler-eg*)))

(defparameter *sail-2d-spin-normal-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :pointfun #'sailpointingnormal
   :lightness 0.1d0
   :x0 (to-spinor 0 
		  (slot-value *sail-2d-cart-normal-eg* 'x0)
		  *sail-2d-cart-normal-eg*)))

(defparameter *sail-2d-spin-fixed-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :pointfun #'sailpointingfixed
   :lightness 0.1d0
   :x0 (to-spinor 0 (slot-value *sail-2d-cart-fixed-eg* 'x0) *sail-2d-cart-fixed-eg*)
   :rs (slot-value *sail-2d-cart-fixed-eg* 'rs))
  "Two-dimensional spinor EOM with fixed pointing sail example")

(defparameter *sail-2d-spin-table-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :pointfun #'sailpointingtable
   :lightness 0.1d0
   :x0 (to-spinor 0 (slot-value *sail-2d-cart-fixed-eg* 'x0) *sail-2d-cart-fixed-eg*)
   :rs (slot-value *sail-2d-cart-table-eg* 'rs))
  "Two-dimensional spinor EOM with lookup table pointing sail example")

;; Kustaanheimo-Stiefel Orbit Element equations of motion

(defclass ksstate ()
  ((alpha :initarg :alpha :documentation "U at s=0")
   (beta :initarg :beta :documentation "dU/ds / w0 at s=0")
   (e :initarg :e :documentation "Keplerian orbital energy")
   (tm :initarg :tm :documentation "Time"))
  (:documentation "Kustaanheimo-Stiefel orbital element state"))

(let ((slots '(alpha beta e tm)))
  (defmethod print-object ((x ksstate) s)
    (format s "#<KSSTATE ~{~a~^ ~}>"
	    (loop for slot in slots
	       collect slot
	       collect (slot-value x slot)))))

(defstatearithmetic ksstate (alpha beta e tm))

(defmethod eom (s (x ksstate) &optional sc)
  "Kustaanheimo-Stiefel orbit element equations of motion from Arakida and Fukushima"
  (with-slots (alpha beta e tm) x
    (with-slots (basis mu accfun x0) sc
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
	       (fv (funcall accfun s rv vv sc))
	       (ff (- (/ (*g fv rv u) 2) (* w u))))
	  (make-instance
	   'ksstate
	   :alpha (- (* ff (/ (sin (* w0 s)) w0)))
	   :beta (* ff (/ (cos (* w0 s)) w0))
	   :e (* rm (scalar (*i vv fv)))
	   :tm rm))))))

(defmethod to-spinor (s (x ksstate) &optional sc)
  "Convert KS state to spinor state"
  (with-slots (alpha beta e tm) x
    (with-slots (mu x0) sc
      (let* ((e0 (energy x0 sc))
	     (hk0 (- e0))
	     (w0 (sqrt (/ hk0 2)))
	     (hk (- e)))
	(make-instance
	 'spinorstate
	 :u (+ (* alpha (cos (* w0 s)))
	       (* beta (sin (* w0 s))))
	 :duds (* w0
		  (- (* beta (cos (* w0 s)))
		     (* alpha (sin (* w0 s)))))
	 :tm tm)))))

(defmethod to-cartesian (s (x ksstate) &optional sc)
  "Convert KS state to cartesian"
  (to-cartesian s (to-spinor s x sc) sc))

(defmethod energy ((x ksstate) sc)
  (slot-value x 'e))

(defmethod to-ks (s (x spinorstate) &optional sc)
  "Convert spinor state to KS"
  (with-slots (u duds tm) x
    (with-slots (mu x0) sc
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

(defmethod to-ks (tm (x cartstate) &optional sc)
  "Convert cartesian state to KS"
  (to-ks 0 (to-spinor 0 x sc) sc))

(defparameter *sail-2d-ks-kepler-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :accfun #'(lambda (s r v sc) (ve2))
   :lightness 0d0
   :x0 (to-ks 0 (slot-value *sail-2d-spin-kepler-eg* 'x0) *sail-2d-spin-kepler-eg*))
  "Two-dimensional Kustaanheimo-Stiefel EOM Keplerian (unperturbed) orbit example")

(defparameter *sail-2d-ks-normal-eg*
  (make-instance
   'sail
   :tf (* 8 pi)
   :pointfun #'sailpointingnormal
   :lightness 0.1d0
   :x0 (to-ks 0 (slot-value *sail-2d-cart-normal-eg* 'x0) *sail-2d-cart-normal-eg*))
  "Two-dimensional Kustaanheimo-Stiefel EOM with sun-pointing sail example")
  
(defparameter *sail-2d-ks-fixed-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :pointfun #'sailpointingfixed
   :lightness 0.1d0
   :x0 (to-ks 0 (slot-value *sail-2d-cart-fixed-eg* 'x0) *sail-2d-cart-fixed-eg*)
   :rs (slot-value *sail-2d-cart-fixed-eg* 'rs))
  "Two-dimensional Kustaanheimo-Stiefel EOM with fixed pointing sail example")

(defparameter *sail-2d-ks-table-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :pointfun #'sailpointingtable
   :lightness 0.1d0
   :x0 (to-ks 0 (slot-value *sail-2d-cart-fixed-eg* 'x0) *sail-2d-cart-fixed-eg*)
   :rs (slot-value *sail-2d-cart-table-eg* 'rs))
  "Two-dimensional Kustaanheimo-Stiefel EOM with fixed pointing sail example")

;; Export functions

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
     collect (multiple-value-list (to-cartesian s x sc))))

;;; Classical orbital element conversion

;; position & velocity to classical orbit elements
(defun energy-rv (r v mu)
  "Orbit energy from position, velocity, & gravitational parameter"
  (- (/ (* v v) 2)
     (/ mu r)))
(defun sma-rv (r v mu)
  "semi-major axis from position, velocity, & gravitational parameter"
  (/ mu (* -2 (energy-rv r v mu))))
(defun eccv-rv (rv vv mu)
  "eccentricity vector from position, velocity, & gravitational parameter"
  (/ (- (* rv (- (norme2 vv) (/ mu (norme rv))))
	(* vv (scalar (*i rv vv))))
       mu))
(defun mombv-rv (rv vv)
  "orbital momentum bivector from position & velocity"
  (*o rv vv))
(defun nodev-rv (mombv basis)
  "ascending node vector from momentum and basis"
  (- (*i (third basis) mombv)))
(defun inc-rv (mombv basis)
  "inclination from momentum bivector and basis"
  (acos (scalar (dual (*o (third basis) (unitg mombv))))))
(defun lan-rv (nodev basis)
  "Longitude of ascending node from node vector and basis"
  (let ((tmp (atan (scalar (*i nodev (second basis)))
		   (scalar (*i nodev (first basis))))))
    (if (< tmp 0) (+ tmp (* 2 pi)) tmp)))	 
(defun aop-rv (nodev eccv mombv)
  "argument of perigee from node vector, eccentricity vector, and momentum bivector"
  (let* ((n (unitg nodev))
	 (e (unitg eccv))
	 (h (unitg mombv))
	 (ov (*i n h))
	 (tmp (atan (scalar (*i e ov)) (scalar (*i e n)))))
    (if (< tmp 0) (+ tmp (* 2 pi)) tmp)))
(defun truan-rv (eccv rv mombv)
  "true anomaly from eccentricity, position, and velocity vectors"
  (let* ((e (unitg eccv))
	 (h (unitg mombv))
	 (ruv (unitg rv))
	 (ov (*i e h))
	 (tmp (atan (scalar (*i ruv ov))
		    (scalar (*i ruv e)))))
    (if (< tmp 0) (+ tmp (* 2 pi)) tmp)))
(defun rv2coe (rv vv mu basis)
  "Return classical orbital elements given position & velocity vectors, gravitational parameter, and a list of 3 basis vectors"
  (let* ((mombv (mombv-rv rv vv))
	 (nodev (nodev-rv mombv basis))
	 (eccv (eccv-rv rv vv mu)))
    (make-hash
     :sma (sma-rv (norme rv) (norme vv) mu)
     :ecc (norme eccv)
     :inc (inc-rv mombv basis)
     :lan (raan-rv nodev basis)
     :aop (aop-rv nodev eccv mombv)
     :truan (truan-rv eccv rv mombv))))

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
	       collect slot
	       collect (slot-value x slot)))))

(defmethod to-cartesian (f (x coestate) &optional sc)
  (with-slots (a e i lan aop tm) x
    (with-slots (mu basis) sc
      (multiple-value-bind (r v) (coe2rv a e i lan aop f mu basis)
	(values
	 tm
	 (make-instance
	  'cartstate
	  :r r
	  :v v))))))

(defmethod to-coe (tm (x cartstate) &optional sc)
  "Convert cartesian to classical orbital element state"
  (with-slots (r v) x
    (with-slots (mu basis) sc
      (let* ((mombv (mombv-rv r v))
	     (nodev (nodev-rv mombv basis))
	     (eccv (eccv-rv r v mu)))
	(values
	 (truan-rv eccv r mombv)
	 (make-instance
	  'coestate
	  :a (sma-rv (norme r) (norme v) mu)
	  :e (norme eccv)
	  :i (inc-rv mombv basis)
	  :lan (lan-rv nodev basis)
	  :aop (aop-rv nodev eccv mombv)
	  :tm tm))))))

(defmethod energy ((x coestate) sc)
  (with-slots (a) x
    (with-slots (mu) sc
      (- (/ mu 2 a)))))

(defmethod eom (f (x coestate) &optional sc)
  "Classical orbital element equations of motion. Watch for singularities."
  (with-slots (a e i lan aop tm) x
    (with-slots (mu accfun basis) sc
      (multiple-value-bind (r v) (coe2rv a e i lan aop f mu basis) ; position & velocity
	(let* ((rm (norme r)) ; radius
	       (p (* a (- 1 (expt e 2)))) ; semi latus rectum
	       (fv (funcall accfun f r v sc)) ; force vector
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
			 (* tf (+ 1 (/ rm p)) (sin f))))))))))))

;;; 3D examples

(defparameter *sail-3d-coe-kepler-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
   :accfun #'(lambda (s r v sc) (ve3))
   :lightness 0d0
   :x0 (make-instance
	'coestate
	:a 1.1
	:e 0.1
	:i 0.1
	:lan 0.1
	:aop 0.1
	:tm 0)
   :basis (list (ve3 :e1 1) (ve3 :e2 1) (ve3 :e3 1))
   :rs (re3 :s 1)))

(defparameter *sail-3d-cart-kepler-eg*
  (let ((sc *sail-3d-coe-kepler-eg*))
    (with-slots (t0 x0 basis rs) sc
      (make-instance
       'sail
       :t0 0
       :tf (* 4 pi)
       :accfun #'(lambda (s r v sc) (ve3))
       :lightness 0d0
       :x0 (second (multiple-value-list (to-cartesian t0 x0 sc)))
       :basis basis
       :rs rs))))

(defparameter *sail-3d-spin-kepler-eg*
  (let ((sc *sail-3d-cart-kepler-eg*))
    (with-slots (x0 basis rs) sc
      (make-instance
       'sail
       :t0 0
       :tf (* 4 pi)
       :accfun #'(lambda (s r v sc) (ve3))
       :lightness 0d0
       :x0 (to-spinor 0 x0 sc)
       :basis basis
       :rs rs))))

(defparameter *sail-3d-ks-kepler-eg*
  (let ((sc *sail-3d-cart-kepler-eg*))
    (with-slots (x0 basis rs) sc
      (make-instance
       'sail
       :t0 0
       :tf (* 4 pi)
       :accfun #'(lambda (s r v sc) (ve3))
       :lightness 0d0
       :x0 (to-ks 0 x0 sc)
       :basis basis
       :rs rs))))
