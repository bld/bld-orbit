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

;; Define infinity norm method for BLD-ODE Runge-Kutta method on geometric algebra objects
(defmethod norminfx ((x g))
  (norminf x))

;; Variable to access problem specific data
(defvar *sc*)

(defgeneric eom (s x) (:documentation "Equations of motion given independent variable S & state X"))

(defgeneric to-cartesian (s x) (:documentation "Convert a state to cartesian state"))

;; Acceleration functions

(defun gravity (tm x)
  "Gravitational acceleration"
  (let ((xcart (to-cartesian tm x)))
    (with-slots (r) xcart
      (with-slots (mu) *sc*
	(- (* (/ mu (norme2 r)) (unitg r)))))))

(defun sailidealacc (s x)
  "Sail acceleration with fixed orientation relative to sail position frame"
  (let ((xcart (to-cartesian s x)))
    (with-slots (r) xcart
      (with-slots (lightness mu pointfun) *sc*
	(let ((n (funcall pointfun s x)))
	  (* lightness mu (/ (norme2 r))
	     (expt (scalar (*i (unitg r) n)) 2)
	     n))))))

;; Pointing functions

(defun sailpointingfixed (s x)
  "Return sail normal vector"
  (with-slots (rs basis) *sc*
    (let ((xcart (to-cartesian s x)))
      (with-slots (r v) xcart
	(let* ((rvbasis (rvbasis r v))
	       (rrv (recoverrotor3d rvbasis basis)))
	  (rot (first basis) (*g rrv rs)))))))

(defun sailpointingnormal (tm x)
  "Sail normal to sunlight"
  (let ((xcart (to-cartesian tm x)))
    (unitg (slot-value xcart 'r))))

(defun sailpointingtable (tm x)
  "Sail pointing from lookup table"
  (with-slots (r) x
    (with-slots (rs) *sc*
      (rot (unitg r) 
	   (second (find tm rs :test #'<= :key #'first))))))

;; Sail class
(defclass sail ()
  ((eom :initarg :eom :initform #'eom :documentation "Equations of motion")
   (mu :initarg :mu :initform 1d0 :documentation "Gravitational parameter of central body at position (0 0 0) wrt spacecraft")
   (accfun :initarg :accfun :initform #'sailidealacc :documentation "Sail acceleration function")
   (pointfun :initarg :pointfun :initform #'sailpointingfixed :documentation "Sail pointing function")
   (lightness :initarg :lightness :initform 0d0 :documentation "Ratio of solar to gravitational acceleration")
   (basis :initarg :basis :initform (list (ve2 :c1 1) (ve2 :c10 1)))
   (t0 :initarg :t0 :initform 0d0 :documentation "Initial time")
   (tf :initarg :tf :initform (* 2 pi) :documentation "Final time")
   (x0 :initarg :x0 :initform (make-instance 'cartstate :r (ve2 :c1 1) :v (ve2 :c10 1)))
   (rs :initarg :rs :initform (re2 :c0 1d0) :documentation "Sail orientation rotor wrt orbital position frame")))

;; Propagate a trajectory
(defmethod propagate ((sc sail))
  "Propagate sailcraft trajectory"
  (let ((*sc* sc)) ; bind *sc* variable for use by equations of motion and acceleration functions
    (with-slots (eom t0 tf x0) sc
      (rka eom t0 tf x0))))

;; Cartesian state

(defclass cartstate ()
  ((r :initarg :r :documentation "Position vector")
   (v :initarg :v :documentation "Velocity vector"))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a>" (slot-value x 'r) (slot-value x 'v)))

(defstatearithmetic cartstate (r v))

(defmethod eom (tm (x cartstate))
  "Solar sail cartesian equations of motion"
  (with-slots (r v) x
    (with-slots (lightness mu accfun) *sc*
      (make-instance
       'cartstate
       :r v
       :v (+ (funcall accfun tm x)
	     (gravity tm x))))))

(defparameter *sail-2d-cart-kepler-eg*
  (make-instance 'sail))

(defparameter *sail-2d-cart-normal-eg*
  (make-instance
   'sail
   :pointfun #'sailpointingnormal
   :lightness 0.1d0))

(defparameter *sail-2d-cart-fixed-eg*
  (make-instance
   'sail
   :lightness 0.1
   :rs (rotor (bve2 :c11 1) (* 35.5 (/ pi 180)))))

(defparameter *sail-2d-cart-table-eg*
  (make-instance
   'sail
   :lightness 0.1
   :tf (* 4 pi)
   :pointfun #'sailpointingtable
   :rs (list (list (* 2 pi) (rotor (bve2 :c11 1) (/ pi 2)))
	     (list (* 4 pi) (rotor (bve2 :c11 1) (/ pi 2 3))))))

;; Kustaanheimo-Stiefel spinor states

(defclass spinorstate ()
  ((u :initarg :u :documentation "Spinor of position: r = u e1 (revg u)")
   (duds :initarg :duds :documentation "Spinor derivative wrt s")
   (tm :initarg :tm :documentation "Time"))
  (:documentation "Spinor state"))

(defmethod print-object ((x spinorstate) stream)
  (format stream "#<SPINORSTATE :U ~a :DUDS ~a :TM ~a>" (slot-value x 'u) (slot-value x 'duds) (slot-value x 'tm)))

(defstatearithmetic spinorstate (u duds tm))

(defmethod eom (s (x spinorstate))
  "Spinor equations of motion: 2 dU/ds - E U = f r U, dt/ds = |r|"
  (with-slots (u duds tm) x
    (with-slots (mu basis accfun) *sc*
      (let* ((rm (norme2 u))
	     (e (/ (- (* 2 (norme2 duds)) mu) rm))
	     (r (spin (first basis) u))
	     (f (funcall accfun s x)))
      (make-instance
       'spinorstate
       :u duds
       :duds (/ (+ (*g3 f r u) (* e u)) 2)
       :tm (norme2 u))))))

(defparameter *sail-2d-spin-kepler-eg*
  (make-instance
   'sail
   :accfun #'(lambda (s x) (ve2))
   :lightness 0.0d0
   :x0 (make-instance 
	'spinorstate
	:u (re2 :c0 1d0)
	:duds (re2 :c11 -0.5d0)
	:tm 0)))

(defparameter *sail-2d-spin-normal-eg*
  (make-instance
   'sail
   :pointfun #'sailpointingnormal
   :lightness 0.1d0
   :x0 (make-instance 
	'spinorstate
	:u (re2 :c0 1d0)
	:duds (re2 :c11 -0.5d0)
	:tm 0)))

(defparameter *sail-2d-spin-fixed-eg*
  (make-instance
   'sail
   :lightness 0.1d0
   :x0 (make-instance 
	'spinorstate
	:u (re2 :c0 1d0)
	:duds (re2 :c11 -0.5d0)
	:tm 0)
   :rs (rotor (bve2 :c11 1d0) (* 35.5 (/ pi 180d0)))))

;; Kustaanheimo-Stiefel Orbit Element equations of motion
;; WORK IN PROGRESS - DOESN'T WORK YET

(defclass ksstate ()
  ((tm :initarg :tm :documentation "Time")
   (alpha :initarg :alpha :documentation "U_0")
   (beta :initarg :beta :documentation "dU/ds_0 / omega")
   (e :initarg :e :documentation "Keplerian orbital energy"))
  (:documentation "Kustaanheimo-Stiefel Orbit Element State"))

(defstatearithmetic ksstate (tm alpha beta e))

(defmethod print-object ((x ksstate) stream)
  (with-slots (tm alpha beta e) x
    (format stream "#<KSSTATE :TM ~a :ALPHA ~a :BETA ~a :E ~a>" tm alpha beta e)))

(defmethod eom (s (x ksstate))
  "Kustaanheimo-Stiefel Keplerian orbit element equations of motion"
  (with-slots (tm alpha beta e) x
    (with-slots (basis accfun) *sc*
      (let* ((w0 (sqrt (/ e -2d0)))
	     (u (+ (* alpha (cos (* w0 s))) (* beta (sin (* w0 s)))))
	     (up (* w0 (- (* beta (cos (* w0 s))) (* alpha (sin (* w0 s))))))
	     (r (spin (first basis) u))
	     (rm (norme2 u))
	     (v (* 2 (/ rm) (*g3 up (first basis) (revg u))))
	     (f (funcall accfun s x))
	     (ff (/ (*g3 f r u) 2d0)))
	(make-instance
	 'ksstate
	 :tm rm
	 :alpha (* -1 ff (/ w0) (sin (* w0 s)))
	 :beta (* ff (/ w0) (cos (* w0 s)))
	 :e (* rm (scalar (*i f v))))))))

(defparameter *sail-2d-ks-kepler-eg*
  (make-instance
   'sail
   :accfun #'(lambda (s x) (ve2))
   :lightness 0.0d0
   :x0 (make-instance
	'ksstate
	:tm 0d0
	:alpha (re2 :c0 1d0)
	:beta (re2 :c11 -1d0)
	:e -0.5d0))
  "Unperturbed Kepler problem")

(defparameter *sail-2d-ks-normal-eg*
  (make-instance
   'sail
   :pointfun #'sailpointingnormal
   :lightness 0.1d0
   :x0 (make-instance
	'ksstate
	:tm 0d0
	:alpha (re2 :c0 1d0)
	:beta (re2 :c11 -1d0)
	:e -0.5d0))
  "Sail pointed normal to sunlight")

(defparameter *sail-2d-ks-fixed-eg*
  (make-instance
   'sail
   :lightness 0.1d0
   :x0 (make-instance
	'ksstate
	:tm 0d0
	:alpha (re2 :c0 1d0)
	:beta (re2 :c11 -1d0)
	:e -0.5d0)
   :rs (rotor (bve2 :c11 1d0) (* 35.5 (/ pi 180d0))))
  "Sail pointed at fixed orientation to sunlight")

;; State conversion functions

(defgeneric to-spinor (s x) (:documentation "Convert a state to spinor state"))

(defgeneric to-ks (s x &key mu basis) (:documentation "Convert to a Kustaanheimo-Stiefel state"))

(defmethod to-cartesian (tm (x cartstate))
  (values x tm))

(defmethod to-cartesian (s (x spinorstate))
  "Convert spinor state to cartesian coordinates given S and X"
  (with-slots (u duds tm) x
    (with-slots (basis) *sc*
      (let* ((r (norme2 u))
	     (dudt (/ duds r)))
	(values 
	 (make-instance
	  'cartstate
	  :r (spin (first basis) u)
	  :v (* 2 (*g3 dudt (first basis) (revg u))))
	 tm)))))

(defmethod to-cartesian (s (x ksstate))
  "Convert Kustaanheimo-Stiefel orbit element state to cartesian"
  (with-slots (tm alpha beta e) x
    (with-slots (basis) *sc*
      (let* ((w0 (sqrt (/ e -2d0)))
	     (u (+ (* alpha (cos (* w0 s))) (* beta (sin (* w0 s)))))
	     (up (* w0 (- (* beta (cos (* w0 s))) (* alpha (sin (* w0 s)))))))
	(values
	 (make-instance
	  'cartstate
	  :r (spin (first basis) u)
	  :v (* 2d0 (/ (norme2 u)) (*g3 up (first basis) (revg u))))
	 tm)))))

(defmethod to-spinor (tm (x cartstate))
  "Convert cartesian state to spinor given time and X"
  (with-slots (r v) x
    (with-slots (basis) *sc*
      (let ((u (recoverspinor3d (norme r) (rvbasis r v) basis)))
	(make-instance
	 'spinorstate
	 :u u
	 :duds (* 0.5d0 (*g3 v u (first basis)))
	 :tm tm)))))

(defmethod to-ks (tm (x cartstate) &key mu basis)
  (with-slots (r v) x
    (let* ((rvbasis (rvbasis r v))
	   (u (recoverspinor3d (norme r) rvbasis basis))
	   (up (/ (*g3 v u (first basis)) 2d0))
	   (e (/ (- (* 2 (norme2 up)) mu) (norme2 u)))
	   (w0 (sqrt (/ e -2d0))))
      (make-instance
       'ksstate
       :tm tm
       :e e
       :alpha u
       :beta (/ up w0)))))

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
					(gref r #b1) (gref r #b10) (aif (gref r #b100) it 0)
					(gref v #b1) (gref v #b10) (aif (gref v #b100) it 0))))))))

(defun spinor-to-cartesian-traj (traj)
  "Convert results of PROPAGATE for spinor problem to cartesian trajectory"
  (loop for (s x) in traj
     collect (list (slot-value x 'tm) (to-cartesian s x))))

;; OLDER STUFF: NEED TO REWORK

;; Kustaanheimo-Stiefel-Hestenes (KSH) state & equations of motion

(defclass kshstate ()
  ((alpha :initarg :alpha :documentation "Initial orbit spinor")
   (beta :initarg :beta :documentation "Initial ")
   (e :initarg :e :documentation "Specific Kepler orbit energy")
   (tm :initarg :tm :documentation "Time")))

(defmethod print-object ((x kshstate) stream)
  (with-slots (alpha beta e tm) x
    (format stream "#<KSHSTATE :alpha ~a :beta ~a :e ~a :tm ~a>" alpha beta e tm)))

(defstatearithmetic kshstate (alpha beta e tm))

;; Sail force functions
(defvar *lightness* 0.1 "Sail lightness number")
(defvar *mu* 1 "Gravitational parameter")
(defmethod sail-normal-force ((rv g) mu lightness)
  "Solar sail force function when sail normal to the sun. RV is the position vector from the sun to the sail."
  (* (unitg rv) (/ (* lightness mu) (norme2 rv))))

;; Parameters used in equations of motion
(defvar *ksh-forcefun* #'(lambda (s x) (ve3)) "KSH force function")
(defvar *ksh-sigma0* (ve3 :c1 1) "KSH reference unit position vector")

;; Convenience functions for KSH EOM
(defun w0 (e)
  "Average orbit angular velocity given energy"
  (sqrt (- (/ e 2))))
(defmethod u ((alpha g) (beta g) w0 s)
  "Orbit spinor given ALPHA, BETA, W0, and S"
  (+ (* alpha (cos (* w0 s)))
     (* beta (sin (* w0 s)))))
(defmethod alpha-from-u-duds-w0-s ((u g) (duds g) w0 s)
  "Alpha (U0) given U, DUDS, W0, and S"
  (- (* u (cos (* w0 s)))
     (* duds (/ (sin (* w0 s))
		w0))))
(defmethod beta-from-u-duds-w0-s ((u g) (duds g) w0 s)
  "Beta (dU0/ds/w0) given U, DUDS, w0, and s"
  (+ (* u (sin (* w0 s)))
     (* duds (/ (cos (* w0 s))
		w0))))
(defmethod duds (alpha beta w0 s)
  "s-derivative of spinor given BETA and W0"
  (* (- (* beta (cos (* w0 s)))
	(* alpha (sin (* w0 s))))
     w0))
(defmethod dalphads ((ff g) w0 s)
  "s-derivative of alpha given FF, W0, and S"
  (* ff (- (/ (sin (* w0 s))
	      w0))))
(defmethod dbetads ((ff g) w0 s)
  "s-derivative of beta given FF, W0, and S"
  (* ff (/ (cos (* w0 s))
	   w0)))
(defmethod deds ((f g) (duds g) (sigma g) (u g))
  "s-derivative of energy"
  (scalar (*i f (*g3 duds sigma (revg u)))))
(defmethod dtmds ((u g))
  "s-derivative of time"
  (norme2 u))

;; Hash table state
(defmethod ksheom ((s number) (x hash-table))
  "KS equations of motion in Hestenes GA form. 
Orbit elements are: alpha (U0), beta (dU0/ds/w0), e (orbit energy), tm (time) 
Expects a variable *data* containing sigma0 (initial orbit frame vector) and forcefun (force function of s, x, and *data*)"
  (with-keys (alpha beta e tm) x
    (let* ((w0 (w0 e))
	   (u (u alpha beta w0 s))
	   (duds (duds alpha beta w0 s))
	   (sigma (spin *ksh-sigma0* alpha))
	   (r (spin sigma u))
	   (rm (norme2 u))
	   (f (funcall *ksh-forcefun* s x))
	   (ff (/ (*g3 f r u) 2)))
      (make-hash
       :alpha (dalphads ff w0 s)
       :beta (dbetads ff w0 s)
       :e (deds f duds sigma u)
       :tm (dtmds u)))))

(defmethod ksheom ((s number) (x kshstate))
  "KS equations of motion in Hestenes GA form. 
Orbit elements are: alpha (U0), beta (dU0/ds/w0), e (orbit energy), tm (time) 
Expects a variable *data* containing sigma0 (initial orbit frame vector) and forcefun (force function of s, x, and *data*)"
  (with-slots (alpha beta e tm) x
    (let* ((w0 (w0 e))
	   (u (u alpha beta w0 s))
	   (duds (duds beta w0))
	   (sigma (spin *ksh-sigma0* alpha))
	   (r (spin sigma u))
	   (f (funcall *ksh-forcefun* s x))
	   (ff (/ (*g3 f r u) 2)))
      (make-instance 
       'kshstate
       :alpha (dalphads ff w0 s)
       :beta (dbetads ff w0 s)
       :e (deds f duds sigma u)
       :tm (dtmds u)))))

;; Spinor/position/velocity conversion functions
(defmethod rv2u ((rv g) (vv g) (basis list))
  "Spinor (U) from position, velocity, and basis vectors"
  (recoverspinor3d (norme rv) (rvbasis rv vv) basis))
(defmethod rv2dudt ((vv g) (u g) (basis1 g))
  "Time derivative of spinor (dU/dt) from velocity, spinor (U), and 1st basis vector"
  (* (*g3 vv u basis1)
     (/ (* 2 (norme2 u)))))
(defun rv2spinors (rv vv basis)
  "Convert position and velocity vectors to spinor (U) and spinor time derivative (dU/dt)"
  (let ((u (rv2u rv vv basis)))
    (make-hash
     :u u
     :dudt (rv2dudt vv u (first basis)))))
(defun duds2dt (duds u)
  "Spinor time derivative (dU/dt) from spinor s-derivative (dU/ds) and spinor (U)"
  (/ duds (norme2 u)))
(defun dudt2ds (dudt u)
  "Spinor s-derivative (dU/ds) from spinot time derivative (dU/dt) and spinor (U)"
  (* dudt (norme2 u)))
(defun spinor2r (u basis1)
  "Position vector from spinor (U) and 1st basis vector"
  (spin basis1 u))
(defmethod spinors2v ((dudt g) (u g) (basisx g))
  "Velocity vector from spinor time derivative (dU/dt), spinor (U), and 1st basis vector"
  (graden (* (*g3 dudt basisx (revg u)) 2d0) 1))
(defun spinors2energy (u duds mu)
  "Specific orbital energy (E) from spinor (U), spinor s-derivative (dU/ds), and gravitational parameter (mu)"
  (/ (- (* 2 (norme2 duds)) mu) (norme2 u)))

;; KSH state conversion functions
(defun rv2ksh (rv vv basis mu)
  "KSH state initialized at time=0 and s=0 from position and velocity vector, list of basis vectors, and gravitational parameter (mu)"
  (with-keys (u dudt) (rv2spinors rv vv basis)
    (let* ((duds (dudt2ds dudt u))
	   (e (spinors2energy u duds mu)))
      (make-hash
       :alpha u 
       :beta (/ duds (sqrt (- (/ e 2))))
       :e e
       :tm 0))))
(defun ksh2spinors (x s)
  "Spinor (U) and spinor s-derivative (dU/ds) from KSH state and S"
  (with-keys (alpha beta e tm) x
    (let ((w0 (w0 e)))
      (make-hash
       :u (u alpha beta w0 s) ; spinor
       :duds (duds alpha beta w0 s))))) ; spinor s-derivative
(defun ksh2rv (x s sigma0)
  "Position and velocity vectors from KSH state, s, and initial orbit position unit vector (sigma0)"
  (with-keys (alpha beta e tm) x
    (with-keys (u duds) (ksh2spinors x s)
      (let ((sigma (spin sigma0 alpha)))
	(make-hash
	 :r (spinor2r u sigma) ; position
	 :v (spinors2v (duds2dt duds u) u sigma)))))) ; velocity

;; Classical orbit elements
(defun coe2rv (sma ecc inc raan aop truan mu basis)
  "Convert cassical orbital elements to position & velocity vectors.
SMA semi major axis
ECC eccentricity
INC inclination
RAAN right ascension of ascending node
AOP argument of perigee
TRUAN true anomaly
MU gravitational parameter
BASIS list of 3 orthogonal basis vectors to express position & velocity in"
  (let* ((r-raan (rotor (*o (first basis) (second basis)) raan))
	 (basis-raan (mapcar #'(lambda (x) (rot x r-raan)) basis))
	 (r-inc (rotor (*o (second basis-raan) (third basis-raan)) inc))
	 (basis-inc (mapcar #'(lambda (x) (rot x r-inc)) basis-raan))
	 (r-aop (rotor (*o (first basis-inc) (second basis-inc)) aop))
	 (basis-aop (mapcar #'(lambda (x) (rot x r-aop)) basis-inc))
	 (r-truan (rotor (*o (first basis-aop) (second basis-aop)) truan))
	 (basis-truan (mapcar #'(lambda (x) (rot x r-truan)) basis-aop))
	 (p (* sma (- 1 (* ecc ecc))))
	 (r (/ p (+ 1 (* ecc (cos truan)))))
	 (ruv (first basis-truan))
	 (rv (* ruv r))
	 (angm (sqrt (/ p mu)))
	 (angmbv (* (*o (first basis-truan) (second basis-truan)) angm))
	 (eccuv (first basis-aop))
	 (eccv (* eccuv ecc))
	 (vv (* (*g (invv angmbv) (+ eccv ruv))
		  mu)))
    (values (graden rv 1) (graden vv 1))))

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
  (acos (/ (scalar (*i (third basis) (dual mombv)))
	   (norme mombv))))
(defun raan-rv (nodev basis)
  "right ascension of ascending node from node vector and basis"
  (let ((tmp (acos (scalar (*i (first basis) (unitg nodev))))))
    (if (< 0 (scalar (*i (second basis) nodev)))
	(- (* 2 pi) tmp)
	tmp)))
(defun aop-rv (nodev eccv basis)
  "argument of perigee from node vector, eccentricity vector, and basis"
  (let ((tmp (acos (scalar (*i (unitg nodev) (unitg eccv))))))
    (if (< 0 (scalar (*i (third basis) eccv)))
	(- (* 2 pi) tmp)
	tmp)))
(defun truan-rv (eccv rv vv)
  "true anomaly from eccentricity, position, and velocity vectors"
  (let ((tmp (acos (scalar (*i (unitg eccv) (unitg rv))))))
    (if (< 0 (scalar (*i rv vv)))
	(- (* 2 pi) tmp)
	tmp)))
(defun rv2coe (rv vv mu basis)
  "Return classical orbital elements given position & velocity vectors, gravitational parameter, and a list of 3 basis vectors"
  (let* ((mombv (mombv-rv rv vv))
	 (nodev (nodev-rv mombv basis))
	 (eccv (eccv-rv rv vv mu)))
    (make-hash
     :sma (sma-rv (norme rv) (norme vv) mu)
     :ecc (norme eccv)
     :inc (inc-rv mombv basis)
     :raan (raan-rv nodev basis)
     :aop (aop-rv nodev eccv basis)
     :truan (truan-rv eccv rv vv))))

;; Test KSH equations of motion
(defparameter *kshtest*
  (make-hash*
   basis (list (ve3 :c1 1) (ve3 :c10 1) (ve3 :c100 1))
   sigma0 (first basis)
   forcefun #'(lambda (s x) (ve3))
   r0 (ve3 :c1 1)
   v0 (ve3 :c10 1.1)
   mu 1
   x0 (rv2ksh r0 v0 basis mu)
   s0 0
   sf (* pi 2))
  "Test data for KSH equations of motion")

(defun testksh (data)
  "Propagate an orbit from test data"
  (with-keys (s0 sf x0 forcefun sigma0 mu) data
    (let ((*ksh-sigma0* sigma0)
	  (*ksh-forcefun* forcefun)
	  (*mu* mu))
      (rka #'ksheom s0 sf x0))))

;; NOTE: Something is very off with the sail forcing
(defparameter *kshsailtest*
  (make-hash*
   basis (list (ve3 :c1 1) (ve3 :c10 1) (ve3 :c100 1))
   sigma0 (first basis)
   alpha 0 ; (* -35.5 (/ pi 180))
   delta 0
   mu 1
   beta 0.1
   forcefun #'(lambda (s x) 
		(with-keys (r v) (ksh2rv x s sigma0)
		  (destructuring-bind (posuv tanuv orbuv) (rvbasis r v)
		    (let ((normuv (+ (* posuv (cos alpha))
				     (* orbuv (* (sin alpha) (cos delta)))
				     (* tanuv (* (sin alpha) (sin delta))))))
		      (* normuv
			 (/ (* beta mu (expt (scalar (*i posuv normuv)) 2))
			    (norme2 r)))))))
;;   forcefun #'sail-ideal-forcefun
   r0 (ve3 :c1 1)
   v0 (ve3 :c10 1)
   x0 (rv2ksh r0 v0 basis mu)
   s0 0
   sf (* pi 2))
  "Test data for KSH equations of motion with solar sail")

(defparameter *kshsail2dtest*
  (make-hash*
   basis (list (ve2 :c1 1) (ve2 :c10 1))
   sigma0 (first basis)
   alpha 0
   beta 0.1
   mu 1
   forcefun_0 #'(lambda (s x) (ve2))
   forcefun_radial #'(lambda (s x)
		       (with-keys (r v) (ksh2rv x s sigma0)
			 (* (unitg r) beta mu (/ (norme2 r)))))
   forcefun #'(lambda (s x)
		       (with-keys (r v) (ksh2rv x s sigma0)
			 (let ((n (rot (unitg r) (rotor (bve2 :c11 1) alpha))))
			   (* n beta mu (expt (cos alpha) 2) (/ (norme2 r))))))
   r0 (ve2 :c1 1)
   v0 (ve2 :c10 1)
   x0 (rv2ksh r0 v0 basis mu)
   s0 0
   sf (* pi 2))
  "Test data for KSH EOM with solar sail in 2D")
		    

(defparameter *ksh-forcefun-sail-normal*
  #'(lambda (s x)
      (with-keys (r v) (ksh2rv x s *ksh-sigma0*)
	(sail-normal-force r *mu* *lightness*))))

(defvar *ephemeris*
  (make-hash
   :alpha (re2 :c0 1)
   :beta (re2 :c11 -1)
   :e -0.5)
  "Ephemeris of a planet")

(defun planet-ksh-eom (s x)
  "Equations of motion of a planet. Only the time propagates as a derivative of s. Used for fast integration of planetary orbits."
  (with-keys (alpha beta e) *ephemeris*
    (dtmds (u alpha beta (w0 e) s))))

(defun ksh-traj-to-cart (result sigma0)
  "Return cartesian coordinate trajectory given KSH trajectory results and orbit reference vector"
  (loop for (s x) in result
     collect (list (gethash :tm x)
		   (ksh2rv x s sigma0))))

