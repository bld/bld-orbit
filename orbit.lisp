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

(defgeneric to-spinor (s x) (:documentation "Convert a state to spinor state"))

(defgeneric to-ks (s x) (:documentation "Convert to a Kustaanheimo-Stiefel orbit element state"))

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

(defun gravity (s r)
  "Gravitational acceleration"
  (with-slots (mu) *sc*
    (- (* (/ mu (norme2 r)) (unitg r)))))

(defun sailidealacc (s r v)
  "Sail acceleration with fixed orientation relative to sail position frame"
  (with-slots (lightness mu pointfun) *sc*
    (let ((n (funcall pointfun s r v)))
      (* lightness mu (/ (norme2 r))
	 (expt (scalar (*i (unitg r) n)) 2)
	 n))))

;; Pointing functions

(defun sailpointingfixed (s r v)
  "Sail pointed at fixed orientation with respect to position/velocity frame defined by rotor RS"
  (with-slots (rs basis) *sc*
    (let* ((rvbasis (rvbasis r v))
	   (rrv (recoverrotor3d rvbasis basis)))
      (rot (first basis) (*g rrv rs)))))

(defun sailpointingnormal (s r v)
  "Sail normal to sunlight"
  (unitg r))

(defun sailpointingtable (s r v)
  "Return sail normal vector"
  (with-slots (rs basis) *sc*
    (let* ((rvbasis (rvbasis r v))
	   (rrv (recoverrotor3d rvbasis basis))
	   (rsi (second (find s rs :test #'<= :key #'first))))
      (rot (first basis) (*g rrv rsi)))))

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
   (rs :initarg :rs :initform (re2 :c0 1d0) :documentation "Sail orientation rotor wrt orbital position frame"))
  (:documentation "Solar sail orbit problem"))

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
       :v (+ (funcall accfun tm r v)
	     (gravity tm r))))))

(defmethod to-cartesian (tm (x cartstate))
  (values x tm))

(defparameter *sail-2d-cart-kepler-eg*
  (make-instance 
   'sail
   :tf (* 8 pi)
   :x0 (make-instance
	'cartstate
	:r (ve2 :c1 1)
	:v (ve2 :c10 1.1))))

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
	     (v (* 2 (*g3 duds (first basis) (revg u))))
	     (f (funcall accfun s r v)))
      (make-instance
       'spinorstate
       :u duds
       :duds (/ (+ (*g3 f r u) (* e u)) 2)
       :tm (norme2 u))))))

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

(defparameter *sail-2d-spin-kepler-eg*
  (make-instance
   'sail
   :tf (* 8 pi)
   :accfun #'(lambda (s r v) (ve2))
   :lightness 0.0d0
   :x0 (let ((*sc* *sail-2d-cart-kepler-eg*))
	 (to-spinor 0 (slot-value *sc* 'x0)))))

(defparameter *sail-2d-spin-normal-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
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
   :tf (* 4 pi)
   :lightness 0.1d0
   :x0 (make-instance 
	'spinorstate
	:u (re2 :c0 1d0)
	:duds (re2 :c11 -0.5d0)
	:tm 0)
   :rs (rotor (bve2 :c11 1d0) (* 35.5 (/ pi 180d0)))))

;; Kustaanheimo-Stiefel Orbit Element equations of motion
;; WORK IN PROGRESS - PERTURBED CASE DIVERGES FROM CARTESIAN & SPINOR FORMULATIONS

(defclass ksstate ()
  ((tm :initarg :tm :documentation "Time")
   (alpha :initarg :alpha :documentation "U_0")
   (beta :initarg :beta :documentation "dU/ds_0 / omega")
   (e :initarg :e :documentation "Work done by the perturbing force"))
  (:documentation "Kustaanheimo-Stiefel Orbit Element State"))

(defstatearithmetic ksstate (tm alpha beta e))

(defmethod print-object ((x ksstate) stream)
  (with-slots (tm alpha beta e) x
    (format stream "#<KSSTATE :TM ~a :ALPHA ~a :BETA ~a :E ~a>" tm alpha beta e)))

(defmethod eom (s (x ksstate))
  "Kustaanheimo-Stiefel Keplerian orbit element equations of motion"
  (with-slots (tm alpha beta e) x
    (with-slots (basis accfun x0) *sc*
      (let* ((w (/ (- (slot-value x0 'e) e) 2d0)) ; work done by disturbing force
	     (omega (sqrt (/ e -2d0))) ; average orbit angular velocity
	     (u (+ (* alpha (cos (* omega s))) (* beta (sin (* omega s))))) ; spinor of orbit position
	     (up (* omega (- (* beta (cos (* omega s))) (* alpha (sin (* omega s)))))) ; du/ds
	     (r (spin (first basis) u)) ; position vector
	     (rm (norme2 u)) ; position vector magnitude (radius)
	     (v (* 2 (/ rm) (*g3 up (first basis) (revg u)))) ; velocity vector
	     (f (funcall accfun s r v)) ; perturbing force
	     (ff (/ (- (*g3 f r u) (* w u)) 2d0))) ; KS perturbing force
	(make-instance
	 'ksstate
	 :tm rm
	 :alpha (* -1 ff (/ omega) (sin (* omega s)))
	 :beta (* ff (/ omega) (cos (* omega s)))
	 :e (* rm (scalar (*i f v))))))))

(defmethod to-cartesian (s (x ksstate))
  "Convert Kustaanheimo-Stiefel orbit element state to cartesian"
  (with-slots (tm alpha beta e) x
    (with-slots (basis) *sc*
      (let* ((omega (sqrt (/ e -2d0)))
	     (u (+ (* alpha (cos (* omega s))) (* beta (sin (* omega s)))))
	     (up (* omega (- (* beta (cos (* omega s))) (* alpha (sin (* omega s)))))))
	(values
	 (make-instance
	  'cartstate
	  :r (spin (first basis) u)
	  :v (* 2d0 (/ (norme2 u)) (*g3 up (first basis) (revg u))))
	 tm)))))

(defmethod to-ks (tm (x cartstate))
  "Make KS element state from time & cartesian state"
  (with-slots (r v) x
    (with-slots (mu basis) *sc*
      (let* ((rvbasis (rvbasis r v))
	     (u (recoverspinor3d (norme r) rvbasis basis))
	     (up (/ (*g3 v u (first basis)) 2d0))
	     (e (/ (- (* 2 (norme2 up)) mu) (norme2 u)))
	     (omega (sqrt (/ e -2d0))))
	(make-instance
	 'ksstate
	 :tm tm
	 :e e
	 :alpha u
	 :beta (/ up omega))))))

(defparameter *sail-2d-ks-kepler-eg*
  (make-instance
   'sail
   :tf (* 8 pi)
   :accfun #'(lambda (s r v) (ve2))
   :lightness 0.0d0
   :x0 (let ((*sc* *sail-2d-cart-kepler-eg*))
	 (to-ks 0 (slot-value *sc* 'x0))))
  "Unperturbed Kepler problem")

(defparameter *sail-2d-ks-normal-eg*
  (make-instance
   'sail
   :tf (* 4 pi)
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
   :tf (* 4 pi)
   :lightness 0.1d0
   :x0 (make-instance
	'ksstate
	:tm 0d0
	:alpha (re2 :c0 1d0)
	:beta (re2 :c11 -1d0)
	:e -0.5d0)
   :rs (rotor (bve2 :c11 1d0) (* 35.5 (/ pi 180d0))))
  "Sail pointed at fixed orientation to sunlight")

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

