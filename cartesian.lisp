(in-package :bld-orbit)

(defclass cartesian-state ()
  ((r :initarg :r :accessor r)
   (v :initarg :v :accessor v)))

(defmethod print-object ((x cartesian-state) stream)
  (format stream "#<CARTESIAN-STATE :R ~a :V ~a>" (r x) (v x)))

(defmethod norminfx ((g g)) (norminf g))

(defstatearithmetic cartesian-state (r v))

(defclass cartesian-problem ()
  ((central-body :initarg :central-body :documentation "Central body")
   (ref :initarg :ref :documentation "Reference frame")
   (utc0 :initarg :utc0 :documentation "Initial UTC time")
   (utcf :initarg :utcf :documentation "Final UTC time")
   (eom :initarg :eom :documentation "Equation of motion function")
   (x0 :initarg :x0 :documentation "Initial state at UTC0")
   (hmin :initarg :hmin :documentation "Minimum step size")
   (hmax :initarg :hmax :documentation "Maximum step size")
   (tol :initarg :tol :documentation "Integration tolerance")))

(defgeneric to-cartesian (utc b &key))

(defmethod to-cartesian (utc (b body)
			 &key (ref (slot-value b 'ref))
			   (observer (slot-value b 'center)))
  "Make a cartesian state from UTC time and body"
  (with-slots (name ephemeris) b
    (with-kernel ephemeris
      (let ((rv (spk-ezr name (utc-to-epht utc) observer :ref ref)))
	(make-instance
	 'cartesian-state
	 :r (ve3 :e1 (aref rv 0) :e2 (aref rv 1) :e3 (aref rv 2))
	 :v (ve3 :e1 (aref rv 3) :e2 (aref rv 4) :e3 (aref rv 5)))))))

(defmethod eom (utc (x cartesian-state) (p cartesian-problem))
  "Cartesian equations of motion"
  (with-slots (central-body ref) p
    (with-slots (r v) x
      (make-instance
       'cartesian-state
       :r v
       :v (gravity r central-body)))))

(defmethod propagate ((p cartesian-problem))
  "Propagate a cartesian problem"
  (with-slots (eom utc0 utcf x0 hmin hmax tol central-body) p
    ;; Pre-load central body kernel to speed up calculations
    (with-kernel (slot-value central-body 'ephemeris)
      (rka eom utc0 utcf x0 :param p :hmin hmin :hmax hmax :tol tol))))

(defmethod convert-results (results (problem cartesian-problem) &key (body (slot-value problem 'central-body)) (time-fn #'utc-to-epht) (observer (slot-value body 'center)) (ref :eclipj2000))
  "Convert results to alternative format given results and problem used to generate them. Keys specify:
BODY: Body to reference coordinates to (default: central body)
TIME-FN: Function to convert time from UTC (default: UTC-TO-EPHT)
OBSERVER: Frame to reference body to (default: body's center)
REF: reference frame (default: :eclipj2000)"
  (loop for (utc x) in results
     for x-b = (to-cartesian utc body :observer observer :ref ref)
     collect (list (funcall time-fn utc) (+ x x-b))))
