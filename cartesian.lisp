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
   (nbodies :initarg :nbodies :documentation "Additional bodies to treat as perturbing point masses")
   (spk :initarg :spk :documentation "SPK kernels to load for problem")
   (lsk :initarg :lsk :documentation "LSK kernels to load for problem")
   (ref :initarg :ref :documentation "Reference frame")
   (et0 :initarg :et0 :documentation "Initial ephemeris time (seconds past J2000)")
   (etf :initarg :etf :documentation "Final ephemeris time (seconds past J2000)")
   (eom :initarg :eom :documentation "Equation of motion function")
   (x0 :initarg :x0 :documentation "Initial state at UTC0")
   (hmin :initarg :hmin :documentation "Minimum step size")
   (hmax :initarg :hmax :documentation "Maximum step size")
   (tol :initarg :tol :documentation "Integration tolerance")))

(defgeneric to-cartesian (et b &key))

(defmethod to-cartesian (et (b body)
			 &key (ref (slot-value b 'ref))
			   (observer (slot-value b 'center)))
  "Make a cartesian state from UTC time and body"
  (with-slots (name) b
    (let ((rv (spk-ezr name et observer :ref ref)))
      (make-instance
       'cartesian-state
       :r (ve3 :e1 (aref rv 0) :e2 (aref rv 1) :e3 (aref rv 2))
       :v (ve3 :e1 (aref rv 3) :e2 (aref rv 4) :e3 (aref rv 5))))))

(defmethod eom-kepler (et (x cartesian-state) (p cartesian-problem))
  "Cartesian Keplerian equations of motion given ephemeris time, cartesian state, and cartesian problem object"
  (with-slots (central-body) p
    (with-slots (r v) x
      (make-instance
       'cartesian-state
       :r v
       :v (gravity r central-body)))))

(defmethod nbody-accel (et (r ve3) p)
  (with-slots (central-body nbodies ref) p
    (loop with obs = (slot-value central-body 'name) ; central body name
       with f = (ve3) ; force vector to increment by each body
       ;; iterate over n-bodies
       for nb in nbodies
       ;; Position of n-body relative to central body
       for r-nb = (position-vector et nb :ref ref :observer obs)
       ;; Position of spacecraft relative to n-body
       for r-sc-nb = (- r r-nb)
       ;; Increment total force
       do (setf f (+ f
		     (gravity r-sc-nb nb) ; Direct n-body force
		     (gravity r-nb nb))) ; Indirect n-body force
       ;; Return total force
       finally (return f))))

(defmethod nbody-battin (et (r ve3) p)
  (with-slots (central-body nbodies ref) p
    (with-slots ((obs name)) central-body
      (loop with acc = (ve3)
	 for nb in nbodies
	 for rho = (position-vector et nb :ref ref :observer obs)
	 for d = (- r rho)
	 for q = (scalar (*i (/ r (norme2 rho)) (- r (* 2 rho))))
	 for f = (* q (/ (+ 3 (* 3 q) (expt q 2))
			 (+ 1 (expt (+ 1 q) 3/2))))
	 do (setf acc
		  (+ acc
		     (with-slots (mu) nb
		       (- (* (/ mu (expt (norme d) 3))
			     (+ r (* f rho)))))))
	 finally (return acc)))))

(defmethod eom-nbody (et (x cartesian-state) (p cartesian-problem))
  "Cartesian equations of motion with central body plus other bodies as perturbing n-body forces"
  (with-slots (central-body nbodies ref) p
    (with-slots (r v) x
      (make-instance
       'cartesian-state
       :r v
       :v (+ (gravity r central-body)
	     (nbody-accel et r p))))))
	     ;;(nbody-battin et r p))))))
       
(defmethod propagate ((p cartesian-problem))
  "Propagate a cartesian problem"
  (with-slots (eom et0 etf x0 hmin hmax tol central-body spk lsk) p
    ;; Pre-load SPK and LSK kernels
    (dolist (k spk) (unless (kernelp k) (furnsh k)))
    (unless (kernelp lsk) (furnsh lsk))
    ;; Integrate
    (unwind-protect
	 (rka eom et0 etf x0 :param p :hmin hmin :hmax hmax :tol tol)
      ;; Unload kernels
      (progn
	(mapcar #'unload spk)
	(unload lsk)))))

(defmethod convert-results (results (problem cartesian-problem) &key (body (slot-value problem 'central-body)) (time-fn #'utc-to-epht) (observer (slot-value body 'center)) (ref :eclipj2000))
  "Convert results to alternative format given results and problem used to generate them. Keys specify:
BODY: Body to reference coordinates to (default: central body)
TIME-FN: Function to convert time from UTC (default: UTC-TO-EPHT)
OBSERVER: Frame to reference body to (default: body's center)
REF: reference frame (default: :eclipj2000)"
  (loop for (utc x) in results
     for x-b = (to-cartesian utc body :observer observer :ref ref)
     collect (list (funcall time-fn utc) (+ x x-b))))

(defmethod to-csv (results (problem cartesian-problem) &key (separator #\,) stream)
  (write-csv
   (append '(("Time" "X" "Y" "Z" "VX" "VY" "VZ"))
	   (loop for (tm x) in results
	      collect
		(with-slots (r v) x
		  (cons tm
			(append
			 (map 'list #'identity (coef r))
			 (map 'list #'identity (coef v)))))))
   :stream stream
   :separator separator))
