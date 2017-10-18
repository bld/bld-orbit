;;;; Kustaanheimo-Stiefel equations of motion

(in-package :bld-orbit)

;; KS problem classes

(defclass ks-kepler-problem (spinor-problem)
  ((alpha :initarg :alpha :documentation "U at time 0")
   (beta :initarg :beta :documentation "dU/dS at time 0")
   (w0 :initarg :w0 :documentation "specific angular velocity")
   (x0ks :initarg :x0ks :documentation "initial KS state"))
  (:documentation "Kustaanheimo-Stiefel Kepler problem for closed orbits"))

(defclass ks-problem (spinor-problem)
  ((w0 :initarg :w0 :documentation "initial specific angular velocity")
   (x0ks :initarg :x0ks :documentation "initial KS state")
   (accelfn :initarg :accelfn :documentation "perturbation acceleration function"))
  (:documentation "Kustaanheimo-Stiefel perturbed orbit element problem for closed orbits"))

;; KS Kepler state

(defclass ks-kepler-state ()
  ((et :initarg :et)))

(defmethod print-object ((x ks-kepler-state) stream)
  (with-slots (et) x
    (format stream "#<KS-KEPLER-STATE :ET ~a>" et)))

(defstatearithmetic ks-kepler-state (et))

(defmethod eom-kepler (s (x ks-kepler-state) (p ks-kepler-problem))
  "Equations of motion for Kepler orbit using Kustaanheimo-Stiefel equations"
  (with-slots (alpha beta w0) p
    (let ((u (+ (* alpha (cos (* w0 s)))
		(* beta (sin (* w0 s))))))
      (make-instance 'ks-kepler-state :et (norme2 u)))))

;; KS general perturbed state

(defclass ks-state ()
  ((et :initarg :et :documentation "Ephemeris time")
   (alpha :initarg :alpha :documentation "Initial spinor")
   (beta :initarg :beta :documentation "Initial spinor derivative")
   (e :initarg :e :documentation "Specific orbital energy")))

(defmethod print-object ((x ks-state) stream)
  (with-slots (et alpha beta e) x
    (format stream "#<KS-STATE :ET ~a :ALPHA ~a :BETA ~a :E ~a>" et alpha beta e)))

(defstatearithmetic ks-state (et alpha beta e))

(defmethod eom-kepler (s (x ks-state) (p ks-kepler-problem))
  (let* ((xs (to-spinor-state s x :problem p))
	 (xc (to-cartesian-state s xs)))
    (make-instance
     'ks-state
     :et (norme (slot-value xc 'r))
     :alpha (re3)
     :beta (re3)
     :e 0)))

(defmethod eom-nbody (s (x ks-state) (p ks-problem))
  (with-slots (w0 accelfn x0ks) p
    (with-slots (et alpha beta e) x
      (let* ((xs (to-spinor-state s x :problem p)) ; spinor state
	     (xc (to-cartesian-state s xs)) ; cartesian state
	     (r (slot-value xc 'r))
	     (v (slot-value xc 'v))
	     (rm (norme r))
	     (hk (- e))
	     (w (- (/ hk 2) (expt w0 2)))
	     (f (nbody-accel et r p)) ; n-body acceleration
	     (u (slot-value xs 'u))
	     (ff (- (/ (*g f r u) 2) (* w u))))
	(make-instance ; KS state derivative
	 'ks-state
	 :et rm ; ephemeris time
	 :alpha (- (* (/ ff w0) (sin (* w0 s)))) ; initial spinor
	 :beta (* (/ ff w0) (cos (* w0 s))) ; initial spinor derivative
	 :e (* rm (scalar (*i v f)))))))) ; specific orbit energy

;; Convert between KS and other states

(defmethod to-initial-ks-kepler-state (s (x0s spinor-state) (p spinor-problem))
  (with-slots (u duds et) x0s
    (with-slots (central-body) p
      (with-slots (mu) central-body
	(let* ((rm (norme2 u))
	       (e (/ (- (* 2 (norme2 duds)) mu) rm))
	       (hk0 (- e))
	       (w0 (sqrt (/ hk0 2)))
	       (alpha (- (* u (cos (* w0 s)))
			 (* (/ duds w0) (sin (* w0 s)))))
	       (beta (+ (* u (sin (* w0 s)))
			(* (/ duds w0) (cos (* w0 s))))))
	  (values
	   (make-instance 'ks-kepler-state :et et)
	   0 ; s = 0
	   alpha beta w0)))))) ; return calculated orbit constants

(defmethod to-initial-ks-state (s (x0s spinor-state) (p spinor-problem))
  (with-slots (u duds et) x0s
    (with-slots (central-body) p
      (with-slots (mu) central-body
	(let* ((rm (norme2 u))
	       (e (/ (- (* 2 (norme2 duds)) mu) rm))
	       (hk0 (- e))
	       (w0 (sqrt (/ hk0 2)))
	       (alpha (- (* u (cos (* w0 s)))
			 (* (/ duds w0) (sin (* w0 s)))))
	       (beta (+ (* u (sin (* w0 s)))
			(* (/ duds w0) (cos (* w0 s))))))
	  (values
	   (make-instance 'ks-state :et et :alpha alpha :beta beta :e e)
	   0 ; s = 0
	   w0)))))) ; return calculated natural frequency

(defmethod to-spinor-state (s (x ks-kepler-state) &key problem)
  (with-slots (et) x
    (with-slots (alpha beta w0) problem
      (let ((u (+ (* alpha (cos (* w0 s)))
		  (* beta (sin (* w0 s)))))
	    (duds (* w0 (- (* beta (cos (* w0 s)))
			   (* alpha (sin (* w0 s)))))))
	(values
	 (make-instance
	  'spinor-state
	  :u u
	  :duds duds
	  :et et)
	 s)))))

(defmethod to-spinor-state (s (x ks-state) &key problem)
  "Convert KS general state to spinor state"
  (with-slots (et alpha beta e) x
    (with-slots (w0) problem
      (values
       (make-instance
	'spinor-state
	:u (+ (* alpha (cos (* w0 s)))
	      (* beta (sin (* w0 s))))
	:duds (* w0 (- (* beta (cos (* w0 s)))
		       (* alpha (sin (* w0 s)))))
	:et et)
       s))))

(defmethod to-ks-kepler-problem ((p spinor-problem))
  (with-slots
	(;; Cartesian
	 central-body nbodies spk lsk ref et0 etf eom x0 hmin hmax tol
	 ;; Spinor
	 x0s s0 sfmax stopfn stopval stoptest stoptol) p
    (multiple-value-bind (x0ks s0ks alpha beta w0)
	(to-initial-ks-kepler-state s0 x0s p)
      (make-instance
       'ks-kepler-problem
       ;; Cartesian problem slots
       :central-body central-body
       :nbodies nbodies
       :spk spk
       :lsk lsk
       :ref ref
       :et0 et0
       :etf etf
       :eom eom
       :x0 x0
       :hmin hmin
       :hmax hmax
       :tol tol
       ;; Spinor problem slots
       :x0s x0s
       :s0 s0ks
       :sfmax sfmax
       :stopfn stopfn
       :stopval stopval
       :stoptest stoptest
       :stoptol stoptol
       ;; KS problem slots
       :alpha alpha
       :beta beta
       :w0 w0
       :x0ks x0ks))))

(defmethod to-ks-problem ((p spinor-problem))
  (with-slots
	(;; Cartesian
	 central-body nbodies spk lsk ref et0 etf eom x0 hmin hmax tol
	 ;; Spinor
	 x0s s0 sfmax stopfn stopval stoptest stoptol) p
    (multiple-value-bind (x0ks s0ks w0)
	(to-initial-ks-state s0 x0s p)
      (make-instance
       'ks-problem
       ;; Cartesian problem slots
       :central-body central-body
       :nbodies nbodies
       :spk spk
       :lsk lsk
       :ref ref
       :et0 et0
       :etf etf
       :eom eom
       :x0 x0
       :hmin hmin
       :hmax hmax
       :tol tol
       ;; Spinor problem slots
       :x0s x0s
       :s0 s0ks
       :sfmax sfmax
       :stopfn stopfn
       :stopval stopval
       :stoptest stoptest
       :stoptol stoptol
       ;; KS problem slots
       :w0 w0
       :x0ks x0ks))))

(defmethod to-ks-kepler-problem ((p cartesian-problem))
  "Convert cartesian problem to Kustaanheimo-Stiefel problem"
  (to-ks-kepler-problem (to-spinor-problem p)))

(defmethod propagate ((p ks-kepler-problem))
  (with-slots (eom s0 sfmax x0ks stopfn stopval stoptest stoptol tol central-body spk lsk) p
    ;; Pre-load SPK and LSK kernels
    (dolist (k spk) (unless (kernelp k) (furnsh k)))
    (unless (kernelp lsk) (furnsh lsk))
    ;; Integrate
    (unwind-protect
	 (rka-stop-nr
	  eom s0 sfmax x0ks
	  :param p
	  :stopfn stopfn
	  :stopval stopval
	  :stoptest stoptest
	  :stoptol stoptol
	  :tol tol)
      ;; Unload kernels
      (progn
	(mapcar #'unload spk)
	(unload lsk)))))

(defmethod propagate ((p ks-problem))
  (with-slots (eom s0 sfmax x0ks stopfn stopval stoptest stoptol tol central-body spk lsk) p
    ;; Pre-load SPK and LSK kernels
    (dolist (k spk) (unless (kernelp k) (furnsh k)))
    (unless (kernelp lsk) (furnsh lsk))
    ;; Integrate
    (unwind-protect
	 (rka-stop-nr
	  eom s0 sfmax x0ks
	  :param p
	  :stopfn stopfn
	  :stopval stopval
	  :stoptest stoptest
	  :stoptol stoptol
	  :tol tol)
      ;; Unload kernels
      (progn
	(mapcar #'unload spk)
	(unload lsk)))))

(defmethod to-spinor-results (results (p ks-kepler-problem))
  (loop for (s xks) in results
     collect (reverse (multiple-value-list (to-spinor-state s xks :problem p)))))

(defmethod to-spinor-results (results (p ks-problem))
  (loop for (s xks) in results
     collect (reverse (multiple-value-list (to-spinor-state s xks :problem p)))))
