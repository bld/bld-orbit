;;; Kustaanheimo-Stiefel equations of motion for open and closed orbits using Stumpff functions

(in-package :bld-orbit)

;; Problem classes

(defclass ks-kepler-stumpff-problem (spinor-problem)
  ((alpha :initarg :alpha :documentation "U at time 0")
   (beta :initarg :beta :documentation "dU/dS at time 0")
   (a0 :initarg :a0 :documentation "initial negative half of specific orbital energy")
   (x0ks :initarg :x0ks :documentation "initial KS state"))
  (:documentation "Kustaanheimo-Stiefel Kepler problem using Stumpff functions for open or closed orbits"))

(defclass ks-stumpff-problem (spinor-problem)
  ((a0 :initarg :a0 :documentation "initial negative half of specific orbital energy")
   (x0ks :initarg :x0ks :documentation "initial KS state")
   (accelfn :initarg :accelfn :documentation "perturbation acceleration function"))
  (:documentation "Kustaanheimo-Stiefel perturbed orbit element problem for open orbits using Stumpff functions"))

;; States

(defclass ks-kepler-stumpff-state ()
  ((et :initarg :et)))

(defmethod print-object ((x ks-kepler-stumpff-state) stream)
  (with-slots (et) x
    (format stream "#<KS-KEPLER-STATE :ET ~a>" et)))

(defstatearithmetic ks-kepler-stumpff-state (et))

(defclass ks-stumpff-state ()
  ((et :initarg :et :documentation "Ephemeris time")
   (alpha :initarg :alpha :documentation "Negative of initial spinor, -U at s=0")
   (beta :initarg :beta :documentation "Initial spinor derivative, dU/ds at s=0")
   (e :initarg :e :documentation "Specific Kepler orbital energy")))

(defmethod print-object ((x ks-stumpff-state) stream)
  (with-slots (et alpha beta e) x
    (format stream "#<KS-STUMPFF-STATE :ET ~a :ALPHA ~a :BETA ~a :E ~a>" et alpha beta e)))

(defstatearithmetic ks-stumpff-state (et alpha beta e))

;; Equations of motion

(defmethod eom-kepler (s (x ks-kepler-stumpff-state) (p ks-kepler-stumpff-problem))
  "Equations of motion for Kepler orbit using Kustaanheimo-Stiefel equations and Stumpff functions for parabolic & hyperbolic orbits"
  (with-slots (u) (to-spinor-state s x :problem p)
    (make-instance 'ks-kepler-stumpff-state :et (norme2 u))))

(defmethod eom-kepler (s (x ks-stumpff-state) (p ks-stumpff-problem))
  (let* ((xs (to-spinor-state s x :problem p))
	 (xc (to-cartesian-state s xs)))
    (make-instance
     'ks-stumpff-state
     :et (norme (slot-value xc 'r))
     :alpha (re3)
     :beta (re3)
     :e 0)))

(defmethod eom-nbody (s (x ks-stumpff-state) (p ks-stumpff-problem))
  (with-slots (a0 accelfn x0ks) p
    (with-slots (et alpha beta e) x
      (let* ((xs (to-spinor-state s x :problem p)) ; spinor state
	     (xc (to-cartesian-state s xs)) ; cartesian state
	     (r (slot-value xc 'r))
	     (v (slot-value xc 'v))
	     (rm (norme r))
	     (hk (- e))
	     (w (- (/ hk 2) a0))
	     (f (nbody-accel et r p)) ; n-body acceleration
	     (u (slot-value xs 'u))
	     (ff (- (/ (*g f r u) 2) (* w u)))
	     (s2 (expt s 2))
	     (c_1 (c1 (* a0 s2)))
	     (c_0 (c0 (* a0 s2))))
	(make-instance ; KS state derivative
	 'ks-stumpff-state
	 :et rm ; ephemeris time
	 :alpha (* ff c_1 s) ; derivative of negative initial spinor
	 :beta (* ff c_0) ; derivative of initial spinor derivative
	 :e (* rm (scalar (*i v f)))))))) ; specific orbit energy

;; Convert between spinor states to KS states in Stumpff problem

(defmethod to-initial-ks-kepler-stumpff-state (s (x0s spinor-state) (p spinor-problem))
  "Convert to initial KS Kepler Stumpff problem & state from s, initial spinor state, and spinor problem."
  (with-slots (u duds et) x0s
    (with-slots (central-body) p
      (with-slots (mu) central-body
	(let* ((rm (norme2 u)) ; Radius
	       (e (/ (- (* 2 (norme2 duds)) mu) rm)) ; Specific Kepler orbit energy
	       (a0 (- (/ e 2))) ; -E/2, called alpha_0 in Bond 1973
	       (s2 (expt s 2)) ; s^2
	       (c_0 (c0 (* a0 s2))) ; c0 Stumpff function
	       (c_1 (c1 (* a0 s2))) ; c1 Stumpff function
	       (alpha (- (* duds s c_1) (* u c_0))) ; -U at s=0
	       (beta (+ (* duds c_0) (* u s c_1 a0)))) ; dU/ds at s=0
	  (values
	   (make-instance 'ks-kepler-stumpff-state :et et)
	   0 ; set new KS state to s = 0
	   alpha beta a0)))))) ; return calculated orbit constants

(defmethod to-initial-ks-stumpff-state (s (x0s spinor-state) (p spinor-problem))
  (with-slots (u duds et) x0s
    (with-slots (central-body) p
      (with-slots (mu) central-body
	(let* ((rm (norme2 u)) ; Radius
	       (e (/ (- (* 2 (norme2 duds)) mu) rm)) ; Specific Kepler orbit energy
	       (a0 (- (/ e 2))) ; -E/2, called alpha_0 in Bond 1973
	       (s2 (expt s 2)) ; s^2
	       (c_0 (c0 (* a0 s2))) ; c0 Stumpff function
	       (c_1 (c1 (* a0 s2))) ; c1 Stumpff function
	       (alpha (- (* duds s c_1) (* u c_0))) ; -U at s=0
	       (beta (+ (* duds c_0) (* u s c_1 a0)))) ; dU/ds at s=0
	  (values
	   (make-instance 'ks-stumpff-state :et et :alpha alpha :beta beta :e e)
	   0 ; s=0
	   a0)))))) ; return initial negative half Kepler energy

(defmethod to-spinor-state (s (x ks-kepler-stumpff-state) &key problem)
  "Return spinor state and s values given s, KS Kepler Stumpff state, and KS Kepler Stumpff problem"
  (with-slots (et) x
    (with-slots (alpha beta a0) problem
      (let* ((s2 (expt s 2))
	     (c_1 (c1 (* a0 s2)))
	     (c_0 (c0 (* a0 s2)))
	     (u (- (* beta s c_1) (* alpha c_0)))
	     (duds (+ (* beta c_0) (* a0 alpha s c_1))))
	(values (make-instance 'spinor-state :u u :duds duds :et et) s)))))

(defmethod to-spinor-state (s (x ks-stumpff-state) &key problem)
  (with-slots (et alpha beta e) x
    (with-slots (a0) problem
      (let* ((s2 (expt s 2))
	     (c_1 (c1 (* a0 s2)))
	     (c_0 (c0 (* a0 s2)))
	     (u (- (* beta s c_1) (* alpha c_0)))
	     (duds (+ (* beta c_0) (* a0 alpha s c_1))))
	(values (make-instance 'spinor-state :u u :duds duds :et et) s)))))

;; Convert between spinor and KS Stumpff problems

(defmethod to-ks-kepler-stumpff-problem ((p spinor-problem))
  (with-slots (central-body nbodies spk lsk ref et0 etf eom x0 hmin hmax tol ; Cartesian
			    x0s s0 sfmax stopfn stopval stoptest stoptol) p ; Spinor
    ;; Compute initial KS state & required parameters
    (multiple-value-bind (x0ks s0ks alpha beta a0)
	(to-initial-ks-kepler-stumpff-state s0 x0s p)
      ;; Create KS Kepler Stumpff problem
      (make-instance
       'ks-kepler-stumpff-problem
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
       ;; KS Kepler Stumpff problem slots
       :alpha alpha
       :beta beta
       :a0 a0
       :x0ks x0ks))))

(defmethod to-ks-stumpff-problem ((p spinor-problem))
  (with-slots
	(;; Cartesian
	 central-body nbodies spk lsk ref et0 etf eom x0 hmin hmax tol
	 ;; Spinor
	 x0s s0 sfmax stopfn stopval stoptest stoptol) p
    ;; Compute initial KS state & required parameters
    (multiple-value-bind (x0ks s0ks a0)
	(to-initial-ks-stumpff-state s0 x0s p)
      ;; Create KS Kepler Stumpff problem
      (make-instance
       'ks-stumpff-problem
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
       ;; KS Kepler Stumpff problem slots
       :a0 a0
       :x0ks x0ks))))

(defmethod to-ks-kepler-stumpff-problem ((p cartesian-problem))
  "Convert Cartesian problem to KS Kepler Stumpff problem"
  (to-ks-kepler-stumpff-problem (to-spinor-problem p)))

(defmethod to-ks-stumpff-problem ((p cartesian-problem))
  (to-ks-stumpff-problem (to-spinor-problem p)))

;; Propagate KS Stumpff problem

(defmethod propagate ((p ks-kepler-stumpff-problem))
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
	  :tol tol
	  :a bld-ode::*a-dp78*
	  :bl bld-ode::*bl-dp78*
	  :bh bld-ode::*bh-dp78*
	  :c bld-ode::*c-dp78*)
      ;; Unload kernels
      (progn
	(mapcar #'unload spk)
	(unload lsk)))))

(defmethod propagate ((p ks-stumpff-problem))
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

;; Convert output of KS problem propagation to spinor results

(defmethod to-spinor-results (results (p ks-kepler-stumpff-problem))
  (loop for (s xks) in results
     collect (reverse (multiple-value-list (to-spinor-state s xks :problem p)))))

(defmethod to-spinor-results (results (p ks-stumpff-problem))
  (loop for (s xks) in results
     collect (reverse (multiple-value-list (to-spinor-state s xks :problem p)))))
