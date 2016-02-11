;;;; Kustaanheimo-Stiefel equations of motion

(in-package :bld-orbit)

(defclass ks-kepler-state ()
  ((et :initarg :et)))

(defmethod print-object ((x ks-kepler-state) stream)
  (with-slots (et) x
    (format stream "#<KS-KEPLER-STATE :ET ~a>" et)))

(defstatearithmetic ks-kepler-state (et))

(defmethod eom-kepler (s (x ks-kepler-state) p)
  "Equations of motion for Kepler orbit using Kustaanheimo-Stiefel equations"
  (with-slots (alpha beta w0) p
    (let ((u (+ (* alpha (cos (* w0 s)))
		(* beta (sin (* w0 s))))))
      (make-instance 'ks-kepler-state :et (norme2 u)))))

(defclass ks-problem (spinor-problem)
  ((alpha :initarg :alpha :documentation "U at time 0")
   (beta :initarg :beta :documentation "dU/dS at time 0")
   (w0 :initarg :w0 :documentation "specific angular velocity")
   (x0ks :initarg :x0ks :documentation "initial KS state")))

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

(defmethod to-ks-problem ((p spinor-problem))
  (with-slots
	(;; Cartesian
	 central-body nbodies spk lsk ref et0 etf eom x0 hmin hmax tol
	 ;; Spinor
	 x0s s0 sfmax stopfn stopval stoptest stoptol) p
    (multiple-value-bind (x0ks s0ks alpha beta w0)
	(to-initial-ks-kepler-state s0 x0s p)
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
       :alpha alpha
       :beta beta
       :w0 w0
       :x0ks x0ks))))

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

(defmethod to-spinor-results (results (p ks-problem))
  (loop for (s xks) in results
     collect (reverse (multiple-value-list (to-spinor-state s xks :problem p)))))
