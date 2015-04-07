;;;; Kustaanheimo-Stiefel equations of motion

(defclass ks-state ()
  ((utc :initarg :utc)))

(defstatearithmetic ks-state (utc))

(defmethod eom (s (x ks-state) p)
  "Equations of motion for Kepler orbit using Kustaanheimo-Stiefel equations"
  (with-slots (alpha beta w0) p
    (let ((u (+ (* alpha (cos (* w0 s)))
		(* beta (sin (* w0 s))))))
      (make-instance 'ks-state :utc (norme2 u)))))

(defclass ks-problem (spinor-problem)
  ((alpha :initarg :alpha :documentation "U at time 0")
   (beta :initarg :beta :documentation "dU/dS at time 0")
   (w0 :initarg :w0 :documentation "specific angular velocity")
   (x0ks :initarg :x0ks :documentation "initial KS state")))

(defmethod to-initial-ks-state (s (x0s spinor-state) (p spinor-problem))
  (with-slots (u duds utc) x0s
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
	   (make-instance 'ks-state :utc utc)
	   0 ; s = 0
	   alpha beta w0)))))) ; return calculated orbit constants

(defmethod to-spinor-state (s (x ks-state) &key problem)
  (with-slots (utc) x
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
	  :utc utc)
	 s)))))

(defmethod to-ks-problem ((p spinor-problem))
  (with-slots
	(central-body ; Cartesian
	 ref utc0 utcf eom x0 hmin hmax tol
	 ;; Spinor
	 x0s s0 sfmax stopfn stopval stoptest stoptol) p
    (multiple-value-bind (x0ks s0ks alpha beta w0)
	(to-initial-ks-state s0 x0s p)
      (make-instance
       'ks-problem
       ;; Cartesian problem slots
       :central-body central-body
       :ref ref
       :utc0 utc0
       :utcf utcf
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
  (with-slots (eom s0 sfmax x0ks stopfn stopval stoptest stoptol tol central-body) p
    (with-kernel (slot-value central-body 'ephemeris)
      (rka-stop-nr
       eom s0 sfmax x0ks
       :param p
       :stopfn stopfn
       :stopval stopval
       :stoptest stoptest
       :stoptol stoptol
       :tol tol))))

(defmethod to-spinor-results (results (p ks-problem))
  (loop for (s xks) in results
     collect (reverse (multiple-value-list (to-spinor-state s xks :problem p)))))
