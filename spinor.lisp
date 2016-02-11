;;;; Spinor equations of motion

(in-package :bld-orbit)

(defclass spinor-state ()
  ((u :initarg :u)
   (duds :initarg :duds)
   (et :initarg :et)))

(defmethod print-object ((x spinor-state) stream)
  (with-slots (u duds et) x
    (format stream "#<SPINOR-STATE :U ~a :DUDS ~a :ET ~a>" u duds et)))

(defstatearithmetic spinor-state (u duds et))

(defparameter *iframe*
  (list (ve3 :e1 1)
	(ve3 :e2 1)
	(ve3 :e3 1)))

(defmethod eom-kepler (s (x spinor-state) p)
  "Keplerian spinor equations of motion given independent variable S, spinor state, and spinor problem"
  (with-slots (u duds et) x
    (with-slots (central-body) p
      (with-slots (mu) central-body
	(let* ((rm (norme2 u))
	       (e (/ (- (* 2 (norme2 duds)) mu) rm)))
	  (make-instance
	   'spinor-state
	   :u duds
	   :duds (/ (* e u) 2)
	   :et rm))))))

(defmethod nbody-accel (s (x spinor-state) p)
  "Calculate n-body accelerations"
  (with-slots (u et) x
    (with-slots (central-body nbodies ref) p
      (with-slots (name) central-body
	(loop with r = (spin (first *iframe*) u)
	   with a = (ve3)
	   for nb in nbodies
	   for r-nb = (position-vector et nb :ref ref :observer name)
	   for r-sc-nb = (- r r-nb)
	   do (setf a (+ a
			 (gravity r-sc-nb nb) ; direct
			 (gravity r-nb nb))) ; indirect
	   finally (return a))))))
	     
(defmethod eom-nbody (s (x spinor-state) p)
  (with-slots (u duds et) x
    (with-slots (central-body) p
      (with-slots (mu) central-body
	(let* ((rm (norme2 u))
	       (e (/ (- (* 2 (norme2 duds)) mu) rm))
	       (r (spin (first *iframe*) u))
	       (a (nbody-accel s x p)))
	  (make-instance
	   'spinor-state
	   :u duds
	   :duds (/ (+ (*g3 a r u) (* e u)) 2)
	   :et rm))))))

(defun recover-rotor-3d (fs es)
  (let* ((esr (apply #'recipbvs es))
	 (psi (apply #'+ 1 (mapcar #'(lambda (f er) (*g f er)) fs esr))))
    (unitg psi)))

(defun recover-spinor-3d (rm fs es)
  "Recover a spinor in 3D, given radius, new frame, and old frame"
  (* (recover-rotor-3d fs es) (sqrt rm)))

(defmethod rv-frame (r v)
  "Create a frame from the position and velocity"
  (let ((h (*o r v)))
    (list (unitg r)
	  (unitg (*i r h))
	  (unitg (dual h)))))

(defmethod to-spinor-state (et (x cartesian-state) &key (iframe *iframe*))
  "Given ephemeris time (seconds past J2000), cartesian state, and cartesian problem, return spinor state."
  (with-slots (r v) x
    (let ((u (recover-spinor-3d (norme r) (rv-frame r v) *iframe*)))
      (values
       (make-instance
	'spinor-state
	:u u
	:duds (/ (*g3 v u (first *iframe*)) 2)
	:et et)
       0))))

(defmethod to-cartesian-state (s (x spinor-state) &key (iframe *iframe*))
  (with-slots (u duds et) x
    (values
     (make-instance
      'cartesian-state
      :r (spin (first iframe) u)
      :v (graden (* 2 (*g3 (/ duds (norme2 u)) (ve3 :e1 1) (revg u))) 1))
     et)))

(defclass spinor-problem (cartesian-problem)
  ((x0s :initarg :x0s)
   (s0 :initarg :s0 :initform 0)
   (sfmax :initarg :sfmax)
   (stopfn :initarg :stopfn)
   (stopval :initarg :stopval)
   (stoptest :initarg :stoptest)
   (stoptol :initarg :stoptol)))

(defmethod to-spinor-problem ((p cartesian-problem))
  (with-slots (central-body nbodies spk lsk ref et0 etf eom x0 hmin hmax tol) p
    (make-instance
     'spinor-problem
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
     :tol (/ tol 100)
     ;; New spinor slots
     :x0s (to-spinor-state et0 x0 :iframe *iframe*)
     :s0 0
     :sfmax (/ pi 10)
     :stopfn #'(lambda (s x p) (slot-value x 'et))
     :stopval etf
     :stoptest #'>=
     :stoptol 1d-6)))

(defmethod propagate ((p spinor-problem))
  (with-slots (eom s0 sfmax x0s stopfn stopval stoptest stoptol tol central-body spk lsk) p
    ;; Pre-load SPK and LSK kernels
    (dolist (k spk) (unless (kernelp k) (furnsh k)))
    (unless (kernelp lsk) (furnsh lsk))
    ;; Integrate
    (unwind-protect
	 (rka-stop-nr
	  eom s0 sfmax x0s
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
      
(defmethod to-cartesian-results (results)
  (loop for (s xs) in results
     collect (reverse (multiple-value-list (to-cartesian-state s xs :iframe *iframe*)))))
