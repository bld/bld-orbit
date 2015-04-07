;;;; Spinor equations of motion

(in-package :bld-orbit)

(defclass spinor-state ()
  ((u :initarg :u)
   (duds :initarg :duds)
   (utc :initarg :utc)))

(defmethod print-object ((x spinor-state) stream)
  (with-slots (u duds utc) x
    (format stream "#<SPINOR-STATE :U ~a :DUDS ~a :UTC ~a>" u duds utc)))

(defstatearithmetic spinor-state (u duds utc))

(defparameter *iframe*
  (list (ve3 :e1 1)
	(ve3 :e2 1)
	(ve3 :e3 1)))

(defmethod eom (s (x spinor-state) p)
  "Keplerian spinor equations of motion given independent variable S, spinor state, and spinor problem"
  (with-slots (u duds utc) x
    (with-slots (central-body) p
      (with-slots (mu) central-body
	(let* ((rm (norme2 u))
	       (e (/ (- (* 2 (norme2 duds)) mu) rm)))
	  (make-instance
	   'spinor-state
	   :u duds
	   :duds (/ (* e u) 2)
	   :utc rm))))))

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

(defmethod to-spinor-state (utc (x cartesian-state) &key (iframe *iframe*))
  "Given time (UTC), cartesian state, and cartesian problem, return spinor state."
  (with-slots (r v) x
    (let ((u (recover-spinor-3d (norme r) (rv-frame r v) *iframe*)))
      (values
       (make-instance
	'spinor-state
	:u u
	:duds (/ (*g3 v u (first *iframe*)) 2)
	:utc utc)
       0))))

(defmethod to-cartesian-state (s (x spinor-state) &key (iframe *iframe*))
  (with-slots (u duds utc) x
    (values
     (make-instance
      'cartesian-state
      :r (spin (first iframe) u)
      :v (graden (* 2 (*g3 (/ duds (norme2 u)) (ve3 :e1 1) (revg u))) 1))
     utc)))

(defclass spinor-problem (cartesian-problem)
  ((x0s :initarg :x0s)
   (s0 :initarg :s0 :initform 0)
   (sfmax :initarg :sfmax)
   (stopfn :initarg :stopfn)
   (stopval :initarg :stopval)
   (stoptest :initarg :stoptest)
   (stoptol :initarg :stoptol)))

(defmethod to-spinor-problem ((p cartesian-problem))
  (with-slots (central-body ref utc0 utcf eom x0 hmin hmax tol) p
    (make-instance
     'spinor-problem
     ;; Cartesian problem slots
     :central-body central-body
     :ref ref
     :utc0 utc0
     :utcf utcf
     :eom eom
     :x0 x0
     :hmin hmin
     :hmax hmax
     :tol (/ tol 100)
     ;; New spinor slots
     :x0s (to-spinor-state utc0 x0 :iframe *iframe*)
     :s0 0
     :sfmax (/ pi 10)
     :stopfn #'(lambda (s x p) (slot-value x 'utc))
     :stopval utcf
     :stoptest #'>=
     :stoptol 1d-6)))

(defmethod propagate ((p spinor-problem))
  (with-slots (eom s0 sfmax x0s stopfn stopval stoptest stoptol tol central-body) p
    (with-kernel (slot-value central-body 'ephemeris)
      (rka-stop-nr
       eom s0 sfmax x0s
       :param p
       :stopfn stopfn
       :stopval stopval
       :stoptest stoptest
       :stoptol stoptol
       :tol tol))))

(defmethod to-cartesian-results (results)
  (loop for (s xs) in results
     collect (reverse (multiple-value-list (to-cartesian-state s xs :iframe *iframe*)))))

