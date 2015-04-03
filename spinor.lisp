(in-package :bld-orbit)

(defclass spinor-state ()
  ((u :initarg :u)
   (duds :initarg :duds)
   (utc :initarg :utc)))

(defmethod print-object ((x spinor-state) stream)
  (with-slots (u duds utc) x
    (format stream "#<SPINOR-STATE :U ~a :DUDS ~a :UTC ~a>" u duds utc)))

(defstatearithmetic spinor-state (u duds))

(defclass spinor-problem ()
  (central-body
   ref
   utc0
   utcf
   eom
   x0
   hmin
   hmax
   tol
   iframe))

(defmethod eom (s (x spinor-state) (p spinor-problem))
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
	   :tm rm))))))

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

(defmethod to-spinor (utc (x cartesian-state) &key problem)
  "Given time (UTC), cartesian state, and cartesian problem, return spinor state."
  (with-slots (iframe) problem
    (with-slots (r v) x
      (let ((u (recover-spinor-3d (norme2 r) (rv-frame r v) iframe)))
	(values
	 (make-instance
	  'spinor-state
	  :u u
	  :duds (/ (*g3 v u (first iframe)) 2)
	  :utc utc)
	 0)))))

(defmethod to-cartesian (s (x spinor-state) &key problem)
  (with-slots (u duds utc) x
    (with-slots (iframe) problem
      (values
       (make-instance
	'cartesian-state
	:r (spin (first iframe) u)
	:v (graden (* 2 (*g3 (/ duds (norme2 u)) (ve3 :e1 1) (revg u))) 1))
       utc))))
