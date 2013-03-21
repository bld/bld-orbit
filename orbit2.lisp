(in-package :bld-orbit)

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

(defclass cartesian-state ()
  ((r :initarg :r :documentation "Position vector, 2D or 3D Euclidean")
   (v :initarg :v :documentation "Velocity vector, 2D or 3D Euclidean")))

(defmethod print-object ((x cartesian-state) stream)
  (with-slots (r v) x
    (format stream "#<CARTESIAN-STATE :R ~a :V ~a>" r v)))

(defstatearithmetic cartesian-state (r v))

(defmethod norminfx ((x g))
  (norminf x))

(defun sail-attitude-normal (s r v data)
  (unitg r))

(defun sail-force-ideal (s r v data)
  (with-slots (mu lightness attitudefun) data
    (let ((n (funcall attitudefun s r v data)))
      (* lightness mu (/ (norme2 r))
	 (expt (scalar (*i (unitg r) n)) 2)
	 n))))

(defun sail-force-null (s r v data)
  (* 0d0 r))

(defun gravity-inverse-square (r data)
  (with-slots (mu) data
    (- (* (/ mu (norme2 r)) (unitg r)))))

(defclass sail-data ()
  ((mu :initarg :mu :accessor mu)
   (lightness :initarg :lightness :accessor lightness)
   (attitudefun :initarg :attitudefun :accessor attitudefun)
   (sailforcefun :initarg :sailforcefun :accessor sailforcefun)
   (gravityfun :initarg :gravityfun :accessor gravityfun)
   (x0 :initarg :x0 :accessor x0)
   (s0 :initarg :s0 :accessor s0)
   (sf :initarg :sf :accessor sf)
   (eom :initarg :eom :accessor eom)))

(defmethod initialize-instance :after ((sail sail-data) &key)
  "Sets the DATA variable in the EOM closure to the sail data object"
  (with-pandoric (data) (eom sail)
    (setf data sail)))

(defun eom-cartesian ()
  "Return closure lambda with data variable referring to a data structure"
  (let (data)
    (plambda (s x) (data)
      (with-slots (r v) x
	(make-instance 
	 'cartesian-state
	 :r v
	 :v (+ (funcall (sailforcefun data) s r v data)
	       (funcall (gravityfun data) r data)))))))

(defparameter *sail-data* 
  (make-instance
   'sail-data
   :mu 1d0
   :lightness 0.1d0
   :attitudefun #'sail-attitude-normal
   :sailforcefun #'sail-force-null
   :gravityfun #'gravity-inverse-square
   :x0 (make-instance 'cartesian-state :r (ve2 :c1 1d0) :v (ve2 :c10 1d0))
   :s0 0d0
   :sf 1d0
   :eom (eom-cartesian)))

(defparameter *sail-data2* 
  (make-instance
   'sail-data
   :mu 1d0
   :lightness 0.2d0
   :attitudefun #'sail-attitude-normal
   :sailforcefun #'sail-force-null
   :gravityfun #'gravity-inverse-square
   :x0 (make-instance 'cartesian-state :r (ve2 :c1 1d0) :v (ve2 :c10 1d0))
   :s0 0d0
   :sf 1d0
   :eom (eom-cartesian)))

(defmethod propagate ((data sail-data))
  (with-slots (eom s0 sf x0) data
    (rka eom s0 sf x0)))
