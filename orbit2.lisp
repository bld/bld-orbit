(in-package :bld-orbit)

;;; Rotor, spinor and basis functions

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

;;; Closure 

(defclass cartesian-state ()
  ((r :initarg :r :documentation "Position vector, 2D or 3D Euclidean")
   (v :initarg :v :documentation "Velocity vector, 2D or 3D Euclidean")))

(let ((mu 1d0)
      (lightness 0.1d0)
      (attitudefun #'(lambda (s x) 
  (defun eom (s x)
    (with-slots (r v) x
      (let ((ruv (unitg r))
	    (rmag (norme r))
	    
      (make-instance
       'cartesian-state
       :r v
       :v 
