(in-package :bld-orbit2)

(defparameter *J2000* (list (ve3 :e1 1d0)
			    (ve3 :e2 1d0)
			    (ve3 :e3 1d0)) "J2000 inertial frame")

(defun orbit-frame (x)
  "Orbit reference from from sail-sun vector and inertial reference frame"
  (let* ((o3 (unitg (r x)))
	 (o1 (unitg (*x (third (iframe (sc x))) o3)))
	 (o2 (*x o3 o1)))
    (list o1 o2 o3)))

(defun new-frame (rotor frame)
  "New reference frame from rotor & original frame"
  (mapcar #'(lambda (basisvector) (rotateg basisvector rotor)) frame))

(defun make-vector (coefs basis)
  "Make a vector given a list of coefficients and a list of basis vectors"
  (apply #'+
	 (loop for coef in coefs
	    for bv in basis
	    collect (* coef bv))))

(defun recoverrotor3d (fs es)
  "Recover a basis given new and original basis vectors"
  (let* ((esr (apply #'recipbvs es))
	 (psi (apply #'+ 1 (mapcar #'(lambda (f er) (*g f er)) fs esr))))
    (unitg psi)))

(defun recoverrotor2d (fs es)
  "Recover a basis given new and original basis vectors"
  (let* ((esr (apply #'recipbvs es))
	 (psi (apply #'+ 2 (mapcar #'(lambda (f er) (*g f er)) fs esr))))
    (unitg psi)))

(defun recoverrotor (fs es)
  "Recover rotor from 2D or 3D basis frames given new and original basis vectors"
  (cond
    ((apply #'= 3 (mapcar #'dimension (append fs es))) (recoverrotor3d fs es))
    ((apply #'= 2 (mapcar #'dimension (append fs es))) (recoverrotor2d fs es))
    (t (error "FS and ES must all be dimension 2 or 3"))))

(defun recoverspinor (r fs es)
  "Recover a spinor given orbit radius, new basis vectors, and original basis vectors"
  (* (recoverrotor fs es) (sqrt r)))

(defun rv-frame (x)
  "Form a basis from the position and velocity vectors. First is the unit position vector. Second is the complement of the first in the orbit plane. Third (if it exists) is the orbit plane normal vector."
  (let* ((r (r x))
	 (v (v x))
	 (h (unitg (*o r v)))
	 (x (unitg r))
	 (y (unitg (*i r h))))
    (cond
      ((= (dimension r) (dimension v) 2) (list x y))
      ((= (dimension r) (dimension v) 3)
       (list x y (dual h)))
      (t (error "R and V must be 2 or 3 dimensional")))))
