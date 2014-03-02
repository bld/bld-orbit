(in-package :bld-orbit)

;;; Utility functions

(defun new-frame (rotor frame)
  "New reference frame from rotor & original frame"
  (mapcar #'(lambda (basisvector) (rotateg basisvector rotor)) frame))

(defun make-vector (coefs basis)
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
  (cond
    ((apply #'= 3 (mapcar #'dimension (append fs es))) (recoverrotor3d fs es))
    ((apply #'= 2 (mapcar #'dimension (append fs es))) (recoverrotor2d fs es))
    (t (error "FS and ES must all be dimension 2 or 3"))))

(defun recoverspinor (r fs es)
  "Recover a spinor given orbit radius, new basis vectors, and original basis vectors"
  (* (recoverrotor fs es) (sqrt r)))

;;; Pointing functions

(defun fixed (tm x sc)
  "Pointing frame rotates from orbit frame by fixed RS rotor parameter"
  (with-slots (iframe) sc
    (new-frame (rotor-p tm x sc) iframe)))
