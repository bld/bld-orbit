;;; Generalized Model for Solar Sails

(in-package :bld-orbit)

(defclass material ()
  ((r :initarg :r)
   (s :initarg :s)
   (bf :initarg :bf)
   (ef :initarg :ef)
   (bb :initarg :bb)
   (eb :initarg :eb)))

(defparameter *material*
  (make-instance
   'material
   :r 0.91d0
   :s 0.94d0
   :bf 0.79d0
   :ef 0.025d0
   :bb 0.67d0
   :eb 0.27d0))

(defmethod a1 ((m material))
  (with-slots (r s) m
    (* 2 r s)))

(defmethod a2 ((m material))
  (with-slots (r s bf ef bb eb) m
    (+ (* bf (- 1 s) r)
       (* (- 1 r)
	  (- (* ef bf)
	     (* eb bb))
	  (/ (+ ef eb))))))

(defmethod a3 ((m material))
  (with-slots (r s) m
    (- 1 (* r s))))

(defclass gsm-coef ()
  ((m :initarg :m)
   a1 a2 a3
   (j1 :initarg :j1)
   (j2 :initarg :j2)
   (j3 :initarg :j3)
   j3_m
   (l :initarg :l)
   (k2 :initarg :k2)
   (k3 :initarg :k3)
   k3_m)
  (:documentation "Coefficients for use in the GSM force and moment calculations"))

(defparameter *j3*
  #3a(
      (
       (-4.049695115508793e-06     1.577494758486461e-05     6.139288348204577e-04)
       (1.577494758486461e-05    -2.743513355019449e-06    -5.426699035474121e-05)
       (6.139288348204577e-04    -5.426699035474120e-05    -5.620418135529490e-03))
      (
       (1.577494758486461e-05    -2.743513355019449e-06    -5.426699035474121e-05)
       (-2.743513355019448e-06     1.075726674329548e-05     6.106349796879265e-04)
       (-5.426699035474121e-05     6.106349796879265e-04    -2.398752979682822e-02))
      (
       (6.139288348204577e-04    -5.426699035474120e-05    -5.620418135529490e-03)
       (-5.426699035474121e-05     6.106349796879265e-04    -2.398752979682822e-02)
       (-5.620418135529490e-03    -2.398752979682822e-02     8.583040955919375e+01))))

(defmethod initialize-instance :after ((c gsm-coef) &key)
  (with-slots (m a1 a2 a3 j3 j3_m k3 k3_m) c
    (setf
     ;; Reshape j3 and k3 tensors to 2D
     j3_m (array-operations:reshape j3 '(9 3))
     k3_m (array-operations:reshape k3 '(9 3))
     ;; Calculate derived material properties
     a1 (a1 m)
     a2 (a2 m)
     a3 (a3 m))))

