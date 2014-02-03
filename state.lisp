;;; Orbital state classes

(in-package :bld-orbit2)

(defmethod norminfx ((x g))
  "Infinity norm of geometric algebra state variable"
  (norminf x))

;;; Cartesian state

(defclass cartstate ()
  ((r :initarg :r :accessor r :documentation "Position vector")
   (v :initarg :v :accessor v :documentation "Velocity vector")
   (sc :initarg :sc :accessor sc :documentation "Spacecraft data"))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a :sc ~a>" (r x) (v x) (sc x)))

(defstatearithmetic cartstate (r v) :oslots (sc))
