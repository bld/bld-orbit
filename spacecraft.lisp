;;; Spacecraft data classes

(in-package :bld-orbit2)

(defclass sc ()
  ((cb :initarg :cb :accessor cb :documentation "Central body")
   (gfun :initarg :gfun :accessor gfun :documentation "Gravity function")))

(defclass sail (sc)
  ((area :initarg :area :accessor area)
   (mass :initarg :mass :accessor mass)
   (pfun :initarg :pfun :accessor pfun :documentation "Pointing function")
   (afun :initarg :afun :accessor afun :documentation "Acceleration function")))
