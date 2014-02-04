;;; Spacecraft data classes

(in-package :bld-orbit2)

(defclass sc ()
  ((cb :initarg :cb :accessor cb :documentation "Central body")
   (gfun :initarg :gfun :accessor gfun :documentation "Gravity function")))

(defclass sail (sc)
  ((area :initarg :area :accessor area)
   (mass :initarg :mass :accessor mass)
   (pfun :initarg :pfun :accessor pfun :documentation "Pointing function")
   (afun :initarg :afun :accessor afun :documentation "Acceleration function")
   (bfun :initarg :bfun :accessor bfun :documentation "Basis function for pointing")
   (iframe :initarg :iframe :accessor iframe :documentation "Inertial frame")
   (rs :initarg :rs :accessor rs :documentation "Rotor of sail relative to basis")))

