(in-package :bld-orbit2)

(defclass body ()
  ((mu :initarg :mu :accessor mu :documentation "Gravitational parameter (m^3/s^2)")
   (ls :initarg :ls :accessor ls :documentation "Luminosity (W/m^2) at 1 AU")))

(defparameter *sun* 
  (make-instance 
   'body 
   :mu 1.32712440018d20
   :ls 1361.3477d0))
