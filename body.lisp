(in-package :bld-orbit)

(defclass body ()
  ((name :initarg :name)
   (center :initarg :center)
   (ref :initarg :ref)
   (ephemeris :initarg :ephemeris)
   (mu :initarg :mu)))

(defmethod gravity ((r ve3) (b body))
  "Gravitational acceleration given UTC, position vector of s/c relative to body, and body"
  (with-slots (mu) b
    (- (* (/ mu (expt (norme r) 3)) r))))

