;;; Use GECO to optimize trajectories

(in-package :bld-orbit)

(defclass rs-times-chromosome (sequence-chromosome)
  ()
  (:documentation "RS lookup table times sequence chromosome"))

(defclass rs-sia-chromosome (sequence-chromosome)
  ()
  (:documentation "RS lookup table SIA sequence chromosome"))

(defmethod size ((self rs-times-chromosome))
  10)

(defmethod size ((self rs-sia-chromosome))
  10)

