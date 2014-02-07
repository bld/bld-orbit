(in-package :bld-orbit2)

(defclass problem ()
  ((sc :initarg :sc :accessor sc :documentation "Spacecraft data")
   (t0 :initarg :t0 :accessor t0 :documentation "Initial time")
   (tf :initarg :tf :accessor tf :documentation "Final time")
   (x0 :initarg :x0 :accessor x0 :documentation "Initial state")))


(defmethod propagate ((sc sail) &key outfile (hmax-factor 1))
  "Propagate sailcraft trajectory. Default maximum stepsize is the difference between T0 and TF."
  (with-slots (eom t0 tf x0) sc
    (let ((results (rka eom t0 tf x0 :param sc :hmax (/ (- tf t0) hmax-factor))))
      (when outfile (write-cart-traj outfile (to-cartesian-traj results sc)))
      results)))

(defmethod propagate ((h hash-table) &key)
  "Propagate a hash table of objects & write data file"
  (maphash2 #'(lambda (k v) (propagate v)) h))
