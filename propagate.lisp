(in-package :bld-orbit)

(defgeneric propagate (sc &key))

(defmethod propagate ((sc sail) &key (outfile (slot-value sc 'outfile)) (hmax-factor 1))
  "Propagate sailcraft trajectory. Default maximum stepsize is the difference between T0 and TF."
  (with-slots (eom t0 tf x0) sc
    (let ((results (rka eom t0 tf x0 :param sc :hmax (/ (- tf t0) hmax-factor))))
      (when outfile (write-cart-traj outfile (to-cartesian-traj results sc)))
      results)))

(defmethod propagate ((h hash-table) &key)
  "Propagate a hash table of objects & write data file"
  (maphash2 #'(lambda (k v) (propagate v)) h))

(defmethod propagate-table ((sc sail) &key (outfile (slot-value sc 'outfile)) (hmax-factor 1))
  "Propagate sailcraft trajectory. Default maximum stepsize is the difference between T0 and TF."
  (with-slots (eom t0 x0 rs-table rs) sc
    (loop for (tfi rsi) in rs-table
       for t0i = t0 then tf
       for tf = tfi
       for x0i = x0 then (second (car (last results)))
       do (setf rs rsi)
       append (rka eom t0i tf x0i
		   :param sc 
		   :hmax (/ (- tf t0i) hmax-factor)) into results
       finally (when outfile 
		 (write-cart-traj outfile (to-cartesian-traj results sc)))
	 (return results))))
