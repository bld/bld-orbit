;;; Orbital state classes

(in-package :bld-orbit2)

(defmethod norminfx ((x g))
  "Infinity norm of geometric algebra state variable"
  (norminf x))

(defmacro def-unbound (slot (var class) &body body)
  "Define slot reader function that executes BODY to establish the value of the slot for later reads"
  `(defmethod ,slot ((,var ,class))
     (if (slot-boundp ,var ',slot)
	 (slot-value ,var ',slot)
	 (setf (slot-value ,var ',slot)
	       ,@body))))

(defmacro def-derived (var class &rest defs)
  "Define derived slot reader functions, establishing values when first called"
  `(progn
     ,@(loop for def in defs
	  collect `(def-unbound ,(first def) (,var ,class) ,@(rest def)))))
