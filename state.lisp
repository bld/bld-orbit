;;; Orbital state classes

(in-package :bld-orbit)

(defmethod norminfx ((x g))
  "Infinity norm of geometric algebra state variable. Used in state equations of motion."
  (norminf x))

(defmacro defderived (name (ivar statedef pvar) docstring &body form)
  "Define an accessor function to get or calculate (via FORM) & set
derived variables. Also set independant (IVAR) & parameter (PVAR)
variables used in form."
  (let ((var (cond ((symbolp statedef) statedef)
		   ((listp statedef) (first statedef)))))
    `(defmethod ,name (,ivar ,statedef ,pvar) 
       ,docstring
       (aif (gethash ,(make-keyword name) (derived ,var)) it
	    (setf (gethash ,(make-keyword name) (derived ,var))
		  ,@form)))))

(defmacro with-derived (slots ivar statevar pvar &body body)
  `(let ,(mapcar #'(lambda (slot) `(,slot (,slot ,ivar ,statevar ,pvar))) slots)
     ,@body))
