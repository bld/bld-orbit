(in-package :sail-lol)

(defmacro defderived (name &body form)
  `(defpan ,name (x) (self derived)
     (aif (gethash ,(make-keyword name) derived) it
	  (setf (gethash ,(make-keyword name) derived) ,@form))))
