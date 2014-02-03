(in-package :bld-orbit2)

(defparameter *eg-sc*
  (make-instance 
   'sc 
   :cb *sun* 
   :gfun #'gravity))
