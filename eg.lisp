(in-package :bld-orbit2)

(defparameter *eg-sc*
  (make-instance 
   'sc 
   :cb *sun*
   :nb nil
   :gfun #'gravity))

(defparameter *eg-sail*
  (make-instance
   'sail
   :cb *sun*
   :gfun #'gravity
   :area 1200
   :mass 50
   :pfun #'fixed
   :afun #'ideal-sail-acc
   :bfun #'orbit-frame
   :iframe *j2000*
   :rs (rotor (bve3 :e1e3 -1) (atan (sqrt (/ 2))))))
