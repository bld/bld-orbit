(in-package :bld-orbit2)

(defmethod eom (tm (x cartstate) &optional p)
  (with-slots (r v sc) x
    (make-instance 
     'cartstate 
     :r (v x) 
     :v (funcall (gfun sc) x) 
     :sc sc)))
