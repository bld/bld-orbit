(in-package :bld-orbit2)

(defmethod eom (tm (x cartstate) &optional p)
  (with-slots (v sc a g) x
    (make-instance
     'cartstate
     :r v
     :v (+ a g)
     :sc sc
     :tm tm
     :derivp t)))
