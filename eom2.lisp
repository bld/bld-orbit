(in-package :bld-orbit2)

(defmethod eom (tm (x cartstate) &optional p)
  (make-instance
   'cartstate
   :r (v x)
   :v (+ (a x) (g x))
   :cb (cb x)
   :gfun (gfun x)
   :area (area x)
   :mass (mass x)
   :pfun (pfun x)
   :afun (afun x)
   :bfun (bfun x)
   :nfun (nfun x)
   :iframe (iframe x)
   :rs (rs x)
   :tm tm))
