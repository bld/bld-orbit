(in-package :bld-orbit2)

(def-unbound rm2 (x cartstate)
  (norme2 (r x)))

(def-unbound rm (x cartstate)
  (sqrt (rm2 x)))

(def-unbound ru (x cartstate)
  (/ (r x) (rm x)))

(def-unbound g (x cartstate)
  (funcall (gfun x) x))

(def-unbound p (x cartstate)
  (solar-pressure x))

(def-unbound bframe (x cartstate)
  (funcall (bfun x) x))

(def-unbound rb (x cartstate)
  (recoverrotor (bframe x) (iframe x)))

(def-unbound rp (x cartstate)
  (funcall (pfun x) x))

(def-unbound pframe (x cartstate)
  (new-frame (rp x) (iframe x)))

(def-unbound n (x cartstate)
  (funcall (nfun x) x))

(def-unbound a (x cartstate)
  (funcall (afun x) x))

