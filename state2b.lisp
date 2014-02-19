(in-package :bld-orbit2)

(def-derived x cartstate
  (rm2 (norme2 (r x)))
  (rm (sqrt (rm2 x)))
  (ru (/ (r x) (rm x)))
  (g (funcall (gfun x) x))
  (p (solar-pressure x))
  (bframe (funcall (bfun x) x))
  (rb (recoverrotor (bframe x) (iframe x)))
  (rp (funcall (pfun x) x))
  (pframe (new-frame (rp x) (iframe x)))
  (n (funcall (nfun x) x))
  (afun (funcall (afun x) x)))
