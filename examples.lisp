(in-package :bld-orbit)

(defparameter *examples*
  (make-hash
   :kepler
   (let ((t0 (coerce (encode-universal-time 0 0 0 27 2 2014 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 27 2 2015 0) 'double-float)))
     (make-instance
      'sc
      :eom #'eom
      :cb *ssb*
      :sun *sun*
      :gfun #'gravity
      :afun #'no-sail
      :pfun #'oframe
      :iframe *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :outfile "example-kepler.dat"))
   :normal
   (let ((t0 (coerce (encode-universal-time 0 0 0 27 2 2014 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 27 2 2015 0) 'double-float)))
     (make-instance
      'sc
      :eom #'eom
      :cb *ssb*
      :sun *sun*
      :gfun #'gravity
      :afun #'sail-ideal-acc-normal
      :spfun #'solar-pressure
      :area 1260d-6
      :mass 52d0
      :iframe *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (re3 :s 1)
      :outfile "example-normal.dat"))
   :fixed
   (let ((t0 (coerce (encode-universal-time 0 0 0 27 2 2014 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 27 2 2015 0) 'double-float)))
     (make-instance
      'sc
      :eom #'eom
      :cb *ssb*
      :sun *sun*
      :gfun #'gravity
      :afun #'sail-ideal-acc-fixed
      :spfun #'solar-pressure
      :pfun #'fixed
      :ofun #'rv-frame
      :nfun #'first
      :area 1260d-6
      :mass 52d0
      :iframe *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2))))
      :outfile "example-fixed.dat"))))
