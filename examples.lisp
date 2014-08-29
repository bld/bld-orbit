(in-package :bld-orbit)

(defparameter *examples*
  (make-hash
   :kepler
   (let* ((t0 (coerce (encode-universal-time 0 0 0 27 2 2014 0) 'double-float))
	  (tf (coerce (encode-universal-time 0 0 0 27 2 2015 0) 'double-float))
	  (sc (make-instance
	       'sc
	       :eom #'eom
	       :cb *ssb*
	       :gfun #'gravity
	       :afun #'no-sail
	       :iframe *j2000*
	       :t0 t0
	       :tf tf
	       :x0 (position-velocity *earth* t0)
	       :outfile "example-kepler.dat")))
     (with-slots (x0) sc
       (with-slots (tm (sc-x0 sc)) x0
	 (setf tm t0)
	 (setf sc-x0 sc)
	 sc)))
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

(defparameter *examples-spinor*
  (make-hash
   :kepler
   (lethash (kepler) *examples*
     (with-slots (eom cb afun iframe t0 x0) kepler
       (make-instance
	'sc
	:eom eom
	:cb cb
	:afun afun
	:iframe iframe
	:t0 0
	:tf (/ pi 14)
	:x0 (to-spinor x0 t0 kepler)
	:outfile "example-kepler-spinor.dat")))
   :normal
   (lethash (normal) *examples*
     (with-slots (eom cb sun afun spfun area mass iframe t0 x0) normal
       (make-instance
	'sc
	:eom eom
	:cb cb
	:sun sun
	:afun afun
	:spfun spfun
	:area area
	:mass mass
	:iframe iframe
	:t0 0
	:tf (/ pi 14)
	:x0 (to-spinor x0 t0 normal)
	:outfile "example-normal-spinor.dat")))
   :fixed
   (lethash (fixed) *examples*
     (with-slots (eom cb sun afun spfun pfun ofun nfun area mass iframe t0 x0 rs) fixed
       (make-instance
	'sc
	:eom eom
	:cb cb
	:sun sun
	:afun afun
	:spfun spfun
	:pfun pfun
	:ofun ofun
	:nfun nfun
	:area area
	:mass mass
	:iframe iframe
	:t0 0
	:tf (/ pi 14)
	:x0 (to-spinor x0 t0 fixed)
	:rs rs
	:outfile "example-fixed-spinor.dat")))))

(defparameter *examples-ks*
  (make-hash
   :kepler
   (lethash (kepler) *examples-spinor*
     (with-slots (eom cb sun afun iframe t0 tf x0) kepler
       (make-instance
	'sc
	:eom eom
	:cb cb
	:afun afun
	:iframe iframe
	:t0 t0
	:tf tf
	:x0 (to-ks x0 t0 kepler)
	:outfile "example-kepler-ks.dat")))
   :normal
   (lethash (normal) *examples-spinor*
     (with-slots (eom cb sun afun spfun area mass iframe t0 tf x0) normal
       (make-instance
	'sc
	:eom eom
	:cb cb
	:sun sun
	:afun afun
	:spfun spfun
	:area area
	:mass mass
	:iframe iframe
	:t0 t0
	:tf tf
	:x0 (to-ks x0 t0 normal)
	:outfile "example-normal-ks.dat")))
   :fixed
   (lethash (fixed) *examples-spinor*
     (with-slots (eom cb sun afun spfun pfun ofun nfun area mass iframe t0 tf x0 rs) fixed
       (make-instance
	'sc
	:eom eom
	:cb cb
	:sun sun
	:afun afun
	:spfun spfun
	:pfun pfun
	:ofun ofun
	:nfun nfun
	:area area
	:mass mass
	:iframe iframe
	:t0 t0
	:tf tf
	:x0 (to-ks x0 t0 fixed)
	:rs rs
	:outfile "example-fixed-ks.dat")))))

