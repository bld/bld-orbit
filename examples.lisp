(in-package :bld-orbit)

;;; EOM2 examples

(defparameter *sail-3d-eom2-examples*
  (make-hash*
   :cart-kepler
   (make-instance
    'sail
    :eom #'eom2
    :cb *earth*
    :sun *sun*
    :accfun #'no-sail
    :pointfun #'(lambda (s x sc) *j2000*)
    :lightness 0d0
    :area 1200d0
    :mass 45d0
    :optical nil
    :basis *j2000*
    :t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float)
    :tf (coerce (encode-universal-time 0 0 0 17 12 2013 0) 'double-float)
    :rs (re3 :s 1)
    :outfile "sail-3d-cart-geo-kepler-eg.dat"
    :x0 (let* ((w (/ (* 2d0 pi) 86164d0))
	       (r (expt (/ (slot-value *earth* 'mu) (expt w 2d0)) 1/3))
	       (v (* w r)))
	  (make-instance 
	   'cartstate 
	   :r (* r (ve3 :e1 1))
	   :v (* v (unitg (ve3 :e2 1 :e3 0.1))))))
   :cart-normal
   (make-instance
    'sail
    :eom #'eom2
    :cb *earth*
    :sun *sun*
    :accfun #'sail-flat-ideal-acc
    :pointfun #'sail-frame-sun-normal
    :lightness 0d0
    :area 1200d-6
    :mass 4.5d0
    :optical nil
    :basis *j2000*
    :t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float)
    :tf (coerce (encode-universal-time 0 0 0 17 12 2013 0) 'double-float)
    :rs (re3 :s 1)
    :outfile "sail-3d-cart-geo-normal-eg.dat"
    :x0 (let* ((w (/ (* 2d0 pi) 86164d0))
	       (r (expt (/ (slot-value *earth* 'mu) (expt w 2d0)) 1/3))
	       (v (* w r)))
	  (make-instance 
	   'cartstate 
	   :r (* r (ve3 :e1 1))
	   :v (* v (unitg (ve3 :e2 1 :e3 0.1))))))
   :cart-solar-kepler
   (let ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 16 12 2014 0) 'double-float)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'no-sail
      :pointfun #'sail-frame-sun-normal
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (re3 :s 1)
      :outfile "sail-3d-cart-solar-kepler-eg.dat"))
   :cart-solar-normal
   (let ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 16 12 2014 0) 'double-float)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-normal
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (re3 :s 1)
      :outfile "sail-3d-cart-solar-normal-eg.dat"))
   :cart-solar-fixed
   (let ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 16 12 2014 0) 'double-float)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-fixed
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2))))
      :outfile "sail-3d-cart-solar-fixed-eg.dat"))
   :cart-solar-table
   (let ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	 (t1 (coerce (encode-universal-time 0 0 0 16 6 2014 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 16 12 2015 0) 'double-float)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-table
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (list 
	   (list t1 (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2)))))
	   (list tf (rotor (bve3 :e1e2 1) (- (atan (/ (sqrt 2)))))))
      :outfile "sail-3d-cart-solar-table-eg.dat"))
   :ks-solar-kepler
   (let* ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	  (tf-approx (coerce (encode-universal-time 0 0 0 16 12 2014 0) 'double-float))
	  (sf (/ (- tf-approx t0) *au*)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'no-sail
      :pointfun #'sail-frame-sun-normal
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 0
      :tf sf
      :x0 (to-initial-ks (position-velocity *earth* t0) t0 (make-instance 'sail :basis *j2000* :cb *sun*))
      :rs nil
      :outfile "sail-3d-ks-solar-kepler-eg.dat"))
   :ks-solar-normal
   (let* ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	  (tf-approx (coerce (encode-universal-time 0 0 0 16 12 2014 0) 'double-float))
	  (sf (/ (- tf-approx t0) *au*)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-normal
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 0
      :tf sf
      :x0 (to-initial-ks (position-velocity *earth* t0) t0 (make-instance 'sail :basis *j2000* :cb *sun*))
      :rs nil
      :outfile "sail-3d-ks-solar-normal-eg.dat"))
   :ks-solar-fixed
   (let* ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	  (tf-approx (coerce (encode-universal-time 0 0 0 16 12 2014 0) 'double-float))
	  (sf (/ (- tf-approx t0) *au*)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-fixed
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 0
      :tf sf
      :x0 (to-initial-ks (position-velocity *earth* t0) t0 (make-instance 'sail :basis *j2000* :cb *sun*))
      :rs (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2))))
      :outfile "sail-3d-ks-solar-fixed-eg.dat"))
   :ks-solar-table
   (let* ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float))
	  (tf-approx (coerce (encode-universal-time 0 0 0 16 12 2015 0) 'double-float))
	  (sf (/ (- tf-approx t0) *au*)))
     (make-instance
      'sail
      :eom #'eom2
      :cb *sun*
      :sun *sun*
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-table
      :lightness 0d0
      :area 1200d-6
      :mass 45d0
      :optical nil
      :basis *j2000*
      :t0 0
      :tf sf
      :x0 (to-initial-ks (position-velocity *earth* t0) t0 (make-instance 'sail :basis *j2000* :cb *sun*))
      :rs (list 
	   (list (/ sf 2) (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2)))))
	   (list sf (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2))))))
      :outfile "sail-3d-ks-solar-table-eg.dat"))))
