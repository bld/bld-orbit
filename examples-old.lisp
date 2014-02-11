(in-package :bld-orbit)

;;; 2D examples

;; 2D Kepler examples

(defparameter *sail-2d-kepler-examples*
  (make-hash*
   :cart (make-instance 'sail)
   :spin (with-slots (x0) cart
	   (make-instance
	    'sail
	    :tf (* 2 pi)
	    :x0 (to-spinor x0 0 cart)
	    :outfile "sail-2d-spin-kepler-eg.dat"))
   :ks (with-slots (x0) cart
	 (make-instance
	  'sail
	  :tf (* 2 pi)
	  :x0 (to-ks x0 0 cart)
	  :outfile "sail-2d-ks-kepler-eg.dat"))))

;; 2D normal examples

(defparameter *sail-2d-normal-examples*
  (make-hash*
   :cart (make-instance
	  'sail
	  :tf (* 4 pi)
	  :pointfun #'sail-pointing-normal
	  :lightness 0.1
	  :outfile "sail-2d-cart-normal-eg.dat")
   :spin (with-slots (x0 pointfun lightness) cart
	   (make-instance
	    'sail
	    :tf (* 4 pi)
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-spinor x0 0 cart)
	    :outfile "sail-2d-spin-normal-eg.dat"))
   :ks (with-slots (x0 pointfun lightness) cart
	 (make-instance
	  'sail
	  :tf (* 8 pi)
	  :pointfun #'sail-pointing-normal
	  :lightness lightness
	  :x0 (to-ks x0 0 cart)
	  :outfile "sail-2d-ks-normal-eg.dat"))))

;; 2D fixed examples

(defparameter *sail-2d-fixed-examples*
  (make-hash*
   :cart (make-instance
	  'sail
	  :accfun #'sail-ideal-acc
	  :pointfun #'sail-pointing-fixed
	  :lightness 0.1
	  :tf (* 4 pi)
	  :rs (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2))))
	  :outfile "sail-2d-cart-fixed-eg.dat")
   :spin (with-slots (pointfun lightness x0 rs) cart
	   (make-instance
	    'sail
	    :tf (* 4 pi)
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-spinor x0 0 cart)
	    :rs rs
	    :outfile "sail-2d-spin-fixed-eg.dat"))
   :ks (with-slots (pointfun lightness x0 rs) cart
	 (make-instance
	  'sail
	  :tf (* 4 pi)
	  :pointfun pointfun
	  :lightness lightness
	  :x0 (to-ks x0 0 cart)
	  :rs rs
	  :outfile "sail-2d-ks-fixed-eg.dat"))))

;; 2D table examples

(defparameter *sail-2d-table-examples*
  (make-hash*
   :cart (make-instance
	  'sail
	  :lightness 0.1
	  :tf (* 4 pi)
	  :pointfun #'sail-pointing-table
	  :rs (list (list (* 2 pi) (rotor (bve2 :e1e2 1) (/ pi 2)))
		    (list (* 4 pi) (rotor (bve2 :e1e2 1) (atan (/ (sqrt 2))))))
	  :outfile "sail-2d-cart-table-eg.dat")
   :spin (with-slots (pointfun lightness x0 rs) cart
	   (make-instance
	    'sail
	    :tf (* 4 pi)
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-spinor x0 0 cart)
	    :rs rs
	    :outfile "sail-2d-spin-table-eg.dat"))
   :ks (with-slots (pointfun lightness x0 rs) cart
	 (make-instance
	  'sail
	  :tf (* 4 pi)
	  :pointfun pointfun
	  :lightness lightness
	  :x0 (to-ks x0 0 cart)
	  :rs rs
	  :outfile "sail-2d-ks-table-eg.dat"))))

;;; 3D examples

;; 3D Kepler examples

(defparameter *sail-3d-kepler-examples*
  (make-hash*
   :coe (make-instance
	 'sail
	 :tf (* 4 pi)
	 :accfun #'(lambda (s x sc) (ve3))
	 :lightness 0
	 :x0 (make-instance
	      'coestate
	      :a 1.1
	      :e 0.1
	      :i 0.1
	      :lan 0.1
	      :aop 0.1
	      :tm 0)
	 :basis *j2000*
	 :rs (re3 :s 1)
	 :outfile "sail-3d-coe-kepler-eg.dat")
   :cart (with-slots (tf accfun lightness x0 basis rs) coe
	   (make-instance
	    'sail
	    :tf tf
	    :accfun accfun
	    :lightness lightness
	    :x0 (to-cartesian x0 0 coe)
	    :basis basis
	    :rs rs
	    :outfile "sail-3d-cart-kepler-eg.dat"))
   :spin (with-slots (tf accfun lightness x0 basis rs) coe
	   (make-instance
	    'sail
	    :tf tf
	    :accfun accfun
	    :lightness lightness
	    :x0 (to-spinor x0 0 coe)
	    :basis basis
	    :rs rs
	    :outfile "sail-3d-spin-kepler-eg.dat"))
   :ks (with-slots (tf accfun lightness x0 basis rs) coe
	 (make-instance
	  'sail
	  :tf tf
	  :accfun accfun
	  :lightness lightness
	  :x0 (to-ks x0 0 coe)
	  :basis basis
	  :rs rs
	  :outfile "sail-3d-ks-kepler-eg.dat"))
   :mee (with-slots (tf accfun lightness x0 basis rs) coe
	  (make-instance
	   'sail
	   :tf tf
	   :accfun accfun
	   :lightness lightness
	   :x0 (to-mee x0 0 coe)
	   :basis basis
	   :rs rs
	   :outfile "sail-3d-mee-kepler-eg.dat"))))

;; 3D sail normal examples

(defparameter *sail-3d-normal-examples*
  (make-hash*
   :coe (make-instance
	 'sail
	 :tf (* 4 pi)
	 :pointfun #'sail-pointing-normal
	 :lightness 0.1
	 :x0 (make-instance
	      'coestate
	      :a 1.1
	      :e 0.1
	      :i 0.1
	      :lan 0.1
	      :aop 0.1
	      :tm 0)
	 :basis *j2000*
	 :outfile "sail-3d-coe-normal-eg.dat")
   :cart (with-slots (tf pointfun lightness x0 basis) coe
	   (make-instance
	    'sail
	    :tf tf
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-cartesian x0 0 coe)
	    :basis basis
	    :outfile "sail-3d-cart-normal-eg.dat"))
   :spin (with-slots (tf pointfun lightness x0 basis) coe
	   (make-instance
	    'sail
	    :tf tf
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-spinor x0 0 coe)
	    :basis basis
	    :outfile "sail-3d-spin-normal-eg.dat"))
   :ks (with-slots (tf pointfun lightness x0 basis) coe
	 (make-instance
	  'sail
	  :tf tf
	  :pointfun pointfun
	  :lightness lightness
	  :x0 (to-ks x0 0 coe)
	  :basis basis
	  :outfile "sail-3d-ks-normal-eg.dat"))
   :mee (with-slots (tf pointfun lightness x0 basis) coe
	  (make-instance
	   'sail
	   :tf tf
	   :pointfun pointfun
	   :lightness lightness
	   :x0 (to-mee x0 0 coe)
	   :basis basis
	   :outfile "sail-3d-mee-normal-eg.dat"))))

;; 3D sail fixed examples

(defparameter *sail-3d-fixed-examples*
  (make-hash*
   :coe (make-instance
	 'sail
	 :tf (* 4 pi)
	 :pointfun #'sail-pointing-fixed
	 :lightness 0.1d0
	 :x0 (make-instance
	      'coestate
	      :a 1.1
	      :e 0.1
	      :i 0.1
	      :lan 0.1
	      :aop 0.1
	      :tm 0)
	 :basis *j2000*
	 :rs (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2))))
	 :outfile "sail-3d-coe-fixed-eg.dat")
   :cart (with-slots (tf pointfun lightness x0 basis rs) coe
	   (make-instance
	    'sail
	    :tf tf
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-cartesian x0 0 coe)
	    :basis basis
	    :rs rs
	    :outfile "sail-3d-cart-fixed-eg.dat"))
   :spin (with-slots (tf pointfun lightness x0 basis rs) coe
	   (make-instance
	    'sail
	    :tf tf
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-spinor x0 0 coe)
	    :basis basis
	    :rs rs
	    :outfile "sail-3d-spin-fixed-eg.dat"))
   :ks (with-slots (tf pointfun lightness x0 basis rs) coe
	 (make-instance
	  'sail
	  :tf tf
	  :pointfun pointfun
	  :lightness lightness
	  :x0 (to-ks x0 0 coe)
	  :basis basis
	  :rs rs
	  :outfile "sail-3d-ks-fixed-eg.dat"))
   :mee (with-slots (tf pointfun lightness x0 basis rs) coe
	  (make-instance
	   'sail
	   :tf tf
	   :pointfun pointfun
	   :lightness lightness
	   :x0 (to-mee x0 0 coe)
	   :basis basis
	   :rs rs
	   :outfile "sail-3d-mee-fixed-eg.dat"))))

;; 3D table examples

(defparameter *sail-3d-table-examples*
  (make-hash*
   :coe (make-instance
	 'sail
	 :tf (* 4 pi)
	 :pointfun #'sail-pointing-table
	 :lightness 0.1d0
	 :x0 (make-instance
	      'coestate
	      :a 1.1
	      :e 0.1
	      :i 0.1
	      :lan 0.1
	      :aop 0.1
	      :tm 0)
	 :basis *j2000*
	 :rs (list (list (* 2 pi) (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2)))))
		   (list (* 4 pi) (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2)))))
		   )
	 :outfile "sail-3d-coe-table-eg.dat")
   :cart (with-slots (tf pointfun lightness x0 basis rs) coe
	   (make-instance
	    'sail
	    :tf tf
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-cartesian x0 0 coe)
	    :basis basis
	    :rs rs
	    :outfile "sail-3d-cart-table-eg.dat"))
   :spin (with-slots (tf pointfun lightness x0 basis rs) coe
	   (make-instance
	    'sail
	    :tf tf
	    :pointfun pointfun
	    :lightness lightness
	    :x0 (to-spinor x0 0 coe)
	    :basis basis
	    :rs rs
	    :outfile "sail-3d-spin-table-eg.dat"))
   :ks (with-slots (tf pointfun lightness x0 basis rs) coe
	 (make-instance
	  'sail
	  :tf tf
	  :pointfun pointfun
	  :lightness lightness
	  :x0 (to-ks x0 0 coe)
	  :basis basis
	  :rs rs
	  :outfile "sail-3d-ks-table-eg.dat"))
   :mee (with-slots (tf pointfun lightness x0 basis rs) coe
	  (make-instance
	   'sail
	   :tf tf
	   :pointfun pointfun
	   :lightness lightness
	   :x0 (to-mee x0 0 coe)
	   :basis basis
	   :rs rs
	   :outfile "sail-3d-mee-table-eg.dat"))))

