(in-package :bld-orbit)

(defclass cartstate ()
  ;; State slots
  ((r :initarg :r :accessor r :documentation "Position vector")
   (v :initarg :v :accessor v :documentation "Velocity vector")
   (derived :accessor derived :initform (make-hash-table) :documentation "Derived parameters"))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a>" (r x) (v x)))

(defstatearithmetic cartstate (r v))

(defmethod energy-cb ((x cartstate) cb)
  (with-slots (mu) cb
    (with-slots (r v) x
      (- (/ (scalar (exptg v 2)) 2) (/ mu (norme r))))))

(defmethod energy ((x cartstate) sc)
  (with-slots (cb) sc
    (energy-cb x cb)))

(defmethod time-of (s (x cartstate))
  s)

(defmethod eom (tm (x cartstate) sc)
  "Solar sail cartesian equations of motion"
  (with-slots (r v) x
    (with-slots (cb accfun) sc
      (with-slots (mu) cb
	(make-instance
	 'cartstate
	 :r v
	 :v (+ (funcall accfun tm x sc)
	       (gravity tm r mu)))))))

(defmethod to-cartesian ((x cartstate) tm sc)
  (values x tm))

(defderived r-cb (tm x sc)
  (with-slots (cb) sc
    (slot-value (position-velocity cb (time-of tm x)) 'r)))

(defderived r-sc (tm x sc)
  (+ (r-cb tm x sc)
     (r x)))

(defderived r-sun (tm x sc)
  (with-slots (cb sun) sc
    (if (eq sun cb) (r-cb tm x sc)
	(slot-value (position-velocity sun (time-of tm x)) 'r))))

(defderived rsc-sun (tm x sc)
  (- (r-sc tm x sc)
     (r-sun tm x sc)))

(defderived sframe (tm x sc)
  (with-slots (pointfun) sc
    (funcall pointfun tm x sc)))

(defderived a (tm x sc)
  (with-slots (accfun sun) sc
    (with-slots (accfun) sc
      (funcall accfun 
	       (rsc-sun tm x sc)
	       (sframe tm x sc)
	       sc
	       sun))))

(defderived g (tm (x cartstate) sc)
  (with-slots (r) x
    (with-slots (gfun cb) sc
      (funcall gfun r cb))))

(defmethod eom3 (tm (x cartstate) sc)
  (with-slots (v) x
    (make-instance 'cartstate :r v :v (+ (a tm x sc) (g tm x sc)))))

(defparameter *eom3-examples*
  (make-hash
   :kepler
   (let ((t0 (coerce (encode-universal-time 0 0 0 27 2 2014 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 27 2 2015 0) 'double-float)))
     (make-instance
      'sail
      :eom #'eom3
      :cb *sun*
      :sun *sun*
      :gfun #'gravity2
      :accfun #'no-sail
      :pointfun #'(lambda (tm x sc) *j2000*)
      :lightness 0d0
      :area 1260d-6
      :mass 52d0
      :optical nil
      :basis *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (re3 :s 1)
      :outfile "eom3-example-kepler.dat"))
   :normal
   (let ((t0 (coerce (encode-universal-time 0 0 0 27 2 2014 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 27 2 2015 0) 'double-float)))
     (make-instance
      'sail
      :eom #'eom3
      :cb *sun*
      :sun *sun*
      :gfun #'gravity2
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-normal
      :lightness 0d0
      :area 1260d-6
      :mass 52d0
      :optical nil
      :basis *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (re3 :s 1)
      :outfile "eom3-example-normal.dat"))
   :fixed
   (let ((t0 (coerce (encode-universal-time 0 0 0 27 2 2014 0) 'double-float))
	 (tf (coerce (encode-universal-time 0 0 0 27 2 2015 0) 'double-float)))
     (make-instance
      'sail
      :eom #'eom3
      :cb *sun*
      :sun *sun*
      :gfun #'gravity2
      :accfun #'sail-flat-ideal-acc
      :pointfun #'sail-frame-sun-fixed
      :lightness 0d0
      :area 1260d-6
      :mass 52d0
      :optical nil
      :basis *j2000*
      :t0 t0
      :tf tf
      :x0 (position-velocity *earth* t0)
      :rs (rotor (bve3 :e1e2 1) (atan (/ (sqrt 2))))
      :outfile "eom3-example-fixed.dat"))))
