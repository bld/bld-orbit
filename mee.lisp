(in-package :bld-orbit)

;;; Modified equinoctial orbital elements

(defclass meestate ()
  ((p :initarg :p :documentation "Semi latus rectum")
   (f :initarg :f :documentation "e cos (AOP + LAN)")
   (g :initarg :g :documentation "e sin (AOP + LAN)")
   (h :initarg :h :documentation "tan (i/2) cos (LAN)")
   (k :initarg :k :documentation "tan (i/2) sin (LAN)")
   (l :initarg :l :documentation "True longitude: AOP + LAN + true anomaly"))
  (:documentation "Equinoctial orbital element state"))

(defstatearithmetic meestate (p f g h k l))

(let ((slots '(p f g h k l)))
  (defmethod print-object ((x meestate) s)
    (format s "#<MEESTATE ~{~a~^ ~}"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value x slot))))))

(defgeneric to-mee (x s sc))

(defmethod to-mee ((x coestate) truan sc)
  (with-slots (a e i lan aop lan tm) x
    (values
     (make-instance
      'meestate
      :p (* a (- 1 (expt e 2)))
      :f (* e (cos (+ aop lan)))
      :g (* e (sin (+ aop lan)))
      :h (* (tan (/ i 2)) (cos lan))
      :k (* (tan (/ i 2)) (sin lan))
      :l (+ lan aop truan))
     tm)))

(defmethod to-cartesian ((x meestate) tm sc)
  "Modified equinoctial element state to cartesian"
  (with-slots (p f g h k l) x
    (with-slots (cb basis) sc
      (with-slots (mu) cb
	(let* ((alpha2 (- (expt h 2) (expt k 2)))
	       (s2 (+ 1 (expt h 2) (expt k 2)))
	       (w (+ 1 (* f (cos l)) (* g (sin l))))
	       (r (/ p w)))
	  (values
	   (make-instance
	    'cartstate
	    :r (make-vector
		(list
		 (* (/ r s2)
		    (+ (cos l) (* alpha2 (cos l)) (* 2 h k (sin l))))
		 (* (/ r s2)
		    (+ (sin l) (- (* alpha2 (sin l))) (* 2 h k (cos l))))
		 (* (/ (* 2 r) s2)
		    (- (* h (sin l)) (* k (cos l)))))
		basis)
	    :v (make-vector
		(list
		 (* (/ -1 s2)
		    (sqrt (/ mu p))
		    (+ (sin l)
		       (* alpha2 (sin l))
		       (* -2 h k (cos l))
		       g
		       (* -2 f h k)
		       (* alpha2 g)))
		 (* (/ -1 s2)
		    (sqrt (/ mu p))
		    (+ (- (cos l))
		       (* alpha2 (cos l))
		       (* 2 h k (sin l))
		       (- f)
		       (* 2 g h k)
		       (* alpha2 f)))
		 (* (/ 2 s2)
		    (sqrt (/ mu p))
		    (+ (* h (cos l))
		       (* k (sin l))
		       (* f h)
		       (* g k))))
		basis))
	   tm))))))

(defmethod to-coe ((x meestate) tm sc)
  "Classical orbital element state from modified equinoctial"
  (with-slots (p f g h k l) x
    (values
     (make-instance
      'coestate
      :a (/ p (- 1 (expt f 2) (expt g 2)))
      :e (sqrt (+ (expt f 2) (expt g 2)))
      :i (atan (* 2 (sqrt (+ (expt h 2) (expt k 2))))
	       (- 1 (expt h 2) (expt k 2)))
      :aop (atan (- (* g h) (* f k))
		 (+ (* f h) (* g k)))
      :lan (atan k h))
     (- l lan aop))))


(defmethod eom (tm (x meestate) sc)
  "Modified equinoctial orbital elements equations of motion"
  (with-slots (p f g h k l) x
    (with-slots (cb basis accfun) sc
      (with-slots (mu) cb
	(with-slots (r v) (to-cartesian x tm sc)
	  (let* ((rvbasis (rvbasis r v))
		 (fv (funcall accfun tm x sc))
		 (fr (scalar (*i fv (first rvbasis))))
		 (ft (scalar (*i fv (second basis))))
		 (fn (scalar (*i fv (third basis))))
		 (w (+ 1 (* f (cos l)) (* g (sin l))))
		 (s2 (+ 1 (expt h 2) (expt k 2))))
	    (make-instance
	     'meestate
	     :p (* 2 p (/ w) (sqrt (/ p mu)) ft)
	     :f (* (sqrt (/ p mu))
		   (+ (* fr (sin l))
		      (* (+ (* (+ w 1) (cos l)) f) (/ ft w))
		      (* (- (* k (cos l)) (* h (sin l))) g fn (/ w))))
	     :g (* (sqrt (/ p mu))
		   (+ (- (* fr (cos l)))
		      (* (+ (* (+ w 1) (sin l)) g) (/ ft w))
		      (- (* (- (* h (sin l)) (* k (cos l))) (/ (* g fn) w)))))
	     :h (* (sqrt (/ p mu))
		   s2 fn (cos l) (/ (* 2 w)))
	     :k (* (sqrt (/ p mu))
		   s2 fn (sin l) (/ (* 2 w)))
	     :l (+ (* (sqrt (* mu p)) (expt (/ w p) 2))
		   (* (/ w) (sqrt (/ p mu)) (- (* h (sin l)) (* k (cos l))) fn))
	     )))))))
