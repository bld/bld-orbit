(in-package :bld-orbit)

(defclass cartstate ()
  ;; State slots
  ((r :initarg :r :accessor r :documentation "Position vector")
   (v :initarg :v :accessor v :documentation "Velocity vector")
   (derived :accessor derived :initform (make-hash-table) :documentation "Derived parameters")
   (tm :initarg :tm :accessor tm :initform nil :documentation "Time")
   (sc :initarg :sc :accessor sc :initform nil :documentation "Spacecraft data"))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a>" (r x) (v x)))

(defstatearithmetic cartstate (r v) :oslots (tm sc))

(defmethod energy-cb ((x cartstate) cb)
  (with-slots (mu) cb
    (with-slots (r v) x
      (- (/ (scalar (exptg v 2)) 2) (/ mu (norme r))))))

(defmethod energy ((x cartstate) sc)
  (with-slots (cb) sc
    (energy-cb x cb)))

(defmethod time-of (s (x cartstate))
  s)

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

(defderived ru (tm (x cartstate) sc)
  (unitg (slot-value x 'r)))

(defderived rm2 (tm (x cartstate) sc)
  (norme2 (r x)))

(defderived rm (tm (x cartstate) sc)
  (sqrt (rm2 tm x sc)))

(defderived p (tm x sc)
  (funcall (slot-value sc 'spfun) tm x sc))

(defderived g (tm (x cartstate) sc)
  (funcall (slot-value sc 'gfun) tm x sc))

(defderived n (tm (x cartstate) sc)
  (funcall (slot-value sc 'nfun) tm x sc))

(defun n-normal (tm (x cartstate) sc)
  (ru tm x sc))

(defun n-fixed (tm (x cartstate) sc)
  

(defun gravity3 (tm (x cartstate) sc)
  (with-slots (cb) sc
    (with-slots (mu) cb
      (- (* (/ mu (rm2 tm x sc)) 
	    (ru tm x sc))))))

(defun sp-inverse-square (tm (x cartstate) sc)
  "Solar pressure inverse square function"
  (with-slots (sun) sc
    (with-slots (ls) sun
      (/ (* ls (expt (/ *au* (rm tm x sc)) 2)) *c*))))

(defun sail-ideal-acc2 (tm (x cartstate) sc)
  (with-slots (area mass) sc
    (let ((p (solar-pressure
    (/ (* 2 p area (expt (scalar (*i ru n)) 2) n) mass))

(defmethod eom3 (tm (x cartstate) sc)
  (with-slots (v) x
    (make-instance 
     'cartstate 
     :r v 
     :v (+ (a tm x sc) (g tm x sc)))))

;; Old stuff

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

