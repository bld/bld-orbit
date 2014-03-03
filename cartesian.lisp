(in-package :bld-orbit)

(defclass cartstate ()
  ;; State slots
  ((r :initarg :r :accessor r :documentation "Position vector")
   (v :initarg :v :accessor v :documentation "Velocity vector")
   (derived :accessor derived :initform (make-hash-table) :documentation "Derived parameters")
   (sc :initarg :sc)
   (tm :initarg :tm))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a>" (r x) (v x)))

(defstatearithmetic cartstate (r v) :oslots (tm sc))

(defmethod energy-cb ((x cartstate) cb)
  (with-slots (r v) x
    (- (/ (scalar (exptg v 2)) 2) (/ (mu cb) (norme r)))))

(defmethod energy ((x cartstate) sc)
  (energy-cb x (cb sc)))

(defmethod time-of (tm (x cartstate))
  tm)

(defmethod to-cartesian ((x cartstate) tm sc)
  (values x tm))

(defmethod eom (tm (x cartstate) sc)
  "Equations of motion for cartesian state"
  (with-derived (a g) tm x sc
    (with-slots (v (tm-x tm)) x
      (setf tm-x tm)
      (make-instance 
       'cartstate 
       :r v 
       :v (+ a g)
       :tm tm
       :sc sc))))
