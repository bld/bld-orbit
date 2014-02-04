;;; Orbital state classes

(in-package :bld-orbit2)

(defmethod norminfx ((x g))
  "Infinity norm of geometric algebra state variable"
  (norminf x))

;;; Cartesian state

(defclass cartstate ()
  (;; State slots
   (r :initarg :r :accessor r :documentation "Position vector")
   (v :initarg :v :accessor v :documentation "Velocity vector")
   ;; Parameter slots
   (sc :initarg :sc :accessor sc :documentation "Spacecraft data")
   ;; Derived value slots
   (rm2 :accessor rm2 :documentation "Radius squared")
   (rm :accessor rm :documentation "Position vector magnitude")
   (ru :accessor ru :documentation "Position unit vector")
   (g :accessor g :documentation "Gravitational acceleration")
   (p :accessor p :documentation "Solar pressure")
   (bframe :accessor bframe :documentation "Basis frame")
   (rb :accessor rb :documentation "Rotor of basis frame")
   (pframe :accessor pframe :documentation "Pointing frame")
   (rp :accessor rp :documentation "Rotor of s/c pointing frame")
   (a :accessor a :documentation "External acceleration"))
  (:documentation "Cartesian coordinate state"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a :sc ~a>" (r x) (v x) (sc x)))

(defstatearithmetic cartstate (r v) :oslots (sc))

(defmethod initialize-instance :after ((x cartstate) &key)
  (with-slots (r sc rm2 rm ru g p bframe rb pframe rp a) x
    (setf rm2 (norme2 r))
    (setf rm (sqrt rm2))
    (setf ru (unitg r))
    (setf g (funcall (gfun sc) x))
    (setf p (solar-pressure x))
    (setf bframe (funcall (bfun sc) x))
    (setf rb (recoverrotor bframe (iframe sc)))
    (multiple-value-bind (pframe-tmp rp-tmp) (funcall (pfun sc) x)
      (setf pframe pframe-tmp)
      (setf rp rp-tmp))
    (setf a (funcall (afun sc) x))))
