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
   (cb :initarg :cb :accessor cb :documentation "Central body")
   (gfun :initarg :gfun :accessor gfun :documentation "Gravity function")
   (area :initarg :area :accessor area :documentation "Sail area")
   (mass :initarg :mass :accessor mass :documentation "Spacecraft mass")
   (pfun :initarg :pfun :accessor pfun :documentation "Pointing function")
   (afun :initarg :afun :accessor afun :documentation "Acceleration function")
   (bfun :initarg :bfun :accessor bfun :documentation "Basis function for pointing")
   (nfun :initarg :nfun :accessor nfun :documentation "Function to call getting sail normal vector from sail body frame")
   (iframe :initarg :iframe :accessor iframe :documentation "Inertial frame")
   (rs :initarg :rs :Accessor rs :documentation "Rotor wrt basis frame")
   (tm :initarg :tm :initform 0 :accessor tm :documentation "Stores time variable in state")
   ;; Derived slots
   (rm2 :documentation "Norm-squared of position vector")
   (rm :documentation "Norm of position vector")
   (ru :documentation "Unit position vector")
   (g :documentation "Gravitational acceleration")
   (p :documentation "Solar pressure")
   (bframe :documentation "Basis frame")
   (rb :documentation "Rotor of basis frame")
   (rp :documentation "Pointing rotor")
   (pframe :documentation "Pointing frame")
   (n :documentation "Sail normal vector")
   (a :documentation "External acceleration"))
  (:documentation "Cartesian state class"))

(defmethod print-object ((x cartstate) stream)
  (format stream "#<CARTSTATE :r ~a :v ~a>" (r x) (v x)))

(defstatearithmetic cartstate (r v) :oslots (cb gfun area mass pfun afun bfun nfun iframe rs tm))

(defmacro def-unbound (slot (var class) &body body)
  "Define slot reader function that executes BODY to establish the value of the slot for later reads"
  `(defmethod ,slot ((,var ,class))
     (if (slot-boundp ,var ',slot)
	 (slot-value ,var ',slot)
	 (setf (slot-value ,var ',slot)
	       ,@body))))
