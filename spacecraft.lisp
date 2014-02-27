;;; Spacecraft data classes

(in-package :bld-orbit)

(defclass sail ()
  ((eom :initarg :eom :initform #'eom :documentation "Equations of motion")
   (cb :initarg :cb :initform (make-instance 'body :mu 1d0) :documentation "Central body")
   (sun :initarg :sun :documentation "Sun (if it's different from CB)")
   (gfun :initarg :gfun :documentation "Gravity function")
   (accfun :initarg :accfun :initform #'sail-ideal-acc :documentation "Sail acceleration function")
   (basisfun :initarg :basisfun :initform #'rvbasis :documentation "Orbital basis function")
   (pointfun :initarg :pointfun :initform #'sail-pointing-fixed :documentation "Sail pointing function")
   (nfun :initarg :nfun :initform #'first :documentation "Function to call on sail body frame to get normal vector")
   (lightness :initarg :lightness :initform 0)
   (mass :initarg :mass :documentation "Mass of spacecraft (kg)")
   (area :initarg :area :documentation "Sail area (km^2, m^2*1d-6)")
   (optical :initarg :optical)
   (basis :initarg :basis :initform (list (ve2 :e1 1) (ve2 :e2 1)))
   (t0 :initarg :t0 :initform 0d0 :documentation "Initial time")
   (tf :initarg :tf :initform (* 2 pi) :documentation "Final time")
   (x0 :initarg :x0 :initform (make-instance 'cartstate :r (ve2 :e1 1) :v (ve2 :e2 1)))
   (rs :initarg :rs :initform (re2 :s 1d0) :documentation "Sail orientation rotor wrt orbital position frame")
   (rs-table :initarg :rs-table :documentation "Sail orientation rotor wrt orbital position frame")
   (outfile :initarg :outfile :initform "sail-2d-cart-kepler-eg.dat" :documentation "Output filename"))
  (:documentation "Solar sail orbit problem. Default is 2D Kepler with circular orbit, and units of AU and TU."))

(let ((slots '(eom cb accfun pointfun lightness basis t0 tf x0 rs)))
  (defmethod print-object ((sail sail) s)
    (format s "#<SAIL ~{~a~^~%  ~}>"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value sail slot))))))
