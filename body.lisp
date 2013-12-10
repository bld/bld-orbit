;;; Model bodies

(in-package :bld-orbit)

(defclass bodystate ()
  ((tm :initarg :tm :initform 0 :documentation "Time")))

(defmethod print-object ((x bodystate) s)
  (format s "#<BODYSTATE :TM ~a>" (slot-value x 'tm)))

(defstatearithmetic bodystate (tm))

(defclass body ()
  ((eom :initarg :eom :initform #'eom :documentation "Equations of motion")
   (cb :initarg :cb :documentation "Central body")
   (mu :initarg :mu :initform 1 :documentation "Gravitational parameter of body")
   (ls :initarg :ls :initform 0 :documentation "Luminosity of body")
   (x0 :initarg :x0 :initform
       (make-instance 'ksstate :tm 0 :e -0.5 :alpha (re2 :s 1) :beta (re2 :e1e2 -1))
       :documentation "Ephemeris of body's orbit")
   (basis :initarg :basis :initform (list (ve2 :e1 1) (ve2 :e2 1)) :documentation "Basis in which multivectors defined")))

(defmethod eom (s (x bodystate) &optional body)
  (with-slots (x0 muc) body
    (with-slots (u duds tm) (to-spinor x0 s body)
      (make-instance
       'bodystate
       :tm (norme2 u)))))

(defparameter *sun*
  (make-instance
   'body
   :cb nil
   :mu 1.32712440018d11 ; km^3/s^2
   :ls 1367.6 ; W/m^2
   :x0 (make-instance 'ksstate :alpha (re3) :beta (re3) :tm 0 :e 0)
   :basis *j2000*))

(defparameter *earth*
  (let ((b (make-instance
	    'body
	    :cb *sun*
	    :mu 398600.440 ; km^3/s^2
	    :ls 0
	    :basis *j2000*)))
    (setf (slot-value b 'x0)
	  (to-ks
	   (make-instance
	    'cartstate
	    :r (ve3 :e1 3.358497373728686E+07 :e2 1.431508445244126E+08 :e3 -1.825726859413210E+04)
	    :v (ve3 :e1 -2.948987304368306E+01 :e2 6.645626180478944E+00 :e3 1.289207633035708E-04))
	   0 b))
    b))

(defun body-position (tm b)
  (with-slots (x0) b
    (with-slots ((alpha0 alpha) (beta0 beta) (e0 e) (tm0 tm)) x0
      