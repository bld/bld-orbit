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

(let ((slots '(eom cb mu ls x0 basis)))
  (defmethod print-object ((b body) s)
    (format s "#<BODY ~{~a~^~%~}>"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value b slot))))))

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
  (make-instance
   'body
   :cb *sun*
   :mu 398600.440 ; km^3/s^2
   :ls 0
   :basis *j2000*
   :x0 (let* ((xc (make-instance
		   'cartstate
		   :r (ve3 :e1 3.358497373728686d07 :e2 1.431508445244126d08 :e3 -1.825726859413210d04)
		   :v (ve3 :e1 -2.948987304368306d01 :e2 6.645626180478944d00 :e3 1.289207633035708d-04)))
	      (b (make-instance 'body :cb *sun* :basis *j2000* :x0 xc)))
	 (to-ks xc 0 b))))

  #|
(defun body-position (b tm)
  (with-slots (x0 cb) b
    (if cb
	(with-slots (mu) cb
	  (with-slots (a e) (to-coe x0 tm b)
	    (let* ((n (sqrt (/ mu (expt a 3))))
		   (m ( *
|#

(defclass planet ()
  ((sma :initarg :sma)
   (ecc :initarg :ecc)
   (inc :initarg :inc)
   (mlon :initarg :mlon)
   (lper :initarg :lper)
   (lan :initarg :lan)
   (dsma :initarg :dsma)
   (decc :initarg :decc)
   (dinc :initarg :dinc)
   (dmlon :initarg :dmlon)
   (dlper :initarg :dlper)
   (dlan :initarg :dlan)
   (jd :initarg :jd)
   (mu :initarg :mu)))

(defparameter *earth-bary*
  (make-instance
   'planet
   :sma 1.00000018 :dsma -0.00000003
   :ecc 0.01673163 :decc -0.00003661
   :inc -0.00054346 :dinc -0.01337178
   :mlon 100.46691572 :dmlon 35999.37306329
   :lper 102.93005885 :dlper 0.31795260
   :lan -5.11260389 :dlan -0.24123856
   :jd 2451545d0
   :mu 398600.440)) ; km^3/s^2