;;; Model bodies

(in-package :bld-orbit)

(defun rad (deg)
  (* (/ pi 180) deg))

(defun deg (rad)
  (* (/ 180 pi) rad))

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

(defclass planet ()
  ((sma0 :initarg :sma0 :documentation "Semi-major axis (AU)") (dsma :initarg :dsma :documentation "Semi major axis derivative (AU/Cy)")
   (ecc0 :initarg :ecc0 :documentation "Eccentricity") (decc :initarg :decc :documentation "Eccentricity derivative (1/Cy)")
   (inc0 :initarg :inc0 :documentation "Inclination (deg)") (dinc :initarg :dinc :documentation "Inclination derivative (deg/Cy)")
   (mlon0 :initarg :mlon0 :documentation "Mean longitude (deg)") (dmlon :initarg :dmlon :documentation "Mean longitutde derivative (deg/Cy)")
   (lper0 :initarg :lper0 :documentation "Longitude of perihelion (deg)") (dlper :initarg :dlper :documentation "Longitude of perihelion derivative (deg/Cy)")
   (lan0 :initarg :lan0 :documentation "Longitude of ascending node (deg)") (dlan :initarg :dlan :documentation "Longitude of ascending node (deg/Cy)")
   (jd :initarg :jd :documentation "Julian date of ephemeris")
   (mu :initarg :mu :documentation "Gravitational parameter (km^3/s^2")
   (basis :initarg :basis :documentation "Basis of planet")))

(defparameter *earth-barycenter*
  (make-instance
   'planet
   :sma0 1.00000261d0 :dsma 0.00000562d0
   :ecc0 0.01671123d0 :decc -0.00004392d0
   :inc0 -0.00001531d0 :dinc -0.01294668d0
   :mlon0 100.46457166d0 :dmlon 35999.37244981d0
   :lper0 102.93768193d0 :dlper 0.32327364d0
   :lan0 0d0 :dlan 0d0
   :jd 2451545d0
   :mu 398600.440d0 ; km^3/s^2
   :basis *j2000*))

(defun mod-mean-anomaly (m)
  "Modulus mean anomaly to be between -180 and 180 degrees"
  (let ((emr (- m (* 360 (floor m 360)))))
    (cond
      ((< emr -180) (+ emr 360))
      ((> emr 180) (- emr 360))
      (t emr))))

(defun solve-kepler (m es &optional (tol 1d-6))
  "Solve Kepler problem for given mean longitude, eccentricity (as degrees), and tolerance."
  (loop for en = (+ m (* es (sin (rad m)))) then (+ en de)
     for dm = (- m (- en (* es (sin (rad en)))))
     for de = (/ dm (- 1 (* (rad es) (cos (rad en)))))
     until (<= (abs de) tol)
     finally (return en)))

(defun universal-to-julian-date (ut)
  (let ((day-frac (mod (/ ut 60d0 60d0 24d0) 1)))
    (+ (astronomical-julian-date (universal-to-timestamp ut)) day-frac (- 0.5))))

(defgeneric position-vector (p teph &key) (:documentation "Position of object at ephemeris time"))

(defmethod position-vector ((p planet) teph &key (tol 1d-6))
  (with-slots (sma0 dsma ecc0 decc inc0 dinc mlon0 dmlon lper0 dlper lan0 dlan jd basis) p
    (let* ((teph-jd (universal-to-julian-date teph))
	   (t-cy (/ (- teph-jd jd) 36525)) ; time in centuries past Julian date of planet parameters
	   ;; Elements at ephemeris time
	   (sma (+ sma0 (* t-cy dsma))) 
	   (ecc (+ ecc0 (* t-cy decc)))
	   (inc (+ inc0 (* t-cy dinc)))
	   (mlon (+ mlon0 (* t-cy dmlon)))
	   (lper (+ lper0 (* t-cy dlper)))
	   (lan (+ lan0 (* t-cy dlan)))
	   (aop (- lper lan)) ; Argument of periheliion
	   (m (mod-mean-anomaly (- mlon lper))) ; mean anomaly
	   (es (* (/ 180 pi) ecc)) ; eccentricity expressed as degrees
	   (e (solve-kepler m es tol))
	   (orb-coefs (list (* sma (- (cos (rad e)) ecc))
			    (* sma (sqrt (- 1d0 (expt ecc 2))) (sin (rad e)))
			    0))
	   (rtr-lan (rotor (- (dual (third basis))) (rad lan)))
	   (b-lan (new-frame rtr-lan basis))
	   (rtr-inc (rotor (- (dual (first b-lan))) (rad inc)))
	   (b-inc (new-frame rtr-inc b-lan))
	   (rtr-aop (rotor (- (dual (third b-inc))) (rad aop)))
	   (b-aop (new-frame rtr-aop b-inc))
	   (rtr-eqt (rotor (- (dual (first b-aop))) (rad -23.43928)))
	   (b-eqt (new-frame rtr-eqt b-aop))
	   )
      (make-vector orb-coefs b-eqt)
      ;;(make-vector orb-coefs b-aop)
      )))

