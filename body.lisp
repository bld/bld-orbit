;;; Model bodies

(in-package :bld-orbit)

(defun rad (deg)
  (* (/ pi 180) deg))

(defun deg (rad)
  (* (/ 180 pi) rad))

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

(defclass body ()
  ((cb :initarg :cb :documentation "Central body")
   (mu :initarg :mu :documentation "Gravitational parameter")
   (ls :initarg :ls :documentation "Solar luminosity at 1 AU in W/m^2")
   (iframe :initarg :iframe :documentation "Inertial frame of body's orbit")
   (jd :initarg :jd :documentation "Julian date of ephemeris")
   (utc :initarg :utc :documentation "UTC of ephemeris")
   (ec :initarg :ec :documentation "Eccentricity")
   (qr :initarg :qr :documentation "Periapsis distance")
   (in :initarg :in :documentation "Inclination")
   (om :initarg :om :documentation "Longitude of ascending node")
   (w :initarg :w :documentation "Argument of perifocus")
   (tp :initarg :tp :documentation "Time of periapsis passage")
   (n :initarg :n :documentation "Mean motion")
   (ma :initarg :ma :documentation "Mean anomaly")
   (ta :initarg :ta :documentation "True anomaly")
   (a :initarg :a :documentation "Semi-major axis")
   (ad :initarg :ad :documentation "Apoapsis distance")
   (pr :initarg :pr :documentation "Orbital period")
   (xcoe :documentation "COE initial state")
   (xcart :documentation "Cartesian initial state")
   (xspin :documentation "Spinor initial state")
   (xks :documentation "KS initial state")
   (x0 :initarg :x0))
  (:documentation "Body from ephemeris provided by JPL HORIZONS"))

(defmethod initialize-instance :after ((b body) &key)
  (when (every #'(lambda (slot) (slot-boundp b slot))
	       '(a ec in om w jd))
    (with-slots (cb mu jd utc ec qr in om w tp n ma ta a ad pr xcoe xcart xspin xks x0) b
      (setf xcoe (make-instance
		  'coestate
		  :a a
		  :e ec
		  :i (rad in)
		  :lan (rad om)
		  :aop (rad w)
		  :tm jd))
      (setf xcart (to-cartesian xcoe (rad ta) b))
      (setf xspin (to-spinor xcart utc b))
      (setf x0 xspin)
      (setf xks (to-ks xspin 0 b)))))

(defmethod position-velocity ((b body) teph)
  "Position and velocity of body given UTC time"
  (if (slot-boundp b 'cb)
      (with-slots (cb mu iframe jd ec qr in om w tp n ma ta a ad pr xks) b
	(with-slots (alpha beta e tm) xks
	  (let* ((teph-jd (universal-to-julian-date teph))
		 (t-jd (- teph-jd jd))
		 (m (+ ma (* t-jd n 60d0 60d0 24d0)))
		 (ea (solve-kepler m (deg ec)))
		 (w0 (sqrt (/ (- e) 2)))
		 (s (/ (rad ea) 2 w0)))
	    (to-cartesian xks s b))))
      ;; Return 0 position & velocity if there's no central body
      (make-instance 'cartstate :r (ve3) :v (ve3))))

(defmethod eom (tm (x cartstate) (b body))
  (with-slots (v (tm-x tm)) x
    (setf tm-x tm)
    (make-instance
     'cartstate
     :r v
     :v (gravity tm x b)
     :tm tm
     :sc b)))

(defparameter *ssb* (make-instance 'body :mu 1.32712440018d11))

(defparameter *sun* 
  (make-instance 
   'body
   :cb *ssb*
   :mu 1.32712440018d11
   :ls 1361.3477d0
   :iframe *j2000*
   :jd 2456636.500000000d0
   :utc (coerce (encode-universal-time 0 0 0 10 12 2013 0) 'double-float)
   :EC 8.546770937522187d-01 :QR 2.906237974575245d+04 :IN 2.400537289551978d+00
   :OM 5.010638625423464d+01 :W 6.029921568009870d+01 :Tp 2456456.332998391241d0
   :N 1.145648031345176d-05 :MA 1.783364866902255d+02 :TA 1.797489228422292d+02
   :A 1.999848509511635d+05 :AD 3.709073221565746d+05 :PR 3.142326352861636d+07))

(defparameter *earth*
  ;; Units: km-s-deg
  (make-instance
   'body
   :cb *ssb*
   :mu 398600.44d0
   :ls nil
   :iframe *j2000*
   :jd 2457388.500000000d0
   :utc (coerce (encode-universal-time 0 0 0 1 1 2016 0) 'double-float)
   :ec 1.621551291914839d-02
   :qr 1.471327388424463d+08
   :in 1.148321892914738d-02
   :om 1.970427571324530d+02
   :w 2.776747994466758d+02
   :tp 2457403.399663229473d0
   :n 1.141967260079594d-05
   :ma 3.452991025586977d+02
   :ta 3.448181346808943d+02
   :a 1.495578968509943d+08
   :ad 1.519830548595422d+08
   :pr 3.152454650713090d+07))

(defparameter *mars*
  (make-instance
   'body
   :cb *sun*
   :mu 42828.3d0
   :ls nil
   :iframe *j2000*
   :jd 2456636.500000000d0
   :utc (coerce (encode-universal-time 0 0 0 1 1 2016 0) 'double-float)
   :ec 9.656324402943112d-02
   :qr 2.051027822461357d+08
   :in 1.847743229702892d+00
   :om 4.950957526635985d+01
   :w 2.866119228258302d+02
   :tp 2456319.160087579396d0
   :n 6.106027092759354d-06
   :ma 1.674160792865442d+02
   :ta 1.695714105727685d+02
   :a 2.270250583570647d+08
   :ad 2.489473344679938d+08
   :pr 5.895813997073400d+07))

(defparameter *luna*
  (make-instance
   'body
   :cb *earth*
   :mu 4902.798d0
   :ls nil
   :iframe *j2000*
   :jd 2457388.500000000d0
   :utc (coerce (encode-universal-time 0 0 0 1 1 2016 0) 'double-float)
   :EC 6.431655045408705d-02 :QR 3.512367974406679d+05 :IN 5.071049688763354d+00
   :OM 1.746561927360046d+02 :W 2.015670143350679d+02 :Tp  2457376.756737075746d0
   :N 1.553730434340702d-04 :MA 1.576442736350183d+02 :TA 1.602521673187281d+02
   :A 3.753799403111416d+05 :AD 3.995230831816153d+05 :PR 2.317004237306838d+06))

(defun eom-body-ks (s tm b)
  "KS equation of motion for propagating a body state"
  (with-slots (u duds tm) (to-spinor (slot-value b 'xks) s b)
    (norme2 u)))
