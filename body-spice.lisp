(in-package :bld-orbit)

(defclass body ()
  ((ephemeris :initarg :ephemeris)
   (name :initarg :name)
   (ref :initarg :ref)
   (abcorr :initarg :abcorr)
   (center :initarg :center)
   (mu :initarg :mu)
   (ls :initarg :ls)))

(defparameter *j2000-utc-seconds*
  (+ (encode-universal-time 55 58 11 1 1 2000 0) 0.816d0))

(defun utcsec-to-ephemeris-time (utcsec)
  (- utcsec *j2000-utc-seconds*))

(defmethod position-velocity ((b body) teph)
  "Position & velocity of body"
  (with-slots (ephemeris name ref abcorr center) b
    (with-kernel ephemeris
      (let ((pvv (spk-ezr name (utcsec-to-ephemeris-time teph) center :ref ref :abcorr abcorr)))
	(make-instance 
	 'cartstate
	 :r (ve3 :e1 (aref pvv 0) :e2 (aref pvv 1) :e3 (aref pvv 2))
	 :v (ve3 :e1 (aref pvv 3) :e2 (aref pvv 4) :e3 (aref pvv 5)))))))

(defparameter *ssb*
  (make-instance 'body
		 :ephemeris "p:/src/ephemeris/de421.bsp"
		 :name :ssb
		 :center :ssb
		 :ref :j2000
		 :abcorr :none
		 :mu (with-kernel "p:/src/ephemeris/gm_de431.tpc"
		       (body-vrd :sun :gm 1))))

(defparameter *sun*
  (make-instance 'body
		 :ephemeris "p:/src/ephemeris/de421.bsp"
		 :name :sun
		 :center :ssb
		 :ref :j2000
		 :abcorr :none
		 :mu (with-kernel "p:/src/ephemeris/gm_de431.tpc"
		       (body-vrd :sun :gm 1))
		 :ls 1361.3477d0))

(defparameter *earth*
  (make-instance 'body
		 :ephemeris "p:/src/ephemeris/de421.bsp"
		 :name :earth
		 :center :ssb
		 :ref :j2000
		 :abcorr :none
		 :mu (with-kernel "p:/src/ephemeris/gm_de431.tpc"
		       (body-vrd :sun :gm 1))))
