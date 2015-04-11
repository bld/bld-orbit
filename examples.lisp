(in-package :bld-orbit)

;; Ephemeris kernels (SPK)

(defparameter *spk-dir*
  (asdf:system-relative-pathname :bld-orbit "ephemeris/"))
(defparameter *spk* "DE421AllPlanets.bsp")
(defparameter *spk-path*
  (princ-to-string (merge-pathnames-as-file *spk-dir* *spk*)))

;; Leap second kernels (LSK)

(defparameter *lsk-dir*
  (asdf:system-relative-pathname :bld-orbit "ephemeris-time/"))
(defparameter *lsk* "naif0010.tls")
(defparameter *lsk-path*
  (princ-to-string (merge-pathnames-as-file *lsk-dir* *lsk*)))

;; Bodies

(defparameter *ssb*
  (make-instance
   'body
   :name :ssb
   :center nil
   :ref :eclipj2000
   :mu nil))

(defparameter *sun*
  (make-instance
   'body
   :name :sun
   :center :ssb
   :ref :eclipj2000
   :mu 132712440017.99d0))

(defparameter *earth*
  (make-instance
   'body
   :name :earth
   :center :ssb
   :ref :eclipj2000
   :mu 398600.4415d0))

;;; Test 01
;;; Spacecraft with initial condition of Earth at J2000 with only solar gravity

(defun cartesian-test-01 ()
  (make-instance
   'cartesian-problem
   :central-body *sun*
   :nbodies nil
   :spk (list *spk-path*)
   :lsk *lsk-path*
   :ref :eclipj2000
   :et0 0d0
   :etf (convert-unit '(365.25 days) 'sec)
   :eom #'eom-nbody
   :x0 (with-kernels (*spk-path* *lsk-path*)
	 (to-cartesian 0d0 *earth* :observer :sun :ref :eclipj2000))
   :hmin 1d0
   :hmax 27000d0
   :tol 1d-11))

(defun spinor-test-01 ()
  (to-spinor-problem (cartesian-test-01)))

#|
(defun ks-test-01 ()
  (to-ks-problem (spinor-test-01)))
|#

;;; Test 02
;;; Spacecraft with initial condition at Sun-Earth L1 with solar and Earth gravity

(defun cartesian-test-02 ()
  (make-instance
   'cartesian-problem
   :central-body *sun*
   :nbodies (list *earth*)
   :spk (list *spk-path*)
   :lsk *lsk-path*
   :ref :eclipj2000
   :et0 0d0
   :etf (convert-unit '(365.25 days) 'sec)
   :eom #'eom-nbody
   :x0 (make-instance
	'cartesian-state
	:r (ve3 :e1 -26234827.5935499d0 :e2 143254606.4839523d0 :e3 -605.1121386736631d0)
        :v (ve3 :e1 -29.49719930077679d0 :e2 -5.414763857333822d0 :e3 0.000179932481943812d0))
   :hmin 1d2
   :hmax 27000d0
   :tol 1d-12))
