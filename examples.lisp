(in-package :bld-orbit)

(defparameter *ephemeris-dir*
  (asdf:system-relative-pathname :bld-orbit "ephemeris/"))

(defparameter *planetary-ephemeris* "DE421AllPlanets.bsp")

(defparameter *planetary-ephemeris-path*
  (princ-to-string (merge-pathnames-as-file *ephemeris-dir* *planetary-ephemeris*)))

(defparameter *ssb*
  (make-instance
   'body
   :name :ssb
   :center nil
   :ref :eclipj2000
   :ephemeris *planetary-ephemeris-path*
   :mu nil))

(defparameter *sun*
  (make-instance
   'body
   :name :sun
   :center :ssb
   :ref :eclipj2000
   :ephemeris *planetary-ephemeris-path*
   :mu 132712440017.99d0))

(defparameter *earth*
  (make-instance
   'body
   :name :earth
   :center :ssb
   :ref :eclipj2000
   :ephemeris *planetary-ephemeris-path*
   :mu 398600.4415d0))

;;; Test 01
;;; Spacecraft with initial condition of Earth at J2000 with only solar gravity

(defun cartesian-test-01 ()
  (make-instance
   'cartesian-problem
   :central-body *sun*
   :ref :eclipj2000
   :utc0 *utc-j2000*
   :utcf (+ *utc-j2000* (convert-unit '(365.25 day) 'second))
   :eom #'eom-kepler
   :x0 (to-cartesian *utc-j2000* *earth* :observer :sun :ref :eclipj2000)
   :hmin (/ 27000 2)
   :hmax 27000d0
   :tol 1d-11))

(defun spinor-test-01 ()
  (to-spinor-problem (cartesian-test-01)))

(defun ks-test-01 ()
  (to-ks-problem (spinor-test-01)))

;;; Test 02
;;; Spacecraft with initial condition at Sun-Earth L1 with solar and Earth gravity

(defun cartesian-test-02 ()
  (make-instance
   'cartesian-problem
   :central-body *sun*
   :nbodies (list *earth*)
   :ref :eclipj2000
   :utc0 *utc-j2000*
   :utcf (+ *utc-j2000* (convert-unit '(365.25 days) 'sec))
   :eom #'eom-nbody
   :x0 (make-instance
	'cartesian-state
	:r (ve3 :e1 -26234828.18431398d0
		:e2 143254606.368162d0
		:e3 -606.1106522902846d0)
	:v (ve3 :e1 -29.49719927764286
		:e2 -5.414763974791472
		:e3 0.0001800982443111998))
   :hmin 1d0
   :hmax 27000d0
   :tol 1d-13))
