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

(defun cartesian-test-01 ()
  (make-instance
   'cartesian-problem
   :central-body *sun*
   :ref :eclipj2000
   :utc0 *utc-j2000*
   :utcf (+ *utc-j2000* (* 365.25 (convert-unit 'day 'second)))
   :eom #'eom
   :x0 (to-cartesian *utc-j2000* *earth* :observer :sun :ref :eclipj2000)
   :hmin (/ 27000 2)
   :hmax (* 365.25 (convert-unit 'day 'second))
   :tol 1d-11))

(defun spinor-test-01 ()
  (to-spinor-problem (cartesian-test-01)))

(defun ks-test-01 ()
  (to-ks-problem (spinor-test-01)))
