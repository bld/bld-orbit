(defpackage :bld-orbit
  (:use :cl :bld-ga :bld-e3 :bld-e2 :bld-ode :local-time :cl-spice :unit-formulas :cl-fad :cl-csv)
  (:shadowing-import-from :bld-gen
    + - * / expt
    sin cos tan
    atan asin acos
    sinh cosh tanh
    asinh acosh atanh
    log exp sqrt abs
    min max signum)
  (:export
   ;; Time
   *utc-j2000* utc-to-timestamp
   ;; Stumpff functions
   c0 c1
   ;; Body
   body gravity
   ;; Cartesian
   cartesian-state cartesian-problem to-cartesian eom propagate convert-results
   ;; Examples
   *ephemeris-dir* *planetary-ephemeris* *planetary-ephemeris-path*
   *ssb* *sun* *earth*
   cartesian-test-01))
