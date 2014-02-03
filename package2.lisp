(defpackage :bld-orbit2
  (:use :cl :bld-ga :bld-e2 :bld-e3 :bld-ode)
  (:shadowing-import-from :bld-gen
    + - * / expt
    sin cos tan
    atan asin acos
    sinh cosh tanh
    asinh acosh atanh
    log exp sqrt abs
    min max signum))
