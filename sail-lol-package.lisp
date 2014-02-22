(defpackage :sail-lol
  (:use :cl :lol :bld-ga :bld-e2 :bld-e3 :bld-ode :bld-utils :anaphora)
  (:import-from :alexandria make-keyword)
  (:shadowing-import-from
   :bld-gen
   + - * / expt
   sin cos tan
   atan asin acos
   sinh cosh tanh
   asinh acosh atanh
   log exp sqrt abs
   min max signum))
