(ql:quickload :bld-gen)
(ql:quickload :fiveam)
(ql:quickload :bld-orbit)
(ql:quickload :unit-formulas)
(ql:quickload :local-time)

(defpackage :bld-orbit-test
  (:use :cl :bld-orbit :fiveam :unit-formulas :local-time)
  (:shadowing-import-from :bld-gen
			  + - * / expt
			  sin cos tan
			  atan asin acos
			  sinh cosh tanh 
			  asinh acosh atanh
			  log exp sqrt abs
			  min max signum))

(in-package :bld-orbit-test)

(def-suite :bld-orbit)

(in-suite :bld-orbit)

(test utc-to-epht
  (is (= (utc-to-epht *utc-j2000*) 0))
  (is (= (utc-to-epht (+ *utc-j2000* 1000)) 1000)))
