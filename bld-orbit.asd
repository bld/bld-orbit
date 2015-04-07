(asdf:defsystem :bld-orbit
  :name "bld-orbit"
  :author "Ben Diedrich <bldiedrich@gmail.com>"
  :license "MIT"
  :description "Numerical ordinary differential equation solvers"
  :depends-on ("bld-gen" "bld-ga" "bld-e2" "bld-e3" "bld-ode" "local-time" "cl-spice" "unit-formulas" "cl-fad" "cl-csv")
  :serial t
  :components
  ((:file "package")
   (:file "time")
   (:file "body")
   (:file "cartesian")
   (:file "spinor")
   (:file "ks")
   (:file "examples")))
