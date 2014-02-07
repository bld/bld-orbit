(asdf:defsystem :bld-orbit2
  :author "Ben Diedrich"
  :version "0.0.1"
  :license "MIT"
  :description "Orbital mechanics library employing geometric algebra. Currently working to solve solar sail trajectories."
  :depends-on ("bld-ga" "bld-e2" "bld-e3" "bld-ode" "bld-gen" "unit-formulas" "anaphora")
  :serial t
  :components
  ((:file "package2")
   (:file "constants2")
   (:file "frames2")
   (:file "state2")
   (:file "spacecraft2")
   (:file "body2")
   (:file "force2")
   (:file "eom2")
   (:file "export2")
   (:file "eg2")))
