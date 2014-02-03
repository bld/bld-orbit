(asdf:defsystem :bld-orbit2
  :author "Ben Diedrich"
  :version "0.0.1"
  :license "MIT"
  :description "Orbital mechanics library employing geometric algebra. Currently working to solve solar sail trajectories."
  :depends-on ("bld-ga" "bld-e2" "bld-e3" "bld-ode" "bld-gen")
  :serial t
  :components
  ((:file "package2")
   (:file "constants")
   (:file "frames")
   (:file "state")
   (:file "spacecraft")
   (:file "body2")
   (:file "force")
   (:file "eom2")
   (:file "eg")))
