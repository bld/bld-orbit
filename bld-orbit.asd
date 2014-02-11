(asdf:defsystem :bld-orbit
  :author "Ben Diedrich"
  :version "0.0.1"
  :license "MIT"
  :description "Orbital mechanics library employing geometric algebra. Currently working to solve solar sail trajectories."
  :depends-on ("bld-ga" "bld-e2" "bld-e3" "bld-ode" "bld-utils" "bld-gen" "anaphora" "local-time" "geco" "alexandria")
  :serial t
  :components
  ((:file "package")
   (:file "constants")
   (:file "frames")
   (:file "force")
   (:file "state")
   (:file "spacecraft")
   (:file "orbit")
   (:file "body")
   (:file "eom")
   (:file "propagate")
   (:file "export")
   (:file "examples")
   (:file "examples-old")
   (:file "geco-ks")
   ))
