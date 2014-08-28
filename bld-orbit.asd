(asdf:defsystem :bld-orbit
  :author "Ben Diedrich"
  :version "0.0.1"
  :license "MIT"
  :description "Orbital mechanics library employing geometric algebra. Currently working to solve solar sail trajectories."
  :depends-on ("bld-ga" "bld-e2" "bld-e3" "bld-ode" "bld-utils" "bld-gen" "anaphora" "local-time" "alexandria" "cl-spice")
  :serial t
  :components
  ((:file "package")
   (:file "orbit")
   (:file "constants")
   (:file "spacecraft")
   (:file "state")
   (:file "cartesian")
   (:file "spinor")
   (:file "ks")
   (:file "coe")
   (:file "mee")
   (:file "derived")
   (:file "frames")
   (:file "point")
   (:file "force")
   ;;(:file "body")
   (:file "body-spice")
   (:file "propagate")
   (:file "export")
   (:file "examples")
   ))
