(asdf:defsystem :sail-lol
  :depends-on ("lol" "bld-ga" "bld-e2" "bld-e3" "bld-gen" "bld-gen" "bld-ode" "alexandria" "bld-utils" "anaphora")
  :serial t
  :components
  ((:file "sail-lol-package")
   (:file "sail-lol-macros")
   (:file "sail-lol")))
