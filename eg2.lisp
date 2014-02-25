(in-package :bld-orbit2)

(defparameter *eg-sail-fixed-x0*
  (make-instance
   'cartstate
   :r (ve3 :e1 *au*)
   :v (ve3 :e2 (sqrt (/ (mu *sun*) *au*)))
   :cb *sun*
   :gfun #'gravity
   :area 1200
   :mass 50
   :pfun #'fixed
   :afun #'ideal-sail-acc
   :bfun #'rv-frame
   :nfun #'rv-normal
   :iframe *j2000*
   :tm 0
   :rs (rotor (bve3 :e1e2 1) (atan (sqrt (/ 2))))
   ;;:rs (re3 :s 1)
   ))

(defcartesian eg-sail-fixed
    :cb *sun*
    :gfun #'gravity
    :area 1200
    :mass 50
    :pfun #'fixed
    :afun #'ideal-sail-acc
    :bfun #'rv-frame
    :nfun #'rv-normal
    :iframe *j2000*
    :rs (rotor (bve3 :e1e2 1) (atan (sqrt (/ 2)))))

(defparameter *eg-sail-fixed-x0-2* (make-instance 'eg-sail-fixed :tm 0 :r (ve3 :e1 *au*) :v (ve3 :e2 (sqrt (/ (mu *sun*) *au*)))))
x