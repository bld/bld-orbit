(in-package :bld-orbit2)

(defparameter *eg-sc*
  (make-instance 
   'sc 
   :cb *sun*
   :nb nil
   :gfun #'gravity))

(defparameter *eg-sail-fixed*
  (make-instance
   'sail
   :cb *sun*
   :gfun #'gravity
   :area 1200
   :mass 50
   :pfun #'fixed
   :afun #'ideal-sail-acc
   :bfun #'orbit-frame
   :iframe *j2000*
   :rs (rotor (bve3 :e1e3 -1) (atan (sqrt (/ 2))))))

(defparameter *eg-sail-table*
  (make-instance
   'sail
   :cb *sun*
   :gfun #'gravity
   :area 1200
   :mass 50
   :pfun #'table
   :afun #'ideal-sail-acc
   :bfun #'orbit-frame
   :iframe *j2000*
   :rs (let ((year (convert-unit 'years 'seconds)))
	 (list (list year (re3 :s 1))
	       (list (* 2 year) (rotor (bve3 :e1e3 -1) (atan (sqrt (/ 2)))))))))
