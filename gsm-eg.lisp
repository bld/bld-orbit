;;; Examples of the Generalized Sail Model

(in-package :bld-orbit)

(defparameter *material-jplhr*
  (make-instance
   'material
   :r 0.88d0
   :s 0.94d0
   :bf 0.79d0
   :ef 0.05d0
   :bb 0.55d0
   :eb 0.55d0)
  "JPL Halley Rendezvous sail material optical properties")

(defparameter *material-ideal*
  (make-instance
   'material
   :r 1d0
   :s 1d0
   :ef 0d0
   :eb 0d0
   :bf 2/3
   :bb 2/3)
  "Ideal sail material")

(defun coef-flat (&key (m *material-ideal*) (l 1))
  (let ((a (expt l 2))
	(l3-half (/ (expt l 3) 2)))
    (make-instance
     'gsm-coef
     :m m
     :j1 (ve3 :e1 0 :e2 0 :e3 a)
     :j2 (tensor2-to-ve3
	  #2a((0 0 0)
	      (0 0 0)
	      (0 0 a)))
     :j3 (tensor3-to-ve3
	  #3a(((0 0 0)
	       (0 0 0)
	       (0 0 0))
	      ((0 0 0)
	       (0 0 0)
	       (0 0 0))
	      ((0 0 0)
	       (0 0 0)
	       (0 0 a))))
     :l (tensor2-to-ve3
	 #2a((0 0 0)
	     (0 0 0)
	     (0 0 0)))
     :k2 (tensor2-to-ve3
	  #2a((0 0 l3-half)
	      (0 0 (- l3-half))
	      (0 0 0)))
     :k3 (tensor3-to-ve3
	  #3a(((0 0 0)
	       (0 0 0)
	       (0 0 0))
	      ((0 0 0)
	       (0 0 0)
	       (0 0 0))
	      ((0 0 l3-half)
	       (0 0 (- l3-half))
	       (0 0 0)))))))

(defparameter *coef-flat-ideal* (coef-flat))

