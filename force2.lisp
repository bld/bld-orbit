;;; Force related functions

(in-package :bld-orbit2)

;;; Gravity

(defmethod gravity ((x cartstate))
  "Force of gravity"
  (- (* (/ (mu (cb x))
	   (rm2 x))
	(ru x))))

;;; Pointing functions

(defun fixed (x)
  "Return spacecraft pointing rotor in inertial space for a fixed
attitude specified by RS slot"
  (*g (rb x) (rs x))) ; rotor of sail attitude in inertial space

;;; Solar radiation pressure functions

(defun solar-pressure (x)
  "Solar radiation pressure experienced by state X"
  (/ (* (ls (sun x))
	*au2*)
     (rm2 x)
     *c*))

;;; Force functions

(defun ideal-sail-acc (x)
  "Force on an ideal solar sail of given state"
  (/ (* 2 (p x) (area x) (expt (scalar (*i (ru x) (n x))) 2) (n x))
     (mass x)))
