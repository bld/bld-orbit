(in-package :bld-orbit)

(defparameter *J2000* (list (ve3 :e1 1d0)
			    (ve3 :e2 1d0)
			    (ve3 :e3 1d0)) "J2000 inertial frame")

;;; Frame functions 

(defun rv-frame (tm x sc)
  "Form a reference frame from the position and velocity vectors to
CB. First is the unit position vector. Second is the complement of the
first in the orbit plane. Third (if it exists) is the orbit plane
normal vector."
  (with-slots (r v) (to-cartesian x tm sc)
    (with-derived (ru) tm x sc
      (let* ((h (unitg (*o r v)))
	     (x ru)
	     (y (unitg (*i r h))))
	(cond
	  ((= (dimension r) (dimension v) 2) (list x y))
	  ((= (dimension r) (dimension v) 3)
	   (list x y (dual h)))
	  (t (error "R and V must be 2 or 3 dimensional")))))))

(defun orbit-frame (tm x sc)
  "Produce orbit frame from vector to sun"
  (with-derived (rsc-sun) tm x sc
    (with-slots (iframe) sc
      (let* ((o3 (unitg rsc-sun))
	     (o1 (unitg (*i (third iframe) (dual o3))))
	     (o2 (*i o3 (dual o1))))
	(list o1 o2 o3)))))

