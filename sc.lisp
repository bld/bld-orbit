(in-package :bld-orbit)

(defmethod rtn-frame ((x cartesian-state))
  (with-slots (r v) x
    (let ((h (*o r v)))
      (list
       (unitg r)
       (unitg (*i r h))
       (unitg (dual h))))))

(let* ((ac (* 0.1d0 1d-6)) ; characteristic acceleration (mm/s^2 -> km/s^2)
       (au 1.49597870700d8) ; AU (km)
       (sia (atan (sqrt (/ 2))))) ; sun incidence angle
  (defmethod forcefun-eg (et (x cartesian-state) p)
    (with-slots (r v) x
      (destructuring-bind (rh th nh) (rtn-frame x)
	(let* ((q (rotor (*o rh th) sia))
	       (n (rotateg rh q)))
	  (* 2d0 ; double for reflection
	     ac ; acceleration at 1 AU and 0 sun incidence
	     (expt (/ au (norme r)) 2) ; inverse square variation
	     (expt (scalar (*i rh n)) 2) ; cosine-squared SIA variation
	     n)))))) ; normal vector direction of thrust

(defmethod eom-sail-fixed (et (x cartesian-state) p)
  (with-slots (central-body forcefun) p
    (with-slots (r v) x
      (make-instance
       'cartesian-state
       :r v
       :v (+ (gravity r central-body)
	     (forcefun-eg et x p))))))

