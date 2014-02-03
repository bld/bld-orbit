;;; Force related functions

(in-package :bld-orbit2)

;;; Gravity

(defmethod gravity ((x cartstate))
  (let ((r (r x))
	(b (cb (sc x))))
    (- (* (/ (mu b) (norme2 r)) (unitg r)))))

