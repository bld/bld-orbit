;;; Stumpf Functions

(in-package :bld-orbit)

(defmethod c0 ((x number))
  (cond
    ((zerop x) 1)
    ((> x 0) (cos (sqrt x)))
    ((< x 0) (cosh (sqrt (- x))))))

(defmethod c1 ((x number))
  (cond
    ((zerop x) 1)
    ((> x 0) (/ (sin (sqrt x)) (sqrt x)))
    ((< x 0)(/ (sinh (sqrt (- x))) (sqrt (- x))))))
