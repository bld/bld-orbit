;;; Generalized Model for Solar Sails

(in-package :bld-orbit)

(defclass material ()
  ((r :initarg :r)
   (s :initarg :s)
   (bf :initarg :bf)
   (ef :initarg :ef)
   (bb :initarg :bb)
   (eb :initarg :eb))
  (:documentation "Sail material optical properties"))

(defmethod a1 ((m material))
  "A1 derived optical property"
  (with-slots (r s) m
    (* 2 r s)))

(defmethod a2 ((m material))
  "A2 derived optical property"
  (with-slots (r s bf ef bb eb) m
    (+ (* bf (- 1 s) r)
       (* (- 1 r)
	  (- (* ef bf)
	     (* eb bb))
	  (/ (+ ef eb))))))

(defmethod a3 ((m material))
  "A3 derived optical property"
  (with-slots (r s) m
    (- 1 (* r s))))

(defclass gsm-coef ()
  ((m :initarg :m)
   a1 a2 a3
   (j1 :initarg :j1)
   (j2 :initarg :j2)
   (j3 :initarg :j3)
   (l :initarg :l)
   (k2 :initarg :k2)
   (k3 :initarg :k3))
  (:documentation "Coefficients for use in the GSM force and moment calculations"))

(defmethod linear-function-2 (m2 (v ve3))
  "Linear function of a list of VE3 vectors and a VE3"
  (ve3 :e1 (scalar (*i (first m2) v))
       :e2 (scalar (*i (second m2) v))
       :e3 (scalar (*i (third m2) v))))

(defun transpose-ve3-list (l)
  (list
   (make-instance 've3 :coef (map 'vector #'(lambda (vi) (gref vi :e1)) l))
   (make-instance 've3 :coef (map 'vector #'(lambda (vi) (gref vi :e2)) l))
   (make-instance 've3 :coef (map 'vector #'(lambda (vi) (gref vi :e3)) l))))

(defmethod linear-function-3 (m3 (v ve3))
  "Linear function of a list of lists of VE3 vectors with a vector"
  (let* ((m2 (mapcar #'(lambda (vv) (linear-function-2 vv v)) m3))
	 (m2t (transpose-ve3-list m2)))
    (linear-function-2 m2t v)))

(defmethod initialize-instance :after ((c gsm-coef) &key)
  (with-slots (m a1 a2 a3) c
    (setf
     ;; Calculate derived material properties
     a1 (a1 m)
     a2 (a2 m)
     a3 (a3 m))))

(defmethod gsm-force ((p number) (solarvb ve3) (c gsm-coef))
  (let ((sunvb (- solarvb)))
    (with-slots (a1 a2 a3 j1 j2 j3 l k2 k3) c
      (* p (- (* a2 (linear-function-2 j2 sunvb))
	      (* a1 (linear-function-3 j3 sunvb))
	      (* a3 (scalar (*i j1 sunvb)) sunvb))))))

(defmethod gsm-moment ((p number) (solarvb ve3) (c gsm-coef))
  (let ((sunvb (- solarvb)))
    (with-slots (a1 a2 a3 j1 j2 j3 l k2 k3) c
      (* p (- (* a2 (linear-function-2 k2 sunvb))
	      (* a1 (linear-function-3 k3 sunvb))
	      (* a3 (*x (linear-function-2 l sunvb) sunvb)))))))
		       
(defmethod gsm-force-moment ((p number) (solarvb ve3) (c gsm-coef))
  (let ((sunvb (- solarvb)))
    (with-slots (a1 a2 a3 j1 j2 j3 l k2 k3) c
      (let ((f (* p (- (* a2 (linear-function-2 j2 sunvb))
		       (* a1 (linear-function-3 j3 sunvb))
		       (* a3 (scalar (*i j1 sunvb)) sunvb))))
	    (m (* p (- (* a2 (linear-function-2 k2 sunvb))
		       (* a1 (linear-function-3 k3 sunvb))
		       (* a3 (*x (linear-function-2 l sunvb) sunvb))))))
	(values f m)))))

(defmethod vv-to-ve3 ((vv vector))
  (map 'list #'(lambda (v) (make-instance 've3 :coef v)) vv))

(defmethod vvv-to-ve3 ((vvv vector))
  (map 'list #'(lambda (vv) (vv-to-ve3 vv)) vvv))
