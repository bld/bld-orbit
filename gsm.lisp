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
  (with-slots (r s) m
    (* 2 r s)))

(defmethod a2 ((m material))
  (with-slots (r s bf ef bb eb) m
    (+ (* bf (- 1 s) r)
       (* (- 1 r)
	  (- (* ef bf)
	     (* eb bb))
	  (/ (+ ef eb))))))

(defmethod a3 ((m material))
  (with-slots (r s) m
    (- 1 (* r s))))

(defclass gsm-coef ()
  ((m :initarg :m)
   a1 a2 a3
   (j1 :initarg :j1)
   j1_v
   (j2 :initarg :j2)
   j2_v
   (j3 :initarg :j3)
   j3_m
   (l :initarg :l)
   (k2 :initarg :k2)
   (k3 :initarg :k3)
   k3_m)
  (:documentation "Coefficients for use in the GSM force and moment calculations"))


(defmethod linear-function-2 (m2 (v ve3))
  "Linear function of a list of VE3 vectors and a VE3"
  (ve3 :e1 (scalar (*i (first m2) v))
       :e2 (scalar (*i (second m2) v))
       :e3 (scalar (*i (third m2) v))))

(defmethod linear-function-3 (m3 (v ve3))
  "Linear function of a list of lists of VE3 vectors with a vector"
  (let ((m2 (loop for m2i in m3
	       collect (linear-function-2 m2i v))))
    (linear-function-2 m2 v)))

(defmethod initialize-instance :after ((c gsm-coef) &key)
  (with-slots (m a1 a2 a3) c
    (setf
     ;; Calculate derived material properties
     a1 (a1 m)
     a2 (a2 m)
     a3 (a3 m))))

(defmethod skew ((vskew ve3))
  "Form a skew list of VE3 vectors of given VE3 vector"
  (list
   (ve3 :e1 0d0 :e2 (gref vskew :e3) :e3 (- (gref vskew :e2)))
   (ve3 :e1 (- (gref vskew :e3)) :e2 0d0 :e3 (gref vskew :e1))
   (ve3 :e1 (gref vskew :e2) :e2 (- (gref vskew :e1)) :e3 0d0)))

(defmethod gsm-force-moment ((p number) (solarvb ve3) (c gsm-coef))
  (let ((sunvb (- solarvb)))
    (with-slots (a1 a2 a3 j1 j2 j3 l k2 k3) c
      (let ((f (* p (- (* a2 (linear-function-2 j2 sunvb))
		       (* a1 (linear-function-3 j3 sunvb))
		       (* a3 (scalar (*i j1 sunvb)) sunvb))))
	    #+null(m (* p (- (* a2 (linear-function-2 k2 sunvb))
		       (* a1 (linear-function-3 k3 sunvb))
		       (* a3 (linear-function-2 (skew sunvb) (linear-function-2 l sunvb)))))))
	(values f)))))
