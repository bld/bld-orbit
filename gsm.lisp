;;; Generalized Model for Solar Sails

(ql:quickload :cl-json)

(in-package :bld-orbit)

(defclass material ()
  ((a1 :initarg :a1)
   (a2 :initarg :a2)
   (a3 :initarg :a3)
   (r :initarg :r)
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
    (if (and (zerop ef) (zerop eb))
	(* bf (- 1 s) r) ; zero emissivity
	(+ (* bf (- 1 s) r) ; non-zero emissivity
	   (* (- 1 r)
	      (- (* ef bf)
		 (* eb bb))
	      (/ (+ ef eb)))))))

(defmethod a3 ((m material))
  "A3 derived optical property"
  (with-slots (r s) m
    (- 1 (* r s))))

(defmethod initialize-instance :after ((m material) &key)
  (when (not (slot-boundp m 'a1))
    (setf (slot-value m 'a1) (a1 m)))
  (when (not (slot-boundp m 'a2))
    (setf (slot-value m 'a2) (a2 m)))
  (when (not (slot-boundp m 'a3))
    (setf (slot-value m 'a3) (a3 m))))

(defclass gsm-coef ()
  ((m :initarg :m)
   (j1 :initarg :j1)
   (j2 :initarg :j2)
   (j3 :initarg :j3)
   (l :initarg :l)
   (k2 :initarg :k2)
   (k3 :initarg :k3))
  (:documentation "Coefficients for use in the GSM force and moment calculations"))

(defmethod linear-function-2 ((m2 list) (v ve3))
  "Linear function of a list of VE3 vectors and a VE3"
  (make-instance
   've3
   :coef
   (map 'vector #'(lambda (m2v) (scalar (*i m2v v))) m2)))

(defun transpose-ve3-list (l)
  (loop for i below (dimension (first l))
     collect
       (make-instance 've3 :coef (map 'vector #'(lambda (v) (aref (coef v) i)) l))))

(defmethod linear-function-3 (m3 (v ve3))
  "Linear function of a list of lists of VE3 vectors with a vector"
  (let* ((m2 (mapcar #'(lambda (vv) (linear-function-2 vv v)) m3))
	 (m2t (transpose-ve3-list m2)))
    (linear-function-2 m2t v)))

(defmethod gsm-force (p (sunvb ve3) (c gsm-coef))
  (with-slots (m j1 j2 j3 l k2 k3) c
    (with-slots (a1 a2 a3) m
      (* p (- (* a2 (linear-function-2 j2 sunvb))
	      (* a1 (linear-function-3 j3 sunvb))
	      (* a3 (scalar (*i j1 sunvb)) sunvb))))))

(defmethod gsm-moment (p (sunvb ve3) (c gsm-coef))
  (with-slots (m j1 j2 j3 l k2 k3) c
    (with-slots (a1 a2 a3) m
      (* p (- (* a2 (linear-function-2 k2 sunvb))
	      (* a1 (linear-function-3 k3 sunvb))
	      (* a3 (*x (linear-function-2 (transpose-ve3-list l) sunvb) sunvb)))))))

(defmethod gsm-force-moment ((p number) (solarvb ve3) (c gsm-coef))
  (values
   (gsm-force p solarvb c)
   (gsm-moment p solarvb c)))

(defmethod tensor2-to-ve3 ((a array))
  "Convert 2D tensor's columns to list of VE3 vectors"
  (let ((d (array-dimensions a)))
    (loop for i below (second d)
       collect
	 (let ((v (ve3)))
	   (dotimes (j (first d))
	     (setf (aref (coef v) j) (aref a j i)))
	   v))))

#|
(defmethod tensor3-to-ve3 ((a array))
  "Convert columns of 3D tensor to list (of lists) of VE3 vectors"
  (destructuring-bind (d1 d2 d3) (array-dimensions a)
    (loop for i below d1
       collect
	 (loop for j below d3
	    collect
	      (let ((v (ve3)))
		(dotimes (k d2)
		  (setf (aref (coef v) k) (aref a i k j)))
		v)))))
|#

(defun parse-tensor (tdef)
  (let ((size (rest (assoc :--*array-size-- tdef)))
	(type (intern (string-upcase (rest (assoc :--*array-type-- tdef)))))
	(data (rest (assoc :--*array-data-- tdef))))
    (let ((a (make-array size :element-type type)))
      (loop for ai in data
	 for i below (apply #'* size)
	 do (setf (row-major-aref a i) ai))
      a)))

(defun decode-coef-json (f)
  (with-open-file (s f)
    (let ((parsed (json:decode-json s)))
      (list
       (make-array 3 :initial-contents (mapcar #'first (rest (assoc :j-1 parsed))))
       (make-array '(3 3) :initial-contents (rest (assoc :j-2 parsed)))
       (parse-tensor (rest (assoc :j-3 parsed)))
       (make-array '(3 3) :initial-contents (rest (assoc :l parsed)))
       (make-array '(3 3) :initial-contents (rest (assoc :k-2 parsed)))
       (parse-tensor (rest (assoc :k-3 parsed)))))))

(defun load-coef-json (f m)
  (destructuring-bind (j1 j2 j3 l k2 k3)
      (let ((*read-default-float-format* 'double-float))
	(decode-coef-json f))
    (make-instance
     'gsm-coef
     :m m
     :j1 (make-instance 've3 :coef j1)
     :j2 (tensor2-to-ve3 (bld-linalg:transpose j2))
     :j3 (tensor3-to-ve3 j3)
     :l (tensor2-to-ve3 (bld-linalg:transpose l))
     :k2 (tensor2-to-ve3 (bld-linalg:transpose k2))
     :k3 (tensor3-to-ve3 k3))))
