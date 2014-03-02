(in-package :bld-orbit)

;;; Kustaanheimo-Stiefel Orbit Element equations of motion

(defclass ksstate ()
  ((alpha :initarg :alpha :documentation "U at s=0")
   (beta :initarg :beta :documentation "dU/ds / w0 at s=0")
   (e :initarg :e :documentation "Keplerian specific orbital energy")
   (tm :initarg :tm :documentation "Time"))
  (:documentation "Kustaanheimo-Stiefel orbital element state"))

(let ((slots '(alpha beta e tm)))
  (defmethod print-object ((x ksstate) s)
    (format s "#<KSSTATE ~{~a~^ ~}>"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value x slot))))))

(defstatearithmetic ksstate (alpha beta e tm))

(defmethod to-spinor ((x ksstate) s sc)
  "Convert KS state to spinor state"
  (with-slots (alpha beta e tm) x
    (with-slots (x0) sc
      (let* ((e0 (energy x0 sc))
	     (hk0 (- e0))
	     (w0 (sqrt (/ hk0 2)))
	     (hk (- e)))
	(values
	 (make-instance
	  'spinorstate
	  :u (+ (* alpha (cos (* w0 s)))
		(* beta (sin (* w0 s))))
	  :duds (* w0
		   (- (* beta (cos (* w0 s)))
		      (* alpha (sin (* w0 s)))))
	  :tm tm)
	 s)))))

(defmethod to-cartesian ((x ksstate) s sc)
  "Convert KS state to cartesian"
  (to-cartesian (to-spinor x s sc) s sc))

(defmethod energy ((x ksstate) sc)
  (slot-value x 'e))

(defmethod time-of (s (x ksstate))
  (slot-value x 'tm))

(defmethod to-ks ((x spinorstate) s sc)
  "Convert spinor state to KS"
  (with-slots (u duds tm) x
    (with-slots (x0) sc
      (let ((e0 (energy x0 sc)))
	(let* ((hk0 (- e0))
	       (w0 (sqrt (/ hk0 2)))
	       (e (energy x sc))
	       (hk (- e))
	       (alpha (- (* u (cos (* w0 s)))
			 (* (/ duds w0) (sin (* w0 s)))))
	       (beta (+ (* u (sin (* w0 s)))
			(* (/ duds w0) (cos (* w0 s))))))
	  (make-instance
	   'ksstate
	   :alpha alpha
	   :beta beta
	   :e e
	   :tm tm))))))

(defmethod to-ks ((x cartstate) tm sc)
  "Convert cartesian state to KS"
  (to-ks (to-spinor x 0 sc) tm sc))

(defmethod to-initial-ks ((x spinorstate) s sc)
  "Generate initial KS state from spinor state, s (assumed 0), and SC"
  (with-slots (u duds tm) x
    (let* ((e0 (energy x sc))
	   (hk0 (- e0))
	   (w0 (sqrt (/ hk0 2)))
	   (alpha0 u)
	   (beta0 (/ duds w0)))
      (values
       (make-instance
	'ksstate
	:alpha alpha0
	:beta beta0
	:e e0
	:tm tm)
       0))))

(defmethod to-initial-ks ((x cartstate) tm sc)
  (to-initial-ks (to-initial-spinor x tm sc) 0 sc))  

(defmethod to-ks (x s sc)
  "Default: try cartesian first"
  (multiple-value-bind (xc tm) (to-cartesian x s sc)
    (to-ks xc tm sc)))
