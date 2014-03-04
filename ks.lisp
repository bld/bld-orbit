(in-package :bld-orbit)

;;; Kustaanheimo-Stiefel Orbit Element equations of motion

(defclass ksstate ()
  ((alpha :initarg :alpha :documentation "U at s=0")
   (beta :initarg :beta :documentation "dU/ds / w0 at s=0")
   (e :initarg :e :documentation "Keplerian specific orbital energy")
   (tm :initarg :tm :documentation "Time")
   (derived :accessor derived :initform (make-hash-table) :documentation "Derived values")
   (constants :initarg :constants :accessor constants :documentation "Constants used by EOM"))
  (:documentation "Kustaanheimo-Stiefel orbital element state"))

(defun ks-constants (energy0)
  "Constants for use in KS equations calculated from energy of initial state"
  (make-hash*
   :e0 energy0
   :hk0 (- e0)
   :w0 (sqrt (/ hk0 2))))

(let ((slots '(alpha beta e tm)))
  (defmethod print-object ((x ksstate) s)
    (format s "#<KSSTATE ~{~a~^ ~}>"
	    (loop for slot in slots
	       collect (format nil ":~a ~a" slot (slot-value x slot))))))

(defstatearithmetic ksstate (alpha beta e tm) :oslots (constants))

;;; KS state conversion routines

(defmethod to-spinor ((x ksstate) s sc)
  "Convert KS state to spinor state"
  (with-slots (alpha beta e tm constants) x
    (lethash (e0 hk0 w0) constants
      (let ((hk (- e)))
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
      (let ((constants (ks-constants (energy x0 sc))))
	(lethash (e0 hk0 w0) constants
	  (let* ((e (energy x sc))
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
	     :tm tm
	     :constants constants)))))))

(defmethod to-ks ((x cartstate) tm sc)
  "Convert cartesian state to KS"
  (to-ks (to-spinor x 0 sc) tm sc))

(defmethod to-initial-ks ((x spinorstate) s sc)
  "Generate initial KS state from spinor state, s (assumed 0), and SC"
  (with-slots (u duds tm) x
    (let ((constants (ks-constants (energy x sc))))
      (lethash (e0 hk0 w0) constants
	(let ((alpha0 u)
	      (beta0 (/ duds w0)))
	  (values
	   (make-instance
	    'ksstate
	    :alpha alpha0
	    :beta beta0
	    :e e0
	    :tm tm
	    :constants constants)
	   0))))))

(defmethod to-initial-ks ((x cartstate) tm sc)
  (to-initial-ks (to-initial-spinor x tm sc) 0 sc))

(defmethod to-ks (x s sc)
  "Default: try cartesian first"
  (multiple-value-bind (xc tm) (to-cartesian x s sc)
    (to-ks xc tm sc)))

;;; Equation of motion

(defderived xspin (s (x ksstate) sc)
    "Spinor state"
  (to-spinor x s sc))

(defmethod eom (s (x ksstate) sc)
  (with-slots (alpha beta e tm constants) x
    (with-derived (xspin xcart rm a) s x sc
      (lethash (e0 hk0 w0) constants
	(let* ((r (r xcart))
	       (v (v xcart))
	       (u (slot-value xspin 'u))
	       (hk (- e))
	       (w (/ (- hk hk0) 2))
	       (ff (- (/ (*g a r u) 2) (* w u))))
	  (make-instance
	   'ksstate
	   :alpha (- (* ff (/ (sin (* w0 s)) w0)))
	   :beta (* ff (/ (cos (* w0 s)) w0))
	   :e (* rm (scalar (*i v a)))
	   :tm rm
	   :constants constants))))))
