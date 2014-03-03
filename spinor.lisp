(in-package :bld-orbit)

;;; Spinor states

(defclass spinorstate ()
  ((u :initarg :u :documentation "Spinor of position: r = u e1 (revg u)")
   (duds :initarg :duds :documentation "Spinor derivative wrt s")
   (tm :initarg :tm :documentation "Time")
   (derived :accessor derived :initform (make-hash-table) :documentation "Derived values"))
  (:documentation "Spinor state"))

(defmethod print-object ((x spinorstate) stream)
  (format stream "#<SPINORSTATE :U ~a :DUDS ~a :TM ~a>" (slot-value x 'u) (slot-value x 'duds) (slot-value x 'tm)))

(defstatearithmetic spinorstate (u duds tm))

;;; Spinor specific derived slots & methods

(defderived xcart (s (x spinorstate) sc)
    "Cartesian state"
  (to-cartesian x s sc))

(defmethod r ((x spinorstate))
  (r (xcart 0 x nil)))

(defmethod v ((x spinorstate))
  (v (xcart 0 x nil)))

;; Other methods 

(defmethod energy-cb ((x spinorstate) cb)
  (with-slots (u duds) x
    (with-slots (mu) cb
      (/ (- (* 2 (norme2 duds)) mu) (norme2 u)))))

(defmethod energy ((x spinorstate) sc)
  "Keplerian orbit energy"
  (with-slots (u duds) x
    (with-slots (cb) sc
      (with-slots (mu) cb
	(/ (- (* 2 (norme2 duds)) mu) (norme2 u))))))

(defmethod time-of (s (x spinorstate))
  (slot-value x 'tm))

(defderived e (s (x spinorstate) sc)
    "Keplerian orbital energy"
  (with-slots (duds) x
    (with-slots (cb) sc
      (with-slots (mu) cb
	(with-derived (rm) s x sc
	  (/ (- (* 2 (norme2 duds)) mu) rm))))))

(defderived rm (s (x spinorstate) sc)
    "Orbit radius of spinor state"
  (with-slots (u) x
    (norme2 u)))

(defmethod eom (s (x spinorstate) sc)
  "Spinor equations of motion: 2 dU/ds - E U = f r U, dt/ds = |r|"
  (with-slots (u duds tm) x
    (with-slots (cb iframe afun) sc
      (with-slots (mu) cb
	(with-derived (a rm xcart e) s x sc
	(let* ((r (r xcart))
	       (v (v xcart)))
	  (make-instance
	   'spinorstate
	   :u duds
	   :duds (/ (+ (*g3 a r u) (* e u)) 2)
	   :tm rm)))))))

(defmethod to-cartesian ((x spinorstate) s sc)
  "Convert spinor state to cartesian coordinates given S and X"
  (with-slots (u duds tm) x
    (with-slots (iframe) sc
      (let* ((r (norme2 u))
	     (dudt (/ duds r)))
	(values
	 (make-instance
	  'cartstate
	  :r (spin (first iframe) u)
	  :v (graden (* 2 (*g3 dudt (first iframe) (revg u))) 1))
	 tm)))))

(defmethod to-spinor ((x cartstate) tm sc)
  "Convert cartesian state to spinor given time and X"
  (with-slots (r v) x
    (with-slots (iframe) sc
      (let ((u (recoverspinor (norme r) (rv-frame tm x sc) iframe)))
	(values
	 (make-instance
	  'spinorstate
	  :u u
	  :duds (* 0.5d0 (*g3 v u (first iframe)))
	  :tm tm)
	 0)))))

(defmethod to-initial-spinor ((x cartstate) tm sc)
  "Convert to an initial spinor state from cartesian state & time"
  (with-slots (r v) x
    (with-slots (iframe) sc
      (let ((u (recoverspinor (norme r) (rv-frame tm x sc) iframe)))
	(values
	 (make-instance
	  'spinorstate
	  :u u
	  :duds (/ (*g3 v u (first iframe)) 2)
	  :tm tm)
	 0)))))

(defmethod to-spinor (x s sc)
  "Default: try converting to cartesian first"
  (multiple-value-bind (xc tm) (to-cartesian x s sc)
    (to-spinor xc tm sc)))
