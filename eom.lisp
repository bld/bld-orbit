(in-package :bld-orbit)

(defun gravity2 (r b)
  (with-slots (mu) b
    (- (* (/ mu (norme2 r)) (unitg r)))))

(defgeneric dxds (s x sc f g-cb))

(defmethod eom2 (s x sc)
  "Equations of motion for any "
  (with-slots (cb sun pointfun accfun) sc
    ;; position & velocity of SC relative to central body
    (with-slots (r v) (to-cartesian x s sc)
      (let* ((tm (time-of s x)) ; UTC
	     (r-cb (slot-value (position-velocity cb tm) 'r)) ; State of central body
	     (r-sc (+ r-cb r)) ; position of SC in system
	     (r-sun (slot-value (position-velocity sun tm) 'r))
	     (rsc-sun (- r-sc r-sun)) ; position of SC wrt sun
	     (g-cb (gravity2 r cb)) ; gravity of central body
	     (sframe (funcall pointfun s x sc)) ; sail frame
	     (f-sun (funcall accfun rsc-sun sframe sc sun))) ; sun sail force
	(dxds s x sc f-sun g-cb))))) ; derivative of state

(defmethod dxds (tm (x cartstate) sc f g-cb)
  (with-slots (r v) x
    (make-instance
     'cartstate
     :r v
     :v (+ f g-cb))))

(defmethod dxds (s (x ksstate) sc f g-cb)
  (with-slots (alpha beta e tm) x
    (let* ((xspin (to-spinor x s sc))
	   (xcart (to-cartesian xspin s sc))
	   (r (slot-value xcart 'r))
	   (v (slot-value xcart 'v))
	   (x0 (slot-value sc 'x0))
	   (e0 (slot-value x0 'e))
	   (hk0 (- e0))
	   (w0 (sqrt (/ hk0 2)))
	   (rm (norme r))
	   (u (slot-value xspin 'u))
	   (hk (- e))
	   (w (/ (- hk hk0) 2))
	   (ff (- (/ (*g f r u) 2) (* w u))))
      (make-instance
       'ksstate
       :alpha (- (* ff (/ (sin (* w0 s)) w0)))
       :beta (* ff (/ (cos (* w0 s)) w0))
       :e (* rm (scalar (*i v f)))
       :tm rm))))
