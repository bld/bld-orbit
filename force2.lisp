;;; Force related functions

(in-package :bld-orbit2)

;;; Gravity

(defmethod gravity ((x cartstate))
  "Force of gravity"
  (with-slots (sc rm2 rm ru) x
    (with-slots (cb) sc
      (- (* (/ (mu cb) rm2) ru)))))

;;; Pointing functions

(defun fixed (x)
  "Return spacecraft pointing reference frame (and rotor) in inertial
space for a fixed attitude specified by RS in SC object"
  (with-slots (sc bframe rb) x
    (with-slots (iframe rs) sc
      (let ((rp (*g rb rs))) ; rotor of sail attitude in inertial space
	(values (new-frame rp iframe) rp)))))

(defun table-lookup (rs tm)
  "Lookup rotor for spacecraft orientation in table of (time rotor)"
  (second (find tm rs :test #'<= :key #'first)))

(defun table (x)
  "Return spacecraft pointing frame (and rotor) in intertial space for
a sail attitude specified by RS lookup table in SC object"
  (with-slots (sc bframe rb tm) x
    (with-slots (iframe rs) sc
      (let* ((rsi (table-lookup rs tm))
	     (rp (*g rb rsi)))
	(values (new-frame rp iframe) rp)))))

;;; Solar radiation pressure functions

(defun solar-pressure (x)
  "Solar radiation pressure experienced by state X"
  (let ((b (cb (sc x))))
    (/ (* (ls b) *au2*) (rm2 x) *c*)))

;;; Force functions

(defun ideal-sail-acc (x)
  "Force on an ideal solar sail of given state"
  (with-slots (sc ru p pframe) x
    (with-slots (area mass) sc
      (let ((n (third pframe)))
	(/ (* 2 p area (expt (scalar (*i ru n)) 2) n) mass)))))
