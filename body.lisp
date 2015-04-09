(in-package :bld-orbit)

(defclass body ()
  ((name :initarg :name :accessor name)
   (center :initarg :center :accessor center)
   (ref :initarg :ref :accessor ref)
   (mu :initarg :mu :accessor mu)))

(defmethod position-vector (et (b body) &key (observer (slot-value b 'center)) (ref (slot-value b 'ref)) (abcorr :none))
  "Position vector of a body given UTC time, body, and optionally:
:OBSERVER - name of observing body (default: body's CENTER slot)
:REF - name of reference frame (default: body's REF slot)
:ABCORR - aberration correction (default: :NONE)"
  (with-slots (name) b
    (let ((pv (spk-pos name et observer :ref ref :abcorr abcorr)))
      (make-instance 've3 :e1 (aref pv 0) :e2 (aref pv 1) :e3 (aref pv 2)))))

(defmethod gravity ((r ve3) (b body))
  "Gravitational acceleration given UTC, position vector of s/c relative to body, and body"
  (with-slots (mu) b
    (- (* (/ mu (expt (norme r) 3)) r))))
