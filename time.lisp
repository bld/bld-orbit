(in-package :bld-orbit)

(defparameter *utc-j2000* (+ (encode-universal-time 55 58 11 1 1 2000 0) 0.816d0))

(defun utc-to-timestamp (utc)
  "Convert UTC time in seconds to a timestamp"
  (multiple-value-bind (sec rem) (floor utc)
    (universal-to-timestamp sec :nsec (round (* 1d9 rem)))))
