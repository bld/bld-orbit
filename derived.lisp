(in-package :bld-orbit)

;;; Derived parameters used in equation of motion calculations
;; Stored in DERIVED slot of state variables for re-use
;; Calculated & set upon first call

(defderived xcart (tm x sc)
    "Cartesian state"
  (to-cartesian x tm sc))

(defmethod xcart (tm (x cartstate) sc)
  "Just return x for states already cartesian"
  x)

(defderived r-cb (tm x sc)
    "Inertial position of central body"
  (with-slots (cb) sc
    (r (position-velocity cb (time-of tm x)))))

(defderived r-sc (tm x sc)
    "Inertial position of spacecraft"
  (+ (r-cb tm x sc)
     (r (xcart tm x sc))))

(defderived r-sun (tm x sc)
    "Inertial position of sun"
  (with-slots (sun) sc
    (r (position-velocity sun (time-of tm x)))))

(defderived r-sc-sun (tm x sc)
    "Position of spacecraft relative to sun"
  (with-derived (r-sc r-sun) tm x sc
    (- r-sc r-sun)))

(defderived ru-sc-sun (tm x sc)
    "Unit position vector of spacecraft relative to sun"
  (unitg (r-sc-sun tm x sc)))

(defderived rotor-o (tm x sc)
    "Rotor of orbit frame"
  (with-slots (iframe) sc
    (with-derived (oframe) tm x sc
      (recoverrotor oframe iframe))))

(defderived rotor-p (tm x sc)
    "Rotor of pointing frame"
  (with-slots (rs) sc
    (*g (rotor-o tm x sc) rs)))

(defderived oframe (tm x sc)
    "Orbit frame"
  (with-slots (ofun) sc
    (funcall ofun tm x sc)))

(defderived pframe (tm x sc)
  "Pointing frame"
  (with-slots (pfun) sc
    (funcall pfun tm x sc)))

(defderived a (tm x sc)
  "External acceleration"
  (with-slots (afun) sc
    (funcall afun tm x sc)))

(defderived ru (tm x sc)
  "Unit vector of position wrt to CB"
  (unitg (r (xcart tm x sc))))

(defderived rm2 (tm x sc)
  "Radius squared wrt CB"
  (norme2 (r (xcart tm x sc))))

(defderived rm (tm x sc)
  "Radius from CB"
  (with-derived (rm2) tm x sc
    (sqrt rm2)))

(defderived p (tm x sc)
  "Solar pressure"
  (with-slots (spfun) sc
    (funcall spfun tm x sc)))

(defderived g (tm x sc)
  "Central body gravity"
  (with-slots (gfun) sc
    (funcall gfun tm x sc)))

(defderived n (tm x sc)
  "Normal vector as function of PFRAME"
  (with-slots (nfun) sc
    (with-derived (pframe) tm x sc
      (funcall nfun pframe))))
