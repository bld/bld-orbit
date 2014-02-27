#|
Orbital Mechanics Library
=========================

Includes:
* Cartesian coordinate equations of motion
* Kustaanheimo-Stiefel equations of motion cast into geometric algebra by Hestenes
* Simple Keplerian trajectories
* Solar sail trajectories
|#

(in-package :bld-orbit)

;;; Utility functions

(defgeneric energy-cb (x cb) (:documentation "Specific orbital energy from state & central body"))

(defgeneric energy (x sc) (:documentation "Keplerian orbit energy"))

(defgeneric time-of (s x) (:documentation "Universal time"))

(defgeneric eom (s x param) (:documentation "Equations of motion given independent variable S & state X"))

(defgeneric to-cartesian (x s sc) (:documentation "Convert a state to cartesian state"))

(defgeneric to-spinor (x s sc) (:documentation "Convert a state to spinor state"))

(defgeneric to-initial-spinor (x s sc) (:documentation "Convert to an initial spinor state"))

(defgeneric to-ks (x s sc) (:documentation "Convert to a Kustaanheimo-Stiefel orbit element state"))

(defgeneric to-initial-ks (x s sc) (:documentation "Generate initial Kustaanheimo-Stiefel orbit element state from another state & corresponding central body"))

