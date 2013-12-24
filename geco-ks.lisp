;;; Use GECO to optimize trajectories using Kustaanheimo-Stiefel orbital element equations of motion

(in-package :bld-orbit)

(defun sail-3d-Earth-Mars-template (rs)
  (let ((t0 (coerce (+ *t0-rel* (encode-universal-time 0 0 0 16 12 2013 0)) 'double-float)))
    (make-instance
     'sail
     :eom #'eom2
     :cb *sun*
     :sun *sun*
     :accfun #'sail-flat-ideal-acc
     :pointfun #'sail-frame-sun-table
     :lightness 0d0
     :area 1200d-6
     :mass 45d0
     :optical nil
     :basis *j2000*
     :t0 0
     :tf (first (car (last rs)))
     :x0 (to-initial-ks (position-velocity *earth* t0) t0 (make-instance 'sail :basis *j2000* :cb *sun*))
     :rs rs
     :outfile nil)))

(defvar *t0-rel* 0d0)

(defparameter *rs-table-length* 10)

(defparameter *tm-limit* 200)

(defparameter *tm-scale* (/ (* 24 60 60) *au*))

(defparameter *tf-weight* 1d-5)

(defparameter *rf-weight* 1d-5)

(defparameter *xf-weight* (alexandria:plist-hash-table '(:alpha 1d0 :beta 1d0 :e 1d2)))

;; RS-TIMES-CHROMOSOME

(defclass rs-times-chromosome (chromosome)
  ()
  (:documentation "RS lookup table times (as 's') chromosome"))

(defmethod size ((self rs-times-chromosome))
  *rs-table-length*)

(defmethod locus-arity ((self rs-times-chromosome) locus-index)
  *tm-limit*)

(defmethod loci-printable-form ((self rs-times-chromosome))
  (loci self))

;; RS-S-CHROMOSOME

(defclass rs-s-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of scalar quaternion components"))

(defmethod size ((self rs-s-chromosome))
  *rs-table-length*)

(defmethod locus-arity ((self rs-s-chromosome) locus-index)
  101)

(defmethod loci-printable-form ((self rs-s-chromosome))
  (loci self))

;; RS-X-CHROMOSOME

(defclass rs-x-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual x) quaternion components"))

(defmethod locus-arity ((self rs-x-chromosome) locus-index)
  101)

(defmethod size ((self rs-x-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-x-chromosome))
  (loci self))

;; RS-Y-CHROMOSOME

(defclass rs-y-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual y) quaternion components"))

(defmethod locus-arity ((self rs-y-chromosome) locus-index)
  101)

(defmethod size ((self rs-y-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-y-chromosome))
  (loci self))

;; RS-Z-CHROMOSOME

(defclass rs-z-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual z) quaternion components"))

(defmethod locus-arity ((self rs-z-chromosome) locus-index)
  101)

(defmethod size ((self rs-z-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-z-chromosome))
  (loci self))

;; Turn chromosome into RS table

(defun chromosomes-to-rs-table (tm-c s-c x-c y-c z-c)
  "Turn chromosome into list of rotors suitable for table lookup"
  (loop for tm across (loci tm-c)
     for s across (loci s-c)
     for x across (loci x-c)
     for y across (loci y-c)
     for z across (loci z-c)
     sum (* tm *tm-scale*) into tm-accum ; running sum
     collect (list tm-accum
		   (unitg ; normalize rotor
		    (re3 :s (coerce s 'double-float) ; distribute rotor coefficients
			 :e2e3 (coerce x 'double-float)
			 :e1e3 (coerce (- y) 'double-float)
			 :e2e3 (coerce z 'double-float))))))

(defclass rs-table-organism (organism)
  ()
  (:documentation "RS lookup table organism"))

(defmethod chromosome-classes ((self rs-table-organism))
  '(rs-times-chromosome
    rs-s-chromosome
    rs-x-chromosome
    rs-y-chromosome
    rs-z-chromosome))

(defclass rs-table-statistics (population-statistics)
  ()
  (:documentation "RS table population statistics"))

(defclass rs-table-population (generational-population minimizing-score-mixin)
  ()
  (:documentation "Population of RS tables to minimize scores"))

(defmethod organism-class ((self rs-table-population))
  'rs-table-organism)

(defmethod population-statistics-class ((self rs-table-population))
  'rs-table-statistics)

(defclass rs-table-plan (genetic-plan)
  ()
  (:documentation "RS lookup table plan"))

(defmethod propagate-organism ((o rs-table-organism))
  "Turn organism to sail object & propagate"
  (propagate (organism-to-sail o)))

(defmethod evaluate ((self rs-table-organism) (plan rs-table-plan) &aux (chromosomes (genotype self)))
  "Evaluate organism based on minimizing (- tf t0) and difference from Mars state at tf"
  (let* ((sail (organism-to-sail self)))
    (with-slots ((s0 t0) (sf tf)) sail
      (with-slots ((xks-m xks)) *mars*
	(let* ((traj (propagate sail)) ; propagated trajectory
	       (x0 (second (first traj))) ; final state
	       (xf (second (car (last traj)))) ; final state
	       (t0 (slot-value x0 'tm))
	       (tf (slot-value xf 'tm))
	       (tm-cost (* *tf-weight* (- tf t0))) ; Traval time cost
	       (xf-cart (to-cartesian xf sf sail)) ; final sc position/velocity
	       (xf-cart-m (position-velocity *mars* tf)) ; final mars position/velocity
	       (rf (slot-value xf-cart 'r))
	       (rf-m (slot-value xf-cart-m 'r)))
	  (with-slots (alpha beta e tm) xf
	    (with-slots ((alpha-m alpha) (beta-m beta) (e-m e)) xks-m
	      (let* ((dalpha (- alpha alpha-m))
		     (dbeta (- beta beta-m))
		     (alpha-cost (* (norme dalpha) (gethash :alpha *xf-weight*))) ; match alpha
		     (beta-cost (* (norme dbeta) (gethash :beta *xf-weight*))) ; match beta
		     (e-cost (* (abs (- e e-m)) (gethash :e *xf-weight*))) ; match energy
		     (rm-cost (* *rf-weight* (norme (- rf rf-m))))) ; minimize radius to target body
		(values (setf (score self) (+ tm-cost alpha-cost beta-cost e-cost rm-cost))
			tm-cost
			alpha-cost
			beta-cost
			e-cost
			rm-cost)))))))))

(defmethod REGENERATE ((plan rs-table-plan) (old-pop rs-table-population)
		       &AUX (new-pop (make-population (ecosystem old-pop)
						      (class-of old-pop)
						      :size (size old-pop))))
  "Create a new generation from the previous one, and record statistics."
  ;; Show statistics & generation
  (format t "Generation ~a: ~a~%" (generation-number (ecosystem old-pop)) (statistics old-pop))
  (setf (ecosystem new-pop) (ecosystem old-pop))
  ;; selectively reproduce, crossover, and mutate
  (operate-on-population plan old-pop new-pop)
  new-pop)

(defmethod PROB-MUTATE ((self rs-table-plan))
  "This is the probability of mutating an organism, not a single locus as is often used."
  0.03)

(defmethod PROB-CROSS ((self rs-table-plan))
  "The probability of crossover for an organism."
  0.7)

(defmethod OPERATE-ON-POPULATION
    ((plan rs-table-plan) old-population new-population 
     &AUX (new-organisms (organisms new-population))
       (p-cross (prob-cross plan))
       (p-mutate (+ p-cross (prob-mutate plan)))
       (orphan (make-instance (organism-class old-population)))) ; not in any population
  (let ((selector (stochastic-remainder-preselect old-population)))
    (do ((org1 (funcall selector) (funcall selector))
	 org2
	 (random# (geco-random-float 1.0) (geco-random-float 1.0))
	 (i 0 (1+ i)))
	((null org1))
      (cond
	((> p-cross random#)
	 (setf org2 (pick-random-organism old-population))
	 (if (every #'(lambda (c-num)
			(< 1 (hamming-distance (nth c-num (genotype org1)) (nth c-num (genotype org2)))))
		    '(0 1 2 3 4))
	     (uniform-cross-organisms org1 org2 (setf (aref new-organisms i) (copy-organism org1 :new-population new-population)) orphan)
	     (setf (aref new-organisms i) (copy-organism-with-score org1 :new-population new-population))))
	((> p-mutate random#)
	 (mutate-organism
	  (setf (aref new-organisms i)
		(copy-organism org1 :new-population new-population))))
	(t
	 (setf (aref new-organisms i)
	       (copy-organism-with-score org1 :new-population new-population)))))))

(defvar *rs-table-ecosystem* nil "rs table ecosystem")

(defun find-rs-table (&key (pop-size 20) (evaluation-limit 400))
  (setq *rs-table-ecosystem*
	(make-instance 
	 'ecosystem
	 :pop-class 'rs-table-population
	 :pop-size pop-size
	 :plan-class 'rs-table-plan
	 :evaluation-limit evaluation-limit))
  (evolve *rs-table-ecosystem*))

(defun find-min-organism (ecosystem)
  (let ((orgs (organisms (population ecosystem))))
    (find (apply #'min (map 'list #'score orgs)) orgs :key #'score)))

(defmethod organism-to-sail ((organism rs-table-organism) &aux (chromosomes (genotype organism)))
  (let* ((rs (apply #'chromosomes-to-rs-table chromosomes))
	 (sail (sail-3d-earth-mars-template rs)))
    sail))

(defun traj-to-planet (traj planet)
  (loop for (s x) in traj
     for tm = (time-of s x)
     collect (list tm (position-velocity planet tm))))

(defun write-min-traj-and-planets (ecosystem)
  (let* ((sail (organism-to-sail (find-min-organism ecosystem)))
	 (traj-data (propagate sail :hmax-factor 100)))
    (write-cart-traj "earth-mars.dat" (to-cartesian-traj traj-data sail))
    (write-cart-traj "earth.dat" (traj-to-planet traj-data *earth*))
    (write-cart-traj "mars.dat" (traj-to-planet traj-data *mars*))))

(defun final-results (ecosystem)
  (let* ((sail (organism-to-sail (find-min-organism ecosystem)))
	 (traj (propagate sail :hmax-factor 100))
	 (e-data (traj-to-planet traj *earth*))
	 (m-data (traj-to-planet traj *mars*))
	 (xf-m (car (last m-data)))
	 (c-traj (to-cartesian-traj traj sail))
	 (x0 (first c-traj))
	 (xf (car (last c-traj)))
	 (t0 (first x0))
	 (tf (first xf))
	 (days (/ (- tf t0) 24 60 60))
	 (xerr (- (second xf) (second xf-m))))
    (format t "Days: ~a~%Error: ~a~%" days xerr)
    (format t "State error: ~a~%" (- (second (car (last traj))) (slot-value *mars* 'xks)))))

(defun evaluate-min (ecosystem)
  (let ((results (multiple-value-list (evaluate (find-min-organism ecosystem) (plan ecosystem))))
	(keys '(:total :time :alpha :beta :e :rm)))
    (plist-hash-table (loop for result in results
			 for key in keys
			 collect key
			 collect result))))
     