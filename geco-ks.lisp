;;; Use GECO to optimize trajectories using Kustaanheimo-Stiefel orbital element equations of motion

(in-package :bld-orbit)

(defvar *t0-rel* 0d0 "Relative departure time from Earth (seconds)")

(defparameter *inertial-frame* *j2000* "Inirtial reference frame within which multivectors defined")

(defparameter *rs-table-length* 10 "Length of lookup table of sail/orbit frame rotors")

(defparameter *sf-arity* 100 "Number of final independant variable (SF) values")

(defparameter *sf-scale* 1d0 "Real valued maximum scale of SF")

(defparameter *tof-weight* 1d-11 "Weight on time of flight (in days) in cost function")

(defparameter *xf-weight* 1d0 "Weight on final state error from target in cost function")

(defparameter *target* *mars* "Target body")

(defparameter *origin* *earth* "Origin body")

(defparameter *star* *sun* "Star about which bodies orbit")

(defun sail-3d-Earth-Mars-template (rs)
  "Generate sail object given sail rotor (RS) lookup table"
  (let ((t0 (coerce (+ *t0-rel* (encode-universal-time 0 0 0 16 12 2013 0)) 'double-float)))
    (make-instance
     'sail
     :eom #'eom2
     :cb *star*
     :sun *star*
     :accfun #'sail-flat-ideal-acc
     :pointfun #'sail-frame-sun-table
     :lightness 0d0
     :area 1200d-6
     :mass 45d0
     :optical nil
     :basis *inertial-frame*
     :t0 0
     :tf (first (car (last rs)))
     :x0 (to-initial-ks (position-velocity *origin* t0) t0 (make-instance 'sail :basis *inertial-frame* :cb *star*))
     :rs rs
     :outfile nil)))

;; SF-CHROMOSOME

(defclass sf-chromosome (chromosome)
  ()
  (:documentation "SF final chromosome"))

(defmethod size ((self sf-chromosome))
  1)

(defmethod locus-arity ((self sf-chromosome) locus-index)
  *sf-arity*)

(defmethod loci-printable-form ((self sf-chromosome))
  (loci self))

;; RS-S-CHROMOSOME

(defclass rs-s-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of scalar rotor components"))

(defmethod size ((self rs-s-chromosome))
  *rs-table-length*)

(defmethod locus-arity ((self rs-s-chromosome) locus-index)
  101)

(defmethod loci-printable-form ((self rs-s-chromosome))
  (loci self))

;; RS-X-CHROMOSOME

(defclass rs-x-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual x) rotor components"))

(defmethod locus-arity ((self rs-x-chromosome) locus-index)
  101)

(defmethod size ((self rs-x-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-x-chromosome))
  (loci self))

;; RS-Y-CHROMOSOME

(defclass rs-y-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual y) rotor components"))

(defmethod locus-arity ((self rs-y-chromosome) locus-index)
  101)

(defmethod size ((self rs-y-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-y-chromosome))
  (loci self))

;; RS-Z-CHROMOSOME

(defclass rs-z-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual z) rotor components"))

(defmethod locus-arity ((self rs-z-chromosome) locus-index)
  101)

(defmethod size ((self rs-z-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-z-chromosome))
  (loci self))

;; Turn chromosome into RS table

(defmethod to-value ((sf-c sf-chromosome))
  "Convert TF chromosome to TF value"
  (* (/ (aref (loci sf-c) 0) *sf-arity*) *sf-scale*))

(defun chromosomes-to-rs-table (sf-c r-s-c r-x-c r-y-c r-z-c)
  "Turn chromosome into list of rotors suitable for table lookup."
  (loop with sf = (to-value sf-c) ; final time
     with ds = (/ sf *rs-table-length*) ; time step size
     ;; rotor components
     for r-s across (loci r-s-c) ; scalar
     for r-x across (loci r-x-c) ; (- (dual x))
     for r-y across (loci r-y-c) ; (- (dual y))
     for r-z across (loci r-z-c) ; (- (dual z))
     ;; accumulated times for table
     sum ds into s-accum
     collect (list s-accum
		   (unitg ; normalize rotor
		    (re3 :s (coerce r-s 'double-float) ; distribute rotor coefficients
			 :e2e3 (coerce r-x 'double-float)
			 :e1e3 (coerce (- r-y) 'double-float)
			 :e2e3 (coerce r-z 'double-float))))))

(defclass rs-table-organism (organism)
  ()
  (:documentation "RS lookup table organism"))

(defmethod chromosome-classes ((self rs-table-organism))
  '(sf-chromosome
    rs-s-chromosome
    rs-x-chromosome
    rs-y-chromosome
    rs-z-chromosome))

(defmethod pick-random-chromosome-index ((self rs-table-organism))
  (random (length (genotype self))))

(defclass rs-table-statistics (population-statistics)
  ()
  (:documentation "RS table population statistics"))

(defclass rs-table-population (generational-population maximizing-score-mixin)
  ()
  (:documentation "Population of RS tables to maximize scores"))

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
  (let* ((sail (organism-to-sail self)) ; sail object to propagate
	 (traj (propagate sail))) ; propagated trajectory
    (destructuring-bind (s0 x0) (first traj) ; initial state & independant variable
      (destructuring-bind (sf xf) (car (last traj)) ; final state & independant variable
	(let* ((t0 (time-of s0 x0)) ; initial time
	       (tf (time-of sf xf)) ; final time
	       (tof-cost (* (expt (/ (- tf t0) 24d0 60d0 60d0) 2) *tof-weight*))
	       (xf-cart (to-cartesian xf sf sail)) ; final cartesian sc state
	       (xf-cart-m (position-velocity *target* tf))) ; final mars cartesian state
	  (with-slots ((rf r) (vf v)) xf-cart ; sc final pos/vel
	    (with-slots ((rf-m r) (vf-m v)) xf-cart-m ; mars final pos/vel
	      (let ((xf-cost (* (+ (norme2 (- rf rf-m))
				   (norme2 (- vf vf-m)))
				*xf-weight*)))
		(setf (score self) (/ (+ tof-cost xf-cost)))))))))))

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
  0.0015)

(defmethod PROB-CROSS ((self rs-table-plan))
  "The probability of crossover for an organism."
  0.6)

#+null(defmethod OPERATE-ON-POPULATION
    ((plan rs-table-plan) old-population new-population 
     &AUX (new-organisms (organisms new-population))
       (p-cross (prob-cross plan))
       (p-mutate (+ p-cross (prob-mutate plan))))
  (let ((selector (stochastic-remainder-preselect old-population)))
    (do ((org1 (funcall selector) (funcall selector))
	 org2
	 (random# (geco-random-float 1.0) (geco-random-float 1.0))
	 (i 0 (1+ i)))
	((null org1))
      (cond
	((> p-cross random#)
	 (if (and (setq org2 (funcall selector))
		  (every #'(lambda (c-num)
			     (< 1 (hamming-distance (nth c-num (genotype org1)) (nth c-num (genotype org2)))))
			 '(0 1 2 3 4)))
	     (uniform-cross-organisms
	      org1 org2
	      (setf (aref new-organisms i)
		    (copy-organism org1 :new-population new-population))
	      (setf (aref new-organisms (1+ i))
		    (copy-organism org2 :new-population new-population)))
	     (progn
	       (setf (aref new-organisms i) (copy-organism-with-score org1 :new-population new-population))
	       (when org2
		 (setf (aref new-organisms (1+ i))
		       (copy-organism-with-score org2 :new-population new-population)))))
	 (incf i))
	((> p-mutate random#)
	 (mutate-organism
	  (setf (aref new-organisms i)
		(copy-organism org1 :new-population new-population))))
	(t
	 (setf (aref new-organisms i)
	       (copy-organism-with-score org1 :new-population new-population)))))))

(defmethod regenerate ((plan rs-table-plan) (old-pop rs-table-population))
  ;; Show statistics & generation
  (format t "Generation ~a: ~a~%" (generation-number (ecosystem old-pop)) (statistics old-pop))
  ;; Population to store tournament selection
  (let ((tour-pop (make-population (ecosystem old-pop) (class-of old-pop) :size (size old-pop)))
	(p-cross (prob-cross plan))
	(p-mutate (prob-mutate plan)))
    ;; Elitism: include best organism in new population
    (setf (aref (organisms tour-pop) 0) (copy-organism-with-score (best-organism old-pop) :new-population tour-pop))
    ;; Perform tournament selection on the rest of the population
    (dotimes (i (1- (size old-pop)))
      (setf (aref (organisms tour-pop) (1+ i))
	    (copy-organism (tournament-select-organism old-pop 2) :new-population tour-pop)))
    ;; Perform crossover on tour-pop and put in cross-pop
    (let ((cross-pop (make-population (ecosystem tour-pop) (class-of tour-pop) :size (size tour-pop))))
      ;; Copy elite from tour-pop to cross-pop
      (setf (aref (organisms cross-pop) 0) (copy-organism-with-score (aref (organisms tour-pop) 0) :new-population cross-pop))
      (dotimes (i (1- (size tour-pop)))
	;; When crossover happens
	(if (> p-cross (random 1.0))
	    (let* ((o1 (aref (organisms tour-pop) (1+ i)))
		   (i2 (pick-random-organism-index tour-pop))
		   (o2 (aref (organisms tour-pop) i2)))
	      (when (zerop i2)
		(cross-organisms 
		 o1 o2 
		 (setf (aref (organisms cross-pop) (1+ i))
		       (copy-organism o1 :new-population cross-pop))
		 (setf (aref (organisms cross-pop) i2)
		       (copy-organism o2 :new-population cross-pop)))))
	    ;; Else just copy over tour-pop organism
	    (setf (aref (organisms cross-pop) (1+ i)) (copy-organism (aref (organisms tour-pop) (1+ i)) :new-population cross-pop))))
      ;; Mutate cross-pop
      (dotimes (i (1- (size cross-pop)))
	(when (> p-mutate (random 1.0))
	  (mutate-organism (aref (organisms cross-pop) (1+ i)))))
      ;; Return new population
      cross-pop)))

(defvar *rs-table-ecosystem* nil "rs table ecosystem")

(defun find-rs-table (&key (pop-size 20) (evaluation-limit 400) (generation-limit 10))
  (setq *rs-table-ecosystem*
	(make-instance 
	 'ecosystem
	 :pop-class 'rs-table-population
	 :pop-size pop-size
	 :plan-class 'rs-table-plan
	 :generation-limit generation-limit
	 :evaluation-limit evaluation-limit))
  (evolve *rs-table-ecosystem*))

(defun find-best-organism (ecosystem)
  (let* ((fun (typecase (population ecosystem)
		(maximizing-score-mixin #'max)
		(minimizing-score-mixin #'min)
		(t (error "Max or Min score mixin not defined for population"))))
	 (orgs (organisms (population ecosystem))))
    (find (apply fun (map 'list #'score orgs)) orgs :key #'score)))

(defmethod organism-to-sail ((organism rs-table-organism) &aux (chromosomes (genotype organism)))
  (let* ((rs (apply #'chromosomes-to-rs-table chromosomes))
	 (sail (sail-3d-earth-mars-template rs)))
    sail))

(defun traj-to-planet (traj planet)
  (loop for (s x) in traj
     for tm = (time-of s x)
     collect (list tm (position-velocity planet tm))))

(defun write-traj-and-planets (ecosystem)
  (let* ((sail (organism-to-sail (find-best-organism ecosystem)))
	 (traj-data (propagate sail :hmax-factor 100)))
    (write-cart-traj "earth-mars.dat" (to-cartesian-traj traj-data sail))
    (write-cart-traj "earth.dat" (traj-to-planet traj-data *origin*))
    (write-cart-traj "mars.dat" (traj-to-planet traj-data *target*))))

(defun final-results (ecosystem)
  (let* ((sail (organism-to-sail (find-best-organism ecosystem)))
	 (traj (propagate sail :hmax-factor 100))
	 (e-data (traj-to-planet traj *origin*))
	 (m-data (traj-to-planet traj *target*))
	 (xf-m (car (last m-data)))
	 (c-traj (to-cartesian-traj traj sail))
	 (x0 (first c-traj))
	 (xf (car (last c-traj)))
	 (t0 (first x0))
	 (tf (first xf))
	 (days (/ (- tf t0) 24 60 60))
	 (xerr (- (second xf) (second xf-m))))
    (format t "Days: ~a~%Error: ~a~%" days xerr)
    (format t "State error: ~a~%" (- (second (car (last traj))) (slot-value *target* 'xks)))))

(defun evaluate-best (ecosystem)
  (let ((results (multiple-value-list (evaluate (find-best-organism ecosystem) (plan ecosystem))))
	(keys '(:total :time :alpha :beta :e :rm)))
    (plist-hash-table (loop for result in results
			 for key in keys
			 collect key
			 collect result))))
     
