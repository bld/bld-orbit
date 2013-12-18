;;; Use GECO to optimize trajectories

(in-package :bld-orbit)

(defun sail-3d-Earth-Mars-template ()
  (let ((t0 (coerce (encode-universal-time 0 0 0 16 12 2013 0) 'double-float)))
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
     :t0 t0
     :tf 0
     :x0 (position-velocity *earth* t0)
     :rs nil
     :outfile nil)))

(defparameter *rs-table-length* 10)

(defparameter *tm-limit* 1000)

(defparameter *tf-weight* 1d-4)

(defparameter *xf-weight* 1d3)

(defclass rs-times-chromosome (chromosome)
  ()
  (:documentation "RS lookup table times sequence chromosome"))

(defmethod size ((self rs-times-chromosome))
  *rs-table-length*)

(defmethod locus-arity ((self rs-times-chromosome) locus-index)
  *tm-limit*)

(defmethod loci-printable-form ((self rs-times-chromosome))
  (loci self))

(defclass rs-s-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of scalar quaternion components"))

(defmethod size ((self rs-s-chromosome))
  *rs-table-length*)

(defmethod locus-arity ((self rs-s-chromosome) locus-index)
  101)

(defmethod loci-printable-form ((self rs-s-chromosome))
  (loci self))

(defclass rs-x-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual x) quaternion components"))

(defmethod locus-arity ((self rs-x-chromosome) locus-index)
  101)

(defmethod size ((self rs-x-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-x-chromosome))
  (loci self))

(defclass rs-y-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual y) quaternion components"))

(defmethod locus-arity ((self rs-y-chromosome) locus-index)
  101)

(defmethod size ((self rs-y-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-y-chromosome))
  (loci self))

(defclass rs-z-chromosome (chromosome)
  ()
  (:documentation "RS lookup table of (dual z) quaternion components"))

(defmethod locus-arity ((self rs-z-chromosome) locus-index)
  101)

(defmethod size ((self rs-z-chromosome))
  *rs-table-length*)

(defmethod loci-printable-form ((self rs-z-chromosome))
  (loci self))

(defun chromosomes-to-rs-table (t0 tm-c s-c x-c y-c z-c)
  "Turn chromosome into list of rotors suitable for table lookup in"
  (loop for tm across (loci tm-c)
     for s across (loci s-c)
     for x across (loci x-c)
     for y across (loci y-c)
     for z across (loci z-c)
     sum tm into tm-accum
     collect (list (+ t0 (* tm-accum 24d0 60d0 60d0))
		   (unitg 
		    (re3 :s (coerce s 'double-float)
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

(defmethod evaluate ((self rs-table-organism) (plan rs-table-plan) &aux (chromosomes (genotype self)))
  (let* ((sail (sail-3d-earth-mars-template))
	 (rs (apply #'chromosomes-to-rs-table (slot-value sail 't0) chromosomes))
	 (tf (first (car (last rs)))))
    (setf (slot-value sail 'rs) rs)
    (setf (slot-value sail 'tf) tf)
    (let* ((traj (propagate sail))
	   (xf (second (car (last traj))))
	   (tf-cost (* *tf-weight* tf))
	   (x-target (position-velocity *mars* tf))
	   (xf-diff (- xf x-target))
	   (xf-cost (* *xf-weight*
		       (+ (norme (slot-value xf-diff 'r))
			  (norme (slot-value xf-diff 'v))))))
      (setf (score self) (+ tf-cost xf-cost)))))

(defmethod REGENERATE ((plan rs-table-plan) (old-pop rs-table-population)
		       &AUX (new-pop (make-population (ecosystem old-pop)
						      (class-of old-pop)
						      :size (size old-pop))))
  "Create a new generation from the previous one, and record statistics."
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
