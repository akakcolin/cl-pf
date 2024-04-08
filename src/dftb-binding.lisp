(in-package :cl-dftb)

(defgeneric wrap-pointer (pointer)
  (:method ((pointer null)) nil))

(defcstruct DftbPlusAtomList (pDftbPlusAtomList :pointer))
(defcstruct Dftbplus (pDftbplus :pointer))
(defcstruct DftbplusInput)

(defctype p-DftbPlus
  :pointer (:struct DftbPlus))

(defctype p-DftbPlusAtomList
  :pointer (:struct DftbPlusAtomList))

(defctype p-DftbPlusInput
  :pointer (:struct DftbPlusInput))

(defctype size :unsigned-int)

(define-foreign-type c-float ()
  ()
  (:actual-type :float)
  (:simple-parser c-float))

(defmethod translate-to-foreign (value (type c-float)) (float value 0f0))

(defmethod translate-into-foreign-memory (value (type c-float) pointer)
  (translate-into-foreign-memory (translate-to-foreign value type)
                                 (make-instance 'cffi::foreign-built-in-type
                                                :type-keyword :float)
                                 pointer))

(define-foreign-type c-double ()
  ()
  (:actual-type :double)
  (:simple-parser c-double))

(defmethod translate-to-foreign (value (type c-double)) (float value 0d0))

(defmethod translate-into-foreign-memory (value (type c-double) pointer)
  (translate-into-foreign-memory (translate-to-foreign value type)
                                 (make-instance 'cffi::foreign-built-in-type
                                                :type-keyword :double)
                                 pointer))

(define-foreign-type c-int ()
  ()
  (:actual-type :int)
  (:simple-parser c-int))
(defmethod translate-to-foreign (value (type c-int)) (truncate value))
(defmethod translate-into-foreign-memory (value (type c-int) pointer)
  (translate-into-foreign-memory (translate-to-foreign value type)
                                 (make-instance 'cffi::foreign-built-in-type
                                                :type-keyword :int)
                                 pointer))

(define-foreign-type c-ptr ()
  ()
  (:actual-type :pointer)
  (:simple-parser c-ptr))
(defmethod translate-to-foreign (value (type c-ptr))
  (if (or (eql value 0) (eq value nil))
      (null-pointer)
      value))
(defmethod translate-into-foreign-memory (value (type c-ptr) pointer)
  (translate-into-foreign-memory (translate-to-foreign value type)
                                 (make-instance 'cffi::foreign-built-in-type
                                                :type-keyword :pointer)
                                 pointer))


(defcfun (%dftbp-api "dftbp_api") :void
  (major c-ptr)
  (minor c-ptr)
  (patch c-ptr)
  )

(defcfun ("dftbp_is_instance_safe" %dftbp-is-instance-safe) :boolean)

(defcfun ("dftbp_init" dftbp-init) :void
         (instance p-DftbPlus)
         (outputfilename :string))

(defcfun ("dftbp_init_mpi" %dftbp-init-mpi) :void
  (instance p-DftbPlus)
  (outputfilename :string)
  (mpiComm :int))

(defcfun (%dftbp-final "dftbp_final") :void
         (instance p-DftbPlus))

(defcfun (%dftbp-get-atom-list "dftbp_get_atom_list") :void
  (atomListHandler p-DftbPlusAtomList)
  (nAtomC :pointer)
  (nSpeciesC :pointer)
  (elementC  :string)
  (species :pointer))

(defcfun (%dftbp-get-input-from-file "dftbp_get_input_from_file") :void
         (instance p-DftbPlus)
         (filename :string)
         (input p-DftbPlusInput))

(defcfun (%dftbp-process-input "dftbp_process_input") :void
         (instance p-DftbPlus)
         (input p-DftbPlusInput))

(defcfun ("dftbp_set_external_potential" %dftbp-set-external-potential) :void
  (instance p-DftbPlus)
  (extpot :pointer)
  (extpotgrad :pointer))

(defcfun ("dftbp_register_ext_pot_generator" %dftbp-register-ext-pot-generator) :void
  (instance p-DftbPlus)
  (refptr :pointer)
  (extpot :pointer)
  (extpotgrad :pointer))

(defcfun ("dftbp_set_coords" %dftbp-set-coords) :void
  (instance p-DftbPlus)
  (coords :pointer))

(defcfun ("dftbp_set_coords_and_lattice_vecs" %dftbp-set-coords-and-lattice-vecs) :void
  (instance p-DftbPlus)
  (coords :pointer)
  (latvecs :pointer))

(defcfun ("dftbp_set_coords_lattice_origin" %dftbp-set-coords-lattice-origin) :void
  (instance p-DftbPlus)
  (coords :pointer)
  (latvecs :pointer)
  (origin :pointer))

(defcfun ("dftbp_get_nr_atoms" %dftbp-get-nr-atoms) :int
  (instance p-DftbPlus))

(defcfun ("dftbp_nr_kpoints" %dftbp-nr-kpoints) :int
  (instance p-DftbPlus))

(defcfun ("dftbp_get_energy" %dftbp-get-energy) :void
  (instance p-DftbPlus)
  (mermin-energy :pointer))

(defcfun ("dftbp_get_gradients" %dftbp-get-gradients) :void
  (instance p-DftbPlus)
  (gradients :pointer))

(defcfun ("dftbp_get_nr_orbitals" %dftbp-get-nr-orbitals) :void
  (instace p-DftbPlus)
  (nOrbitals :pointer))

(defcfun ("dftbp_get_masses" %dftbp-get-masses) :void
  (instance p-DftbPlus)
  (masses :pointer))

(defcfun ("dftbp_get_stress_tensor" %dftbp-get-stress-tensor) :void
  (instance p-DftbPlus)
  (stresstensor :pointer))

(defcfun ("dftbp_get_gross_charges" %dftbp-get-gross-charges) :void
  (instance p-DftbPlus)
  (charges :pointer))

(defcfun ("dftbp_get_cm5_charges" %dftbp-get-cm5-charges) :void
  (instace p-DftbPlus)
  (charges :pointer))

(defcfun ("dftbp_get_elstat_potential" %dftbp-get-elstat-potential) :void
  (instance p-DftbPlus)
  (nLocations :int)
  (potential :pointer)
  (location :pointer))
