(in-package :cl-dftb)


(defun get-dftbp-api-version ()
  "get the api version of dftbplus"
  (with-foreign-objects ((major :int 1)
                         (minor :int 1)
                         (patch :int 1))
    (%dftbp-api major minor patch)
    (list :major (mem-ref major :int)
          :minor (mem-ref minor :int)
          :patch (mem-ref patch :int))
    )
  )

(defun is-instance-safe ()
  (if (%dftbp-is-instance-safe)
      (print "dftb instance is safe now")
      (print "not safe")
      ))

(defun get-num-atom (instance)
  (%dftbp-get-nr-atoms instance)
  )

(defun get-total-energy (instance energy)
  (%dftbp-get-energy instance energy)
  )

;(defclass Calculator ()
;  ((get-energy :initform 0.0f0 : accessor :get-energy)
;   (get-lattice-vec :initform (make-array 9 :initial-element nil))))
