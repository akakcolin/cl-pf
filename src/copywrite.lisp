(defparameter *c++-copyright-body* "/*
Open Source License
Copyright (c) 2024
*/")

(defparameter *c++-copyright-tail* "/* -^- */")

(defparameter *lisp-copyright-body* "
;; Open Source License
;; Copyright (c) 2024 
")
(defparameter *lisp-copyright-tail* ";; -^-")

(defparameter *fortran-copyright-body* "
! Open Source License
! Copyright (c) 2024
")

(defparameter *fortran-copyright-tail* "!! -^-")


(defstruct copyright-notice
  header-lambda
  body
  tail)

(defparameter *copyright-notices* (make-hash-table :test #'equal))
(defparameter *c++-copyright-notice* (make-copyright-notice :header-lambda #'(lambda (fn) (format nil "/*~%    File: ~a~%*/~%" fn)) :body *c++-copyright-body* :tail *c++-copyright-tail*))
(defparameter *lisp-copyright-notice* (make-copyright-notice :header-lambda #'(lambda (fn) (format nil ";;;~%;;;    File: ~a~%;;;~%" fn)) :body *lisp-copyright-body* :tail *lisp-copyright-tail*))
(defparameter *fortran-copyright-notice* (make-copyright-notice :header-lambda #'(lambda (fn) (format nil "!!!~%!!!    File: ~a~%!!!~%" fn)) :body *fortran-copyright-body* :tail *fortran-copyright-tail*))
(progn
  (setf (gethash "c" *copyright-notices*) *c++-copyright-notice*)
  (setf (gethash "cpp" *copyright-notices*) *c++-copyright-notice*)
  (setf (gethash "h" *copyright-notices*) *c++-copyright-notice*)
  (setf (gethash "f90" *copyright-notices*) *fortran-copyright-notice*)
  (setf (gethash "lsp" *copyright-notices*) *lisp-copyright-notice*)
  (setf (gethash "lisp" *copyright-notices*) *lisp-copyright-notice*)
  )

(defun read-entire-file (pn)
  (with-open-file (fin pn :direction :input)
    (let ((seq (make-string (file-length fin))))
      (read-sequence seq fin)
      seq)))

(defun write-entire-file (contents pn)
  (with-open-file (fout pn :direction :output :if-exists :rename)
    (write contents :escape nil :stream fout)
    (terpri fout)))

(defun remove-copyright (contents notice)
  "Remove the copyright and return the contents.
If the copyright was removed the return the second value t"
  (let* ((cr-tail (copyright-notice-tail notice))
         (pos (search cr-tail contents)))
    (if pos
        (values (subseq contents (+ pos (length cr-tail))) t)
        (values contents nil))))

(defun insert-copyright (contents pathname notice saw-copyright)
  (let ((filename (format nil "~a.~a" (pathname-name pathname) (pathname-type pathname)))
        (body (copyright-notice-body notice))
        (tail (copyright-notice-tail notice)))
    (with-output-to-string (sout)
      (let ((header (funcall (copyright-notice-header-lambda notice) filename)))
        (format sout "~a" header)
        (format sout "~a~%" body)
        (format sout "~a" tail)
        (if (not saw-copyright) (format sout "~%"))
        (format sout "~a~%" contents)))))

(defun copywrite-one-file (path)
  (let* ((pn (pathname path))
         (filename (format nil "~a.~a" (pathname-name pn) (pathname-type pn)))
         (notice (gethash (pathname-type pn) *copyright-notices*))
         (contents (read-entire-file pn)))
    (multiple-value-bind (no-copyright saw-copyright)
        (remove-copyright contents notice)
      (let* ((result (insert-copyright no-copyright pn notice saw-copyright))
             (trimmed (string-right-trim '(#\newline #\linefeed #\nul #\space) result)))
        (write-entire-file trimmed pn)))))


(defun copywrite-all-files (pathname-list)
  (dolist (pn pathname-list)
    (format t "copywriting ~a~%" pn)
    (copywrite-one-file pn)))


(defparameter *all-files* (append
                           (directory #P"~/a.h")
                           (directory #P"~/src/*.cc")
                           (directory #P"~/src/*.f90")
                           (directory #P"~/src/*.lisp")))
(copywrite-all-files *all-files*)
