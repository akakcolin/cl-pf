;; get state
(ql:quickload "uiop")
(defun get-slume-stat ()
  (uiop:run-program (concatenate 'string
                                 "squeue "
                                 " ") :force-shell t :output t)
  )

;; Start the inferior shell, with input and output streams

;; Write a line to the shell
;(progn
;  (defparameter *shell* (uiop:launch-program "bash" :input :stream :output :stream))
;  (write-line "cd /home/dd/FILES/ && ./run_sumpoi.sh /usr/local/MATLAB/MATLAB_Runtime/v96"
;              (uiop:process-info-input *shell*))
  ;; Flush stream
;  (force-output (uiop:process-info-input *shell*))y
  ;(read-line (uiop:process-info-output *shell*))


;  (let ((stream (uiop:process-info-output *shell*)))
;    (loop while (listen stream) do
      ;; Characters are immediately available
;      (princ (read-line stream))
;      (terpri)))
;  (format t "AWE")
  ;(uiop:close-streams *shell*)
;  (uiop:process-alive-p *shell*))
