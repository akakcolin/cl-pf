#/usr/bin/env newlisp

(set 'dirlist (list "0.0" "0.0479" "0.1151" "0.2063" "0.3161" "0.4374" "0.5626" "0.6839" "0.7937" "0.8849" "0.9521" "1.0"))

(define (run_task a i num_frame num)
   (println (string "cd " a " && python parse_amber.py " i " " num_frame " " num))
   (exec (string "cd " a " && python parse_amber.py " i " " num_frame " " num))
)

(define (run_rerun item i numjob)
  (set 'command (string "cd " item " && parallel -j " numjob " -a " i "_pmemd_command"))
  (println command)
  (exec command)
)

(define (string-id i)
  (if 
    (< i 10) (string "00" i)
    (> i 9) (string "0" i)
    (> i 99) (string i)
    )
  )

;; run remd gpu
(exec (string "mpirun -np 12 pmemd.cuda.MPI -ng 12 -groupfile group.dat"))

;; copy mdcrd to each dir
(set 'i 0)
(dolist (item dirlist)
  (println item)
  (exec (string "cp mdcrd." (string-id i) " " item "/mdcrd." i))
  (inc i)
)

;; parse coordinates
(set 'i 0)
(dolist (item dirlist)
  (exec (string "cp parse_amber.py " item ))
  (spawn 'f1 (run_task item i 12 4000))
  (inc i)
)
(println "\n")
(until (sync 1000) (println "."))
(println "\n")

;; create rerun input
(set 'i 0)
(dolist (item dirlist)
  (println item)
  (exec (string "cp " item "/ti.in.rep.* " item  "/ti.in.rep." i))
  (exec (string "cp " item "/prmtop.* "  " ." ))
  (inc i)
)
(exec (string "python rerun_script.py"))

(dolist (item dirlist)
 (println item)
 (exec (string "cp rerun.in.* " item))
 (exec (string "cp prmtop.* " item))
)

;; run rerun 
(set 'i 0)
(dolist (item dirlist)
  (run_rerun item i 20)
  (inc i)
)


(exit)
