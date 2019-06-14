#!/usr/local/bin/newlisp
; only use for shenmong job

(set 'pwd (append (real-path) "/"))


(define (gpuid)
 (amb 2 3 4 5 6 7))

(define (run-task dir item)
 (exec (string "sed -n '/pmemd/,$p' " (append dir item) " > " (append dir "pmemd.sh")))
 (if (file? (append dir "run.sh")) 
   (delete-file (append dir "run.sh")))  
 (println (append dir "run.sh"))
 (append-file (append dir "run.sh") (string "#!/usr/bin/env bash\n")) 
 (append-file (append dir "run.sh") 
 (set 'contens (read-file (append dir "pmemd.sh")))
 (append-file (append dir "run.sh") contens)
 (delete-file (append dir "pmemd.sh"))
 (exec (string "cd " dir " && sh run.sh"))
)

(define (search-tree dir file-name)
 (dolist (item (directory dir {^[^.]}))
   (if (directory? (append dir item))
    ; search the directory
    (search-tree (append dir  item "/") file-name)
    ; or process the file
    ;
    (if (find file-name item))
     (run-task dir item))))
    
(search-tree pwd "ti.sh")

(exit)

