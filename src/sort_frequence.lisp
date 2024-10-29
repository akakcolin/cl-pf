;; parse float number
(ql:quickload "parse-number")
;; grep freq 11.ss.yaml | awk '{print $2}' > infile

(defun load-simpledata-from-file (filename)
  (map 'list #'parse-number:parse-number (uiop:read-file-lines filename))
)

;; infile format just one data each line
;; 保存排序后的数据
(defun save-simple-sort-file (infile outfile)
  (with-open-file (str outfile
                       :direction :output
                       :if-exists :supersede
                       :if-does-not-exist :create)
  (dolist (x (sort (load-simpledata-from-file infile ) #'<))
        (format str "~A~%" x)
    ))
)
