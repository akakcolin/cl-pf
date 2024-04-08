
(defun permutations (bag)
  "Return a list of all the permutations of the input."
  ;; If the input is nil, there is only one permutation:
  ;; nil itself
  (if (null bag)
      '(())
      ;; Otherwise, take an element, e, out of the bag.
      ;; Generate all permutations of the remaining elements,
      ;; And add e to the front of each of these.
      ;; Do this for all possible e to generate all permutations.
      (mapcan #'(lambda (e)
                  (mapcar #'(lambda (p) (cons e p))
                          (permutations
                            (remove e bag :count 1 :test #'eq))))
              bag)))

(defun permutations2 (lst)
  (if (null lst) '(())
    (let (res)
      (map nil
       (lambda (x)
         (setq res
	       (append res
		       (map 'list
			    (lambda (y) (cons x y))
			    (permutations (remove x lst))))))
       lst)
      res)))


(defun map-permutations (fun lst)
  (if (null lst) (funcall fun nil)
    (map nil
       (lambda (x)
         (map-permutations 
          (lambda (l) (funcall fun (cons x l))) 
          (remove x lst)))
       lst)))

;;get all possible combinations of the elements from lists contained on a list
(defun combinations (&rest lists)
  (if (endp lists)
      (list nil)
      (mapcan (lambda (inner-val)
                (mapcar (lambda (outer-val)
                          (cons outer-val
                                inner-val))
                        (car lists)))
              (apply #'combinations (cdr lists)))))
