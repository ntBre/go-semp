(defmacro deftest (name got-form want-form)
  (let ((got (gensym))
	(want (gensym)))
    `(defun ,name ()
       (let ((,got ,got-form)
	     (,want ,want-form))
	 (unless (equal ,got ,want)
	   (format t "~&FAIL: got ~&~a, wanted ~&~a~%" ,got ,want))))))

(deftest test-load-params
    (load-params "opt.out")
  (list (make-instance 'atom-type :label "H")))

(deftest test-parse-params
    (parse-params "DCore=2,3,1.6101301456,3.2139710000,0.0000000000,0.0000000000,0.0000000000")
  ())


(test-parse-params)

(test-load-params)
	  
