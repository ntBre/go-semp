(in-package :cl-user)

(defpackage :com.bwestbro.semp
  (:use :cl :cl-user :uiop :com.bwestbro.strings))

(in-package :com.bwestbro.semp)

(defparameter *charge* 0)

(defparameter *spin* 1)

(defclass atom-type ()
  ((label :initarg :label
	  :accessor label)
   (param-names :initarg :param-names
		:accessor :param-names
		:initform ())
   (param-values :initarg :param-values
		:accessor :param-values
		:initform ())))

(defgeneric push-names (atom-type name-list))

(defmethod push-names ((atom atom-type) name-list)
  (with-slots (param-names) atom
    (setf param-names (append param-names name-list))))

(defmethod push-values ((atom atom-type) value-list)
  (with-slots (param-values) atom
    (setf param-values (append param-values value-list))))

(defun print-atom (atom &optional (stream t))
  (with-slots (label param-names param-values) atom
    (format stream "~& ****~%~2@a~%" label)
    (loop for name in param-names
	  for value in param-values
	  with last = ""
	  do
	     (if (string= name last)
		 (format stream ",~a" value)
		 (format stream "~&~a~a~a" name
			 (if (contains name "DCore") "," "=") value))
	     (setf last name))))

(defun derived-param-p (param)
  "Check if a parameter is derived from other parameters"
  (member (car (uiop:split-string param :separator '(#\=)))
	  (list "PQN" "NValence" "DDN" "KON" "EISol")
	  :test #'string-equal))

(defun parse-params (line)
  "Parse a line of semi-empirical parameters"
  (let ((fields (remove-if #'derived-param-p (fields line)))
	names values)
    (dolist (f fields)
      (destructuring-bind (name value) 
	  (uiop:split-string f :separator '(#\=))
	(when (string= name "DCore")
	  (let ((fields (uiop:split-string value :separator '(#\,))))
	    (setf name (strcat name "=" (pop fields) "," (pop fields)))
	    (setf value (join fields ","))))
	(loop for f in (uiop:split-string value :separator '(#\,))
	      do (push name names)
		 (push f values))))
    (values (nreverse names) (nreverse values))))

(defun load-params (filename)
  "Load semi-empirical parameters from FILENAME, returning them as a
list of ATOM-TYPEs"
  (let ((in-labl nil)
	(in-params nil)
	(atoms ()))
    (with-open-file (infile filename :direction :input)
      (loop for line = (read-line infile nil nil)
	    while line
	    do (cond
		 ((string-equal line " ****")
		  (setf in-labl t))
		 ((and in-labl (string-equal line " "))
		  (setf in-labl nil
			in-params nil))
		 (in-labl
		  (push (make-instance 'atom-type :label
				       (read-from-string line))
			atoms)
		  (setf in-labl nil
			in-params t))
		 (in-params
		  (multiple-value-bind (names values) (parse-params line)
		    (push-names (car atoms) names)
		    (push-values (car atoms) (mapcar #'parse-float values)))))))
    (nreverse atoms)))

(defun dump-params (atoms &optional (filename "params.dat"))
  (with-open-file (outfile filename :direction :output :if-exists :supersede)
    (mapcar #'(lambda (a) (print-atom a outfile)) atoms))
  t)

(defun main ()
  (let ((atoms (load-params "opt.out")))
    (dump-params atoms)))
  
