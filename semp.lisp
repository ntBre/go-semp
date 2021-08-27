(in-package :cl-user)

(defpackage :com.bwestbro.semp
  (:use :cl :uiop))

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

(defmethod print-object ((atom atom-type) stream)
  (print-unreadable-object (atom stream :type t)
    (with-slots (label param-names param-values) atom
      (loop for name in param-names
	    for value in param-values
	    do (format t "~a: ~a=~a~%" label name value)))))

(defun load-params (filename)
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
				       (read-from-string line)) atoms)
		  (setf in-labl nil
			in-params t))
		 (in-params
		  (multiple-value-bind (names values) (parse-params line)
		    (push-names (car atoms) names)
		    (push-values (car atoms) (mapcar #'string->float values)))))))
    (nreverse atoms)))

(defun derived-param-p (param)
  (member (car (uiop:split-string param :separator '(#\=)))
	  (list "PQN" "NValence" "DDN" "KON" "EISol")
	  :test #'string-equal))

;; when dumping params dont print a new one if the name is the same as
;; the last, just put a comma
(defun parse-params (line)
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

(defparameter *mystr*
  "DCore=2,3,1.6101301456,3.2139710000,0.0000000000,0.0000000000,0.0000000000")
  ;; " PQN=2,2,0,0 NValence=4 F0ss=0.4900713060 F0sp=0.4236511293 F0pp=0.3644399818")


