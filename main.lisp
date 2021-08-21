;; maybe this should be a hash so it's easier to check membership?
(defparameter *derived_params*
  (list "PQN" "NValence" "DDN" "KON" "EISol"))

(defparameter *charge* 0)

(defparameter *spin* 1)

(defclass atom-type ()
  ((label :initarg :label
	  :accessor label)
   (params :initarg :params
	   :accessor params
	   :initform ())))

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
		 (in-params (print line)))))
    atoms))

(defun derived-param-p (param)
  (member (car (uiop:split-string param :separator '(#\=)))
	  *derived_params* :test #'string-equal))

(defun parse-params (line)
  (let ((fields (remove-if #'derived-param-p (fields line))))
    fields))

(defparameter *mystr*
  " PQN=2,2,0,0 NValence=4 F0ss=0.4900713060 F0sp=0.4236511293 F0pp=0.3644399818")

(defun new-char-array ()
  (make-array 1 :fill-pointer 0 :adjustable t :element-type 'character))

;; TODO split on something other than whitespace
;; nvm TODO use cl-ppcre, which has a split function
;; nvm nvm use uiop:split-string
;; nvm nvm nvm use this for splitting into fields, and then use uiop for splitting on commas
(defun fields (str)
  "split str on whitespace"
  (flet ((whitespace-p (c)
	   (member c (list #\Newline #\Space #\Tab))))
    (let ((len (1- (length str)))
	  (ret ())
	  (cur (new-char-array)))
      (loop for idx = 0 then (1+ idx)
	    for c across str
	    if (or (whitespace-p c) (= idx len))
	      do (when (> (length cur) 0)
		   (push (coerce cur 'string) ret)
		   (setf cur (new-char-array)))
	    else
	      do (vector-push-extend c cur)
	    end)
      (nreverse ret))))
