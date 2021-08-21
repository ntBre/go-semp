;; maybe this should be a hash so it's easier to check membership?
(defparameter *derived_params*
  (list "PQN" "NValence" "DDN" "KON" "EISol"))

(defparameter *charge* 0)

(defparameter *spin* 1)

(defclass atomo ()
  ((labl :initarg :labl)
   (keys :initarg :keys :initform ())
   (vals :initarg :vals :initform ())))

(defun load-params (filename)
  (let (label param atom atoms (a -1))
    (with-open-file (in filename :direction :input)
      (loop for line = (read-line in nil nil)
	    while line do
	      (cond
		((string-equal line " ****") (setf label t))
		((and label (string-equal line " "))
		 (setf label nil
		       param nil))
		(label
		 (print line)
		 (setf param t
		       label nil))
		(param (print line))
		)))))
