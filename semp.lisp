(in-package :cl-user)

(defpackage :com.bwestbro.semp
  (:use :cl :cl-user :uiop :com.bwestbro.strings))

(in-package :com.bwestbro.semp)

;; TODO way to change these - probably read input file
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
  "Dump the parameters for ATOMS to FILENAME"
  (with-open-outfile filename
    (mapcar #'(lambda (a) (print-atom a outfile)) atoms))
  t)

(defmacro with-open-outfile (filename &body body)
    `(with-open-file (outfile ,filename :direction :output :if-exists :supersede)
       ,@body))

(defun write-com (filename geom &optional (param-file "params.dat"))
  (with-open-outfile filename
    ;; TODO should be able to adjust memory as well
    (format outfile "%mem=1000mb
%nprocs=1
#P PM6=(print,zero)

the title

~a ~a
~a

@~a
" *charge* *spin* geom param-file))
  t)

(defun write-pbs (filename)
  (with-open-outfile filename
    (format outfile "#!/bin/sh
#PBS -N sempirical
#PBS -S /bin/bash
#PBS -j oe
#PBS -m abe
#PBS -l mem=1gb
#PBS -l nodes=1:ppn=1

scrdir=/tmp/$USER.$PBS_JOBID

mkdir -p $scrdir
export GAUSS_SCRDIR=$scrdir
export OMP_NUM_THREADS=1

echo \"exec_host = $HOSTNAME\"

if [[ $HOSTNAME =~~ cn([0-9]{{3}}) ]];
then
  nodenum=${{BASH_REMATCH[1]}};
  nodenum=$((10#$nodenum));
  echo $nodenum

  if (( $nodenum <= 29 ))
  then
    echo \"Using AVX version\";
    export g16root=/usr/local/apps/gaussian/g16-c01-avx/
  elif (( $nodenum > 29 ))
  then
    echo \"Using AVX2 version\";
    export g16root=/usr/local/apps/gaussian/g16-c01-avx2/
  else
    echo \"Unexpected condition!\"
    exit 1;
  fi
else
  echo \"Not on a compute node!\"
  exit 1;
fi

cd $PBS_O_WORKDIR
. $g16root/g16/bsd/g16.profile
g16 {comfile} {outfile}

rm -r $scrdir

")) t)

;; TODO load energies from rel.dat
;; TODO load geometries from file07
(defun main ()
  (let ((atoms (load-params "opt.out")))
    (dump-params atoms)))
  
