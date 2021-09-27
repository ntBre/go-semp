(asdf:load-system :strings)
(asdf:load-system :stats)

(use-package :com.bwestbro.strings)
(use-package :com.bwestbro.stats)

(defparameter ht->cm 219474.5459784)

(defun load-params (filename)
  (with-open-file (in filename :direction :input)
    (loop for line = (read-line in nil nil)
	  while line
	  collect (parse-float line))))

(defparameter ai (load-params "/tmp/rel.dat"))
(defparameter se (load-params "/tmp/semp.dat"))

(defun square (x) (* x x))

(defun norm (a b)
  (* ht->cm (sqrt (reduce #'+ (mapcar #'square (mapcar #'- a b))))))

(defun average (a)
  (/ (reduce #'+ a) (length a)))

(defun rmsd (a b)
  (* ht->cm
     (sqrt (average (mapcar #'square (mapcar #'- a b))))))

(defun pah-scale (x)
  (cond
    ((< x 1111.1) (* x 0.956))
    ((< x 2500.0) (* x 0.952))
    (t            (* x 0.960))))

;; scaling b3lyp/4-31g harmonic freqs from gaussian
(mapcar #'pah-scale
	(list 
	 3324.1282
	 3285.9945
	 1624.5890
	 1257.6891
	 1074.8320
	 1011.0302
	 907.9355
	 858.3195
	 802.9297))

(defparameter b3lyp (list 3191.2 3154.6 1546.6 1197.3 1027.5 966.5 868.0 820.6 767.6))
(defparameter semp (list 3216.3 3187.4 1663.2 1171.3 1048.0 924.9 877.4 850.2 802.0))
(defparameter true (list 3142.6 3120.8 1600.9 1278.8 1061.5 970.2 886.4 878.6 787.4))
