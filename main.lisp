;;; simulation's parameters
(defparameter *max-generation* 1000) ;number of generation
(defparameter *area-sq* 64) ;length of side
(defparameter *cell-number* (* *area-sq* *area-sq*))
(defparameter *dna-length* 512)
(defparameter *space-limit* 1) 
(defparameter *dna-limit* 16) 
(defparameter *mutation-rate* .001)
(defparameter *max-sp* 1)
(defparameter *density-rate* 1.00)
(defparameter *individual-number* (floor (* *cell-number* *density-rate*)))
(defparameter *dispersal-rate* 0.15)
(defparameter *dispersal-num* (floor (* *individual-number* *dispersal-rate*)))
(defconstant *death* -1)
(setf *random-state* (make-random-state t))

;;; utility macros & functions
(defmacro while (test &body body)
                 `(do ()
                    ((not ,test))
                    ,@body))

(defun flatten (lst)
  (labels ((mklist (x) (if (listp x) x (list x))))
    (mapcan #'(lambda (x) (if (atom x) (mklist x) (flatten x))) lst)))

(defun random-choice (lst)
  (nth (random (length lst)) lst))
 
;;; individual structure
(defstruct individual
  (sp-name 1)
  (dna (make-array `(2 ,*dna-length*)
                   :element-type 'bit
                   :initial-element 0)))

;;; field functions
(defun initial-field ()
  (let ((field (make-array *cell-number* 
                           :initial-element nil)))
    (dotimes (i *cell-number*)
          (setf (aref field i) (make-individual)))
    field))

(defun den-initial-field (num)
  (let ((field (initial-field)) (counter 0))
    (while (< counter (- *cell-number* num))
           (setf r (random *cell-number*))
           (when (aref field r)
             (setf (aref field r) nil)
             (incf counter)))
    field))

(defun copy-field (field)
  (let ((fld (initial-field)))
    (dotimes (i (length field))
           (if (null (aref field i))
             (setf (aref fld i) nil)
             (setf (aref fld i) (copy-structure (aref field i)))))
    fld))

;;; mutation functions

(defun dna-mutation (dna)
  (dotimes (i 2)
    (dotimes (j *dna-length*)
      (when (< (random 1.0) *mutation-rate*)
        (setf (aref dna i j) (- 1 (aref dna i j)))))))

(defun population-mutation (field)
  (let ((fld (copy-field field)))
    (dotimes (i (length fld))
      (when (aref fld i)
        (dna-mutation (individual-dna (aref fld i)))))
   fld))

;;; reproduct functions
(defun dna-difference (dna1 dna2)
  (let ((counter 0))
    (dotimes (i *dna-length*)
      (when (and (= (aref dna1 0 i) (aref dna1 1 i))
               (= (aref dna2 0 i) (aref dna2 1 i))
               (/= (aref dna1 0 i) (aref dna2 0 i)))
        (incf counter)))
    counter))

(defun area-points (center rudius)
  (let (lst (l (1+ (* 2 rudius))) 
        (min-p (- center (+ (* rudius *area-sq*) rudius))))
    (dotimes (i l)
      (dotimes (j l)
        (when (= (floor (+ (* i *area-sq*) j min-p) *area-sq*)
                 (floor (+ (* i *area-sq*) min-p rudius) *area-sq*))
          (setf lst (cons (+ (* i *area-sq*) j min-p) lst)))))
    (remove nil (mapcar #'(lambda (x)
                (if (or (< x 0) (<= *cell-number* x))
                  nil
                  x))
            lst))))

(defun fitness (diff)
  (- 1 (expt (/ diff *dna-limit*) 2)))

(defun could-mate? (ind1 ind2)
  (and (= (individual-sp-name ind1) (individual-sp-name ind2))
       (/= (individual-sp-name ind1) *death*)
       (not (eq ind1 ind2))
       (< (dna-difference (individual-dna ind1)
                                (individual-dna ind2))
          *dna-limit*)))

(defun could-mate-dm2? (ind1 ind2)
  (and (could-mate? ind1 ind2)
       (> (fitness (dna-difference (individual-dna ind1)
                                         (individual-dna ind2)))
          (random 1.0))))

(defun area-rep-list (field lst)
  (let (rp-lst)
    (labels ((sea (ls)
               (if (< (length ls) 2)
                 nil
                 (progn
                   (let ((ca (car ls)))
                   (mapcar #'(lambda (x)
                      (when (and 
                              (aref field ca)
                              (aref field x)
                              (could-mate? 
                                (aref field ca)
                                (aref field x))) 
                        (setf rp-lst (cons (list ca x) rp-lst))))
                           (cdr ls)))
                   (sea (cdr ls))))))
      (sea lst))
    rp-lst))

(defun ind-mate (ind1 ind2)
  (let ((ind (make-individual)))
    (setf (individual-sp-name ind) (individual-sp-name ind1))
    (dotimes (i *dna-length*)
      (setf (aref (individual-dna ind) 0 i)
            (aref (individual-dna ind1) (random 2) i)
            (aref (individual-dna ind) 1 i)
            (aref (individual-dna ind2) (random 2) i)))
    ind))

(defun next-gnr-ind (field point)
  (labels ((mate (lim)
       (let ((lst (area-rep-list field (area-points point lim))))
         (if lst
           (let ((pair (random-choice lst)))
             (ind-mate (aref field (car pair))
                        (aref field (cadr pair))))
           (if (< lim *area-sq*)
             (mate (1+ lim))
             (error "error: over field size"))))))
      (mate *space-limit*)))

(defun nil-cell-list (field)
  (let (lst)
    (dotimes (i *cell-number*)
           (when (null (aref field i))
             (setf lst (cons i lst))))
    lst))

(defun ind-cell-list (field)
  (let (lst)
    (dotimes (i *cell-number*)
      (when (individual-p (aref field i))
        (setf lst (cons i lst))))
    lst))

(defun population-reproduction (field)
  (let ((fld (copy-field field)))
    (dotimes (i *cell-number*)
      (when (aref fld i)
        (setf (aref fld i) (next-gnr-ind field i))))
    fld))

;;; speciation functions
(defun could-mate-list? (field lst num)
  (let (flag)
    (if (null lst)
      nil
      (mapcar #'(lambda (x)
         (when (and (aref field x)
                    (could-mate? (aref field x) (aref field num))) 
           (setf flag t)))
        lst))
    flag))
 
(defun com-speciation (field m-list)
  (let (r-list (o-list (copy-list m-list)))
    (mapcar #'(lambda (o)
      (let (con-flag-lst)
        (mapcar #'(lambda (r)
          (when (could-mate-list? field r o)
            (setf con-flag-lst (cons r con-flag-lst)))
          nil )
          r-list)
        (if (zerop (length con-flag-lst))
          (setf r-list (cons (list o) r-list))
          (if (= (length con-flag-lst) 1)
            (setf r-list (subst (cons o (car con-flag-lst))
                                (car con-flag-lst) 
                                r-list))
            (progn 
              (setf r-list (cons (flatten con-flag-lst) r-list))
              (setf r-list (set-difference r-list con-flag-lst)))))))
      o-list)
    (sort r-list #'(lambda (a b) (> (length a) (length b))))))

(defun population-speciation (field)
  (let ((fld (copy-field field)) (sp-list (sp-list field)))
    (mapcar #'(lambda (s)
      (let ((af-s-list (com-speciation field (sp-member fld s))))
        (dotimes (i (length af-s-list))
          (when (/= i 0)
            (mapcar #'(lambda (x)
                        (setf (individual-sp-name (aref fld x))
                              *max-sp*))
                    (nth i af-s-list))
            (setf *max-sp* (1+ *max-sp*))))))
            sp-list)
    fld))

;;; dispersal functions
(defun population-dispersal (field)
  (let ((fld (copy-field field)) (num *dispersal-num*))
    (while (> num 0)
           (let* ((move-ind (random-choice (ind-cell-list fld)))
                  (move-to (random-choice (remove move-ind (area-points move-ind 1)))))
             (if (individual-p (aref fld move-to))
               (let ((c (copy-individual (aref fld move-to))))
                 (setf (aref fld move-to) (copy-individual (aref fld move-ind)))
                 (setf (aref fld move-ind) (copy-individual c))
                 (setf num (- num 2)))
               (progn
                 (setf (aref fld move-to) (copy-individual (aref fld move-ind)))
                 (setf (aref fld move-ind) nil)
                 (decf num)))))
    fld))

;;; profile functions

(defun some-sp-num (field sp)
  (let ((counter 0))
    (dotimes (i *cell-number*)
      (when (and (aref field i)
                 (= (individual-sp-name (aref field i)) sp)) 
        (incf counter)))
    counter))
 
(defun name-of (ind)
  (char "abcdefghijklmnopqrstuvwxyABCDEFGHIJKLMNOPQRSTUVWXYZ" ind))

(defun print-field (field)
  (dotimes (i *cell-number*)
         (progn
           (if (null (aref field i))
             (format t "_ " )
             (format t "~A " 
                     (name-of (mod (individual-sp-name (aref field i)) 50))))
           (when (= (mod i *area-sq*) (- *area-sq* 1))
             (format t "~%")))))

(defun sp-list (field)
  (let (lst)
    (dotimes (i (length field))
      (unless (or (null (aref field i)) 
                  (member (individual-sp-name (aref field i)) lst)) 
        (setf lst (cons (individual-sp-name (aref field i)) lst))))
    lst))

(defun sp-rich (field)
  (length (sp-list field)))
  
(defun ind-num (field)
  (let (r-list (s-list (sp-list field)))
    (mapcar #'(lambda (x)
                (setf r-list (cons (some-sp-num field x) r-list)))
            s-list)
    (sort r-list #'>)))

(defun make-field-sta-file (field generation)
  (setf path (make-pathname :name (concatenate 'string "field" (princ-to-string generation))))
  (setf str1 (open path :direction :output :if-exists :supersede))
  (dotimes (i *cell-number*)
    (progn 
      (when (= (mod i *area-sq*) 0)
        (format str1 "~%"))
      (if (null (aref field i))
        (format str1 "nil ")
        (format str1 "~A " (individual-sp-name (aref field i))))))
  (close str1)
  nil)

(defun nil-cell (field point)
  (let ((p-lst (area-points point 1)) r-lst)
    (mapcar #'(lambda (x)
                (when (null (aref field x))
                  (setf r-lst (cons x r-lst))))
            p-lst)
    r-lst))

(defun ind-cells (field point)
  (let ((p-lst (area-points point 1)) r-lst)
    (mapcar #'(lambda (x)
                (when (individual-p (aref field x))
                  (setf r-lst (cons x r-lst))))
            p-lst)
    r-lst))

(defun sp-member (field sp)
  (let (lst)
    (dotimes (i *cell-number*)
      (if (and 
            (aref field i)
            (= (individual-sp-name (aref field i)) sp))
        (setf lst (cons i lst))))
    lst))

(defun ind-hetero-rate (ind)
  (let ((hetero-sum 0))
    (dotimes (i *dna-length*)
       (when (/= (aref (individual-dna ind) 0 i)
              (aref (individual-dna ind) 1 i))
         (setf hetero-sum (1+ hetero-sum))))
    (/ (float hetero-sum) *dna-length*)))

(defun com-hetero-rate (field)
  (let ((hetero-sum 0))
    (dotimes (i *cell-number*)
           (when (aref field i)
             (setf hetero-sum (+ (ind-hetero-rate (aref field i)) 
                          hetero-sum))))
    (/ hetero-sum *cell-number*)))
 
;;; main function

(defun main-loop ()
  (setf path (make-pathname :name (concatenate 'string  "sp-rich" (princ-to-string (get-universal-time)))))
  (setf str (open path :direction :output
                    :if-exists :supersede))
  (let ((field (den-initial-field *individual-number*)))
    (dotimes (generation *max-generation*)
      (progn
        (setf field (population-reproduction field))
        (setf field (population-dispersal field))
        (population-mutation field)
        (setf field (population-speciation field))
        (print-field field)
        (format t "generation = ~A   hetero-r = ~A  sprich = ~A~%"
                generation (com-hetero-rate field) (sp-rich field))
        (format str "~a~%" (sp-rich field))
        (make-field-sta-file field generation)))
    (mapcar #'print (ind-num field))
    (format str "~A~%"(ind-num field)))
  (format str "area-sq = ~A ,individual-num = ~A ,dna-len = ~A , space-lim = ~A , dna-lim = ~A , mutation-rate = ~A "
          *area-sq* *cell-number* *dna-length* *space-limit*
          *dna-limit* *mutation-rate*)
  (close str)
  nil)

(main-loop)
