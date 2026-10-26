;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; LISP VERSION OF L1N CONSENSUS ;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; IMPLEMENTED USING ALLEGRO CL v11.0 FREE EDITION ;;;;;
;;;;;;;;;; ALL FUNCTIONS MUST BE COMPILED FIRST ;;;;;;;;;;;;;;;;
;;;;;;;;;; THEN USE COMMAND LINES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; The input networks are given in a format similar to extended Newick.
; The 5 networks of Figure 2 example in the paper are encoded as:


(defparameter *INPUT-NETS-1*
  '(
    
    ;Net1
    ((A ((M1 (M2) X1) n3 (X1 M3) n4) n2 ((B (C) X2) n6 (E  (X2 (D1 D2) n9) n8) n7) n5) n1)
    
    ;NET2
    ((A ((M1 M2) n3 M3) n2 ((B (C) X2) n6 (E  (X2 (D1 D2) n9) n8) n7) n5) n1)
    
    ;NET3
    ((A (E ((M1 (C) X1) n4 (M3 X1) n5) n3) n2 (M2 (B (D1 D2) n8) n7) n6) n1)
    
    ;NET4
    ((A ((D1 D2) n5 (C (B E) n8) n7) n6 ((M1 (M2) X1) n3 (X1 M3) n4) n2) n1)
    
    ;NET5
    ((A (C ((D1 D2) n5 (B E) n8) n7) n6 ((M1 (M2) X1) n3 (X1 M3) n4) n2) n1)     
    
    ))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; In this second toy example, Net2 is removed and we have Net3 twice
; See supplementary material ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *INPUT-NETS-2*
  '(
    
    ;Net1
    ((A ((M1 (M2) X1) n3 (X1 M3) n4) n2 ((B (C) X2) n6 (E  (X2 (D1 D2) n9) n8) n7) n5) n1)
    
    ;NET3
    ((A (E ((M1 (C) X1) n4 (M3 X1) n5) n3) n2 (M2 (B (D1 D2) n8) n7) n6) n1)
    
    ;NET3
    ((A (E ((M1 (C) X1) n4 (M3 X1) n5) n3) n2 (M2 (B (D1 D2) n8) n7) n6) n1)
    
    ;NET4
    ((A ((D1 D2) n5 (C (B E) n8) n7) n6 ((M1 (M2) X1) n3 (X1 M3) n4) n2) n1)
    
    ;NET5
    ((A (C ((D1 D2) n5 (B E) n8) n7) n6 ((M1 (M2) X1) n3 (X1 M3) n4) n2) n1)     
    
    ))


; Additional information must be provided:


(defparameter *outgroup* nil)                    ;no outgroup taxon, the input networks are unrooted
(defparameter *outgroup-code* nil)               ;outgroup numerical code, when there is an outgroup taxon


(defparameter *names*                            ;numerical encoding of the taxa
'((1 A) 
  (2 B) 
  (3 C) 
  (4 D1) 
  (5 D2) 
  (6 E) 
  (7 M1) 
  (8 M2) 
  (9 M3)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LAUNCHING L1N CONSENSUS, RESULT WITH FIRST TOY EXAMPLE OF FIG. 2 ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
; 
; CG-USER(11): (ini_part_nets)
; NIL
; CG-USER(12): (extended-network *INPUT-NETS-1*)
;
; ************************************* 
; ENCODING_OF_TAXA 
; ************************************* 
; (1 A) 
; (2 B) 
; (3 C) 
; (4 D1) 
; (5 D2) 
; (6 E) 
; (7 M1) 
; (8 M2) 
; (9 M3) 
; ************************************* 
; * 
; *-----BLOBS----- 
; * 
; ***** 
; (BLOB-SUPPORT 3) 
; (BEST-CYCLE SUPPORT 3 #(2 2 2 2 2 2 1 0 3))         ;vectorial encoding of circular orders; entry 1 is A, entry 2 is B, ... entry 9 is M3
; (HYBRID-SUPPORT 3)                                  ;numbers indicate taxon block; 2 groups all taxa but M1, M2, and M3
; *****                                               ;0 is the hybrid M2 (entry 8); M1 (entry 7, block 1) and M3 (block 3) are the parental taxa
; ***** 
; (BLOB-SUPPORT 2) 
; (BEST-CYCLE SUPPORT 2 #(2 1 0 4 4 3 2 2 2))         ;the hybrid is C (0 in entry 3); the parental taxa are B and (D1, D2)
; (HYBRID-SUPPORT 3)                                  ;the supports are given as counts, e.g. 3 over 5 input networks
; ***** 
; * 
; *--NON-REDUNDANT-BIPARTITIONS-- 
; * 
; *                                                   ;all splits are redundant in this example
; *--IRRESOLUTION-- 
; * 
; (IRRESOLUTION = 0)                                  ;the network is fully resolved, divergence nodes (if any) have degree 3
; *                                                   ;a value of 1 means 1 divergence node with degree 4; 2 means 2 such nodes or 1 of degree 5, etc.
; *--EASYNESS-- 
; * 
; (EASYNESS 1)                                        ;counts the number of occurences of the most frequent input network
; *                                                   ;the higher the easier; with strict majority L1N is ensured to output that network
; (NET-FREQUENCY 1)                                   ;counts the number of occurences of the output network among the inputs
; **********                                          in this example, the ouput is identical to Net1 (see paper) 
; **********
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; LAUNCHING L1N CONSENSUS, RESULT WITH SECOND TOY EXAMPLE - SUP MAT ;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;
; CG-USER(13): (ini_part_nets)
; 
; CG-USER(14): (extended-network *INPUT-NETS-2*)
; 
; ************************************* 
; ENCODING_OF_TAXA 
; ************************************* 
; (1 A) 
; (2 B) 
; (3 C) 
; (4 D1) 
; (5 D2) 
; (6 E) 
; (7 M1) 
; (8 M2) 
; (9 M3) 
; ************************************* 
; * 
; *-----BLOBS----- 
; * 
; ***** 
; (BLOB-SUPPORT 3) 
; (BEST-CYCLE SUPPORT 3 #(2 2 2 2 2 2 1 0 3)) 
; (HYBRID-SUPPORT 3) 
; ***** 
; * 
; *--NON-REDUNDANT-BIPARTITIONS-- 
; * 
; (4 #(0 0 0 1 1 0 0 0 0))                                ;4 non-redundant splits in this second example 
; (2 #(0 1 0 1 1 0 0 0 0)) 
; (2 #(0 1 1 1 1 1 0 0 0)) 
; (1 #(0 1 0 1 1 1 0 0 0)) 
; * 
; *--IRRESOLUTION-- 
; * 
; (IRRESOLUTION = 0)                                      ;fully resolved
; * 
; *--EASYNESS-- 
; * 
; (EASYNESS 2)                                            ;relatively easy (2 / 5)
; * 
; (NET-FREQUENCY 0) 
;                                                        ;however, the output network is different from all input network
; **********
; **********


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;; END OF EXAMPLES ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; UTILITIES ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defmacro nconc-p (l y)
  `(if ,l (nconc ,l ,y) (setf ,l ,y)))

(defun pp (l)
  (dolist (x l)(print x)))

(defun copy-v (v)
  (let* ((li (array-dimension v 0))
         (cv (make-array li)))
(dotimes (i li cv)
  (setf (aref cv i)(aref v i)))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; PARSING L1 NETWORKS IN EXTENDED-NEWICK-LIKE FORMAT ;;
;;;;;;;;;; INTERNAL NODES ARE LABELED ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; NETWORKS BECOME GRAPHS WITH VERTICES AND EDGES ;;;;;;
;;;;;;;;;; NODES AND TAXA ARE ATOMS WITH PROPERTIES ;;;;;;;;;;;;
;;;;;;;;;; E.G. NEIGHBORS, TYPE, ETC ;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defparameter *ltaxa* nil)
(defparameter *lnodes* nil)
(defparameter *lhybrids* nil)
(defparameter *lcycles* nil)
(defparameter *parts* nil)
(defparameter *tax* 0)

(defun ini_part_nets ()                          ;initialization of the data structure when *names* is given
(setf *parts* nil)                               ;to be launched before extended-network
(setf *tax* (length *NAMES*))
(dolist (x *NAMES*)
  (let ((i (car x))
        (a (cadr x)))
    (setf (get a 'num) (- i 1)))))

(defun encodetax ()                              ;used when taxa are not numbered and *names* is nil
  (let ((l (length *ltaxa*)))
    (dotimes (i l)
      (setf (get (nth i *ltaxa*) 'num) i)
      (when (equal *outgroup* (nth i *ltaxa*))
        (setf *outgroup-code* i)))))

(defun printparsing ()
  (print '*************************************)
  (print 'ENCODING_OF_TAXA)
  (print '*************************************)
  (dotimes (i *tax*)
    (print (list (+ 1 i)
                 (dolist (x *ltaxa*)
                   (when (= i (get x 'num))
                     (return x))))))
  (when *outgroup*
    (print '*************************************)
    (print (list 'outgroup (+ 1 *outgroup-code*) *outgroup*)))
  (print '*************************************))

(defun parsenet (l)
  (setf *ltaxa* nil)
  (setf *lnodes* nil)
  (parsaux l)
  (setf *lhybrids* (spe-intersect *ltaxa* *lnodes*))
  (dolist (h *lhybrids*)
    (setf *ltaxa* (remove h *ltaxa*))
    (setf *lnodes* (remove h *lnodes*)))
  (setf *tax* (length *ltaxa*))
  (unless *names* 
    (encodetax)
    (setf *names* T))
  (ascendance l)
  (removesingle)
  (extractcycles)
  (orderedpartitions))

(defun parsaux (l)
  (cond ((null l) T)
        ((atom l)
         (setf *ltaxa* (cons l *ltaxa*)))
        ((atom (car l))
         (parsaux (car l))
         (parsaux (cdr l)))
        (T
         (setf *lnodes* (cons (cadr l) *lnodes*))
         (parsaux (car l))
         (parsaux (cddr l)))))

(defun spe-intersect (l1 l2)                   ;We assume all nodes and types are atoms (not numbers)
  (let ((res nil))
    (dolist (x l1)
      (setf (get x 'inter) T)
      (setf (get x 'type) 'taxon))
    (dolist (x l2)
      (setf (get x 'type) 'node)
      (when (get x 'inter)
        (setf res (cons x res))
        (setf (get x 'type) 'hybrid)))
    (dolist (x l1)
      (setf (get x 'inter) nil))
    res))

(defun ascendance (l)
  (let ((all (append *lnodes* *ltaxa* *lhybrids*)))
    (dolist (x all)
      (setf (get x 'neighbors) nil)
      (setf (get x 'minor) nil)
      (setf (get x 'major) nil))
    (ascaux (car l)(cadr l))))

(defun ascaux (l anc)
  (if (null l) T
    (let ((x (car l)))
      (cond ((not (atom x))
             (let ((n (cadr l)))
               (addneigh n anc)
               (when (equal (get n 'type) 'hybrid)
                 (setf (get n 'major) anc))
               (ascaux (cddr l) anc)
               (ascaux x n)))
            (T
             (addneigh x anc)
             (when (equal (get x 'type) 'hybrid)
               (setf (get x 'minor) anc))
             (ascaux (cdr l) anc))))))

(defun addneigh (n1 n2)
  (setf (get n1 'neighbors)
   (cons n2 (get n1 'neighbors)))
  (setf (get n2 'neighbors)
    (cons n1 (get n2 'neighbors))))

(defun removesingle ()
  (dolist (x *lnodes*)
    (let ((ln (get x 'neighbors)))
      (when (= 2 (length ln))
        (let ((a (car ln))
              (b (cadr ln)))
          (setf (get a 'neighbors)
            (remove x (get a 'neighbors)))
          (setf (get b 'neighbors)
            (remove x (get b 'neighbors)))
          (addneigh a b)
          (setf *lnodes* (remove x *lnodes*))
          (return T))))))

(defun extractcycles ()
  (setf *lcycles* nil)
  (dolist (h *lhybrids*)
    (setf *lcycles* 
      (cons (cycle h) *lcycles*))))

(defun cycle (h)
  (let ((s (get h 'major))
        (e (get h 'minor)))
    (cons h (cycle-aux s h e))))

(defun cycle-aux (point origin destination)
  (cond ((equal point destination)(list destination))
        ((or (equal 'taxon (get point 'type))
             (equal 'hybrid (get point 'type)))                   ;Input netwroks must be L1
         nil)
        (T (let* ((l (othersn point origin))
                  (x (car l))
                  (cx (cycle-aux x point destination))
                  (y (cadr l))
                  (cy (cycle-aux y point destination)))
             (cond (cx (cons point cx))
                   (cy (cons point cy))
                   (T nil))))))    

(defun othersn (point origin)
  (let ((res nil)
(l (get point 'neighbors)))
(dolist (x l res)
  (unless (equal x origin)
    (setf res (cons x res))))))

(defun orderedpartitions ()
  (setf *parts* nil)
  (dolist (cyc *lcycles*)
    (setf *parts* (cons (cycpartition cyc) *parts*))))

(defparameter *vec* nil)

(defun cycpartition (cyc)
  (setf *vec* (make-array (length *ltaxa*) :initial-element nil))
  (let ((ll (length cyc)))
    (dotimes (i ll)
      (cond ((= i 0)
             (encodeclass i (nth i cyc)(nth (+ 1 i) cyc)(nth (- ll 1) cyc)))
            ((= i (- ll 1))
             (encodeclass i (nth i cyc)(nth 0 cyc)(nth (- i 1) cyc)))
            (T
             (encodeclass i (nth i cyc)(nth (+ 1 i) cyc)(nth (- i 1) cyc)))))
    *vec*))

(defun encodeclass (i origin a b)
  (let ((neigh (get origin 'neighbors))	(point nil))
    (setf point
      (dolist (x neigh)
        (unless (or (equal x a)(equal x b))
          (return x))))
    (enclassaux i point origin)))

(defun enclassaux (i point origin)
  (if (equal 'taxon (get point 'type))
      (setf (aref *vec* (get point 'num)) i)
    (let ((l (othersn point origin))
          (ll nil))
      (dolist (x l)
        (unless (or (equal point (get x 'minor))
                    (equal x (get point 'minor)))
          (setf ll (cons x ll))))
      (dolist (x ll)
        (enclassaux i x point)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; EXTRACT NON-TRIVIAL NON-REDUNDANT BIPARTITIONS ;;;;;;
;;;;;;;;;; A REDUNDANT BIPARTITION IS INDUCED BY A BLOB ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defparameter *lbipartitions* nil)

(defun makeblobs ()
  (dolist (x *ltaxa*)
(setf (get x 'blob-neighbors)
 (copy-list (get x 'neighbors))))
(dolist (x *lhybrids*)
  (setf (get x 'blob-neighbors)
    (copy-list (get x 'neighbors))))
(dolist (x *lnodes*)
  (setf (get x 'blob-neighbors)
    (copy-list (get x 'neighbors))))
(dolist (c *lcycles*)
  (blobaux c)))

(defun blobaux (cycle)
  (let ((h (car cycle)))
    (dolist (n (cdr cycle))
      (collapse-blob n h cycle))))

(defun collapse-blob (n h cycle)
  (dolist (x (get n 'blob-neighbors))
    (unless (or (member x (get h 'blob-neighbors))
             (member x cycle))
      (add-edge-blob h x)
      (remove-edge-blob n x)))
  (setf (get h 'blob-neighbors)
    (remove n (get h 'blob-neighbors))))

(defun add-edge-blob (h l)
  (setf (get h 'blob-neighbors)
    (cons l (get h 'blob-neighbors)))
  (setf (get l 'blob-neighbors)
    (cons h (get l 'blob-neighbors))))

(defun remove-edge-blob (h l)
  (setf (get h 'blob-neighbors)
    (remove l (get h 'blob-neighbors)))
  (setf (get l 'blob-neighbors)
    (remove h (get l 'blob-neighbors))))

(defun othersn-blob (point origin)
  (let ((res nil)
        (l (get point 'blob-neighbors)))
    (dolist (x l res)
      (unless (equal x origin)
        (setf res (cons x res))))))

(defun printn (l)
  (dolist (x l)
    (print x)
    (xxx x)))

(defun xxx (x)
  (print (get x 'type))
  (print (get x 'neighbors))
  (print (get x 'blob-neighbors)))

(defun extractbipartitions ()
  (makeblobs)
  (setf *lbipartitions* nil)
  (let* ((o (car *ltaxa*))
(x (car (get o 'blob-neighbors)))
(ln (get x 'blob-neighbors))
(ind (< (length ln) 4)))                     ;ind cheks whether x is a blob center, 
(dolist (n ln)                                    ;in this case the bipartitions x-n are redundant
  (unless (equal n o)
    (extbipaux n x ind)))))

(defun extbipaux (point origin indorigin)
  (if (equal 'taxon (get point 'type))
      (makebiptax point)
    (let* ((ln (get point 'blob-neighbors))
(ind (< (length ln) 4))                          ;ind checks whether point is a blob center
(vec (make-array *tax* :initial-element 0)))     ;in this case the bipartitions point-n are redundant
(dolist (n ln vec)
  (unless (equal n origin)
    (setf vec (orvec vec (extbipaux n point ind)))))
(when (and indorigin ind)                                               ;neither the origin nor the point are blobs
  (setf *lbipartitions* (cons (reorder vec) *lbipartitions*)))          ;we add the bipartition origin-point to the list
vec)))                                                                  

(defun makebiptax (n)
  (let ((v (make-array *tax* :initial-element 0)))
    (setf (aref v (get n 'num)) 1)
    v))

(defun orvec (v1 v2)
  (let ((v (make-array *tax* :initial-element 0)))
    (dotimes (i *tax* v)
      (when (or (= 1 (aref v1 i))
                (= 1 (aref v2 i)))
        (setf (aref v i) 1)))))

(defun reorder (v)
  (let* ((ll (length v))
         (nbc (numbclass v))
         (corsp (make-array nbc :initial-element nil))
         (ptr 0)
         (res (make-array ll)))
    (dotimes (i ll)
(when (null (aref corsp (aref v i)))
  (setf (aref corsp (aref v i)) ptr)
  (setf ptr (+ 1 ptr))))
(dotimes (i ll res)
  (setf (aref res i)(aref corsp (aref v i))))))

(defun numbclass (vec)
  (let ((res 0)
        (lv (length vec)))
    (dotimes (i lv (+ 1 res))
      (setf res (max res (aref vec i))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; EXTRACT AND COUNT NON-REDUNDANT BIPARTITIONS ;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *all_bipartitions* nil)

(defun countbipartitions (lnets)
  (setf *all_bipartitions* nil)
  (let ((res nil))
    (dolist (net lnets (sort res '>car))                          ;(sort (reverse res) '>car))
      (parsenet net)
      (extractbipartitions)
      (nconc-p *all_bipartitions* (list *lbipartitions*))
      (dolist (bip *lbipartitions*)
        (let ((memb (member (list 'bid bip) res :test 'cadr=)))
          (if memb
              (setf (caar memb) (+ 1 (caar memb)))
            (nconc-p res (list (list 1 bip)))))))))

(defun >car (x y)(> (car x)(car y)))

(defun cadr= (x y)
  (let ((v (cadr x))(w (cadr y)))
    (dotimes (i *tax* T)
      (unless (= (aref v i)(aref w i))
        (return nil)))))

(defun all-bipcompatible (bip lbips)
  (dolist (b lbips T)
    (unless (bipcompatible bip (cadr b))
      (return nil))))

(defun bipcompatible (v1 v2)
  (let ((ind01 nil)(ind10 nil)(ind00 nil)(ind11 nil))
    (dotimes (i *tax* T)
(cond ((and (= 0 (aref v1 i))
            (= 1 (aref v2 i)))
       (setf ind01 T))
      ((and (= 1 (aref v1 i))
        (= 0 (aref v2 i)))
(setf ind10 T))
((and (= 0 (aref v1 i))
      (= 0 (aref v2 i)))
 (setf ind00 T))
((and (= 1 (aref v1 i))
      (= 1 (aref v2 i)))
 (setf ind11 T)))
(when (and ind01 ind10 ind00 ind11) (return nil)))))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; EXTRACT, COUNT AND COLLECT COMPATIBLE BLOBS ;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; A BLOB HAS DEGREE >= 4 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; A BLOB IS OBTAINED BY COLLAPSING A CYCLE ;;;;;;;;;;;;
;;;;;;;;;; WE FIRST EXTRACT THE CYCLES ;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; WHICH ARE THEN TURNED INTO BLOBS ;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defparameter *all_cycles* nil)

(defun countcycles (lnets)
  (setf *all_cycles* nil)
  (let ((res nil))
    (dolist (net lnets (sort res '>car))                                        
      (parsenet net)
      (nconc-p *all_cycles* (list *parts*))
      (dolist (cyc *parts*)
        (let ((memb (member (list 'bid cyc) res :test 'cadr=)))
          (if memb
              (setf (caar memb) (+ 1 (caar memb)))
            (nconc-p res (list (list 1 cyc)))))))))

(defun countblobs (lnets)
  (let ((res nil)
        (lcycles (countcycles lnets)))
    (dolist (cyc lcycles)
      (let ((memb (member (list cyc) res :test 'blob=)))
        (if memb
            (nconc-p (car memb)(list cyc))
          (nconc-p res (list (list cyc))))))
    (setf res (makesum res))
    (setf res (sort res '>cadr))                                                            
    (reuni-cycles res)))

(defun blob= (lx ly)
  (let* ((x (cadar lx))
         (rx (reorder x))
         (y (cadar ly))
         (ry (reorder y)))
    (dotimes (i *tax* T)
      (unless (= (aref rx i)(aref ry i))
        (return nil)))))

(defun makesum (llcycles)
  (let ((res nil))
    (dolist (lcyc llcycles res)
      (let ((sum 0))
        (dolist (cyc lcyc)
          (setf sum (+ sum (car cyc))))
        (setf res (cons (cons sum lcyc) res))))))

(defun >cadr (lx ly)
  (> (car lx)(car ly)))

(defun reuni-cycles (lblobs)
  (let ((res nil))
    (dolist (blob lblobs res)
      (nconc-p res (list (cons (car blob)(reuniaux (cdr blob))))))))

(defun reuniaux (lcycles)
  (let ((res nil))
    (do ((lcyc lcycles (cdr lcyc)))
        ((endp lcyc) res)
      (unless (= -1 (caar lcyc))
        (let* ((cyc (car lcyc))
(memb (member cyc (cdr lcyc) :test 'rooted-cyc=)))
(cond (memb
       (nconc-p res (list (cons (+ (car cyc)(caar memb))(cdr cyc))))
       (setf (caar memb) -1))
      (T (nconc-p res (list cyc)))))))))


(defun allcompatible (blob lblobs)
  (dolist (b lblobs T)
    (unless (compatible blob (cadadr b))
      (return nil))))

(defun compatible (b1 b2)                               ; b1 = bip or blob, b2 = blob, NIL when b1 is a redundant bip
  (let* ((nc1 (numbclass b1))                                     
         (nc2 (numbclass b2))
         (tab (make-array (list (+ 1 nC1)(+ 1 nC2)) :initial-element 0)))
    (dotimes (x *tax*)
      (when (= 0 (aref tab (aref b1 x)(aref b2 x)))
        (setf (aref tab (aref b1 x)(aref b2 x)) 1)
        (setf (aref tab (aref b1 x) nc2)
          (+ 1 (aref tab (aref b1 x) nc2)))
        (setf (aref tab nc1 (aref b2 x))
          (+ 1 (aref tab nc1 (aref b2 x))))))
    (check-table tab nc1 nc2)))

(defun check-table (tab nc1 nc2)
  (let ((maxi (lemaxi tab nc1 nc2))
        (maxj (lemaxj tab nc1 nc2)))
    (cond ((or (null maxi)(null maxj)) nil)
          ((and (= (aref tab maxi nc2) nC2)
                (= (aref tab nc1 maxj) nC1))
           T)
          ((and (= (aref tab maxi nc2) (- nC2 1))
                (= (aref tab nc1 maxj) (- nC1 1))
                (= (aref tab maxi maxj) 0)
                (> nc1 2))
T)
          (T nil))))

(defun lemaxi (tab nc1 nc2)
  (let ((ind nil))
    (dotimes (i nC1 ind)
      (let ((x (aref tab i nc2)))
        (cond ((= x 1) nil)
              ((and (null ind)
                    (>= x (- nc2 1)))
               (setf ind i))
              (T (return nil)))))))

(defun lemaxj (tab nc1 nc2)
  (let ((ind nil))
    (dotimes (j nC2 ind)
      (let ((x (aref tab nc1 j)))
        (cond ((= x 1) nil)
              ((and (null ind)
                    (>= x (- nc1 1)))
               (setf ind j))
              (T (return nil)))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; EXTENDED MAJORITY CONSENSUS NETWORK ;;;;;;;;;;;;;;;;;
;;;;;;;;;; WE USE A SIMPLE GREEDY APPROACH ;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; WE FIRST INFER THE MAJORITY-RULE TOB ;;;;;;;;;;;;;;;;
;;;;;;;;;; AND THEN GREEDILY ADD BLOBS/NON-REDUNDANT SPLITS ;;;;
;;;;;;;;;; UNTIL NO MORE CAN BE ADDED ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;; WE CAN ALSO ENFORCE THE ADDITION OF BLOBS ;;;;;;;;;;;
;;;;;;;;;; AND STOP THE ADDITION OF BLOBS AND SPLITS ;;;;;;;;;;;
;;;;;;;;;; FINALLY WE RESOLVE THE BLOBS WITH THE MOST FREQUENT ;
;;;;;;;;;; CIRCULAR ORDERS AND HYBDRID CLADES ;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defparameter *stop-addition* 0)            ;default = 0, we then obtain an extended consensus
;51 with 100 input networks gives the majority-rule consensus


(defparameter *trsh-blobs* 3)               ;default = strict majority; 51 with 100 input networkss; 3 in toy example of Figure 2
;lower values enforce the addition of blobs with frequency >= *trsh-blobs*

(defparameter *allblobs* nil)

(defparameter *allbips* nil)

(defun extended-network (lnets)
  (let* ((maj (/ (length lnets) 2.0))
         (allblobs (countblobs lnets))         
         (extblobs (keepmajblobs allblobs maj))
         (otherblobs (otherblobs allblobs maj))         
         (lbips (countbipartitions lnets))
         (extbips (keepmajbips lbips maj))
         (otherbips (otherbips lbips maj))         
         (lcycs nil)
         (maxtime (+ (length otherblobs)(length otherbips))))
    (setf *allblobs* allblobs)(setf *allbips* lbips)
    (unless (> *stop-addition* maj)
      (dotimes (i maxtime)
        (when (= 0 (resolutionlevel extblobs extbips)) (return T))
        (let* ((maxfreq (themax otherblobs otherbips))
               (blob (if (equal maxfreq 'blob)(car otherblobs) nil))
               (bip (if (equal maxfreq 'bip)(car otherbips) nil)))
          (cond ((and blob 
                      (allcompatible (cadadr blob) extblobs)
                      (allcompat-blobbip (cadadr blob) extbips))
                 (nconc-p extblobs (list blob))
                 (setf otherblobs (cdr otherblobs))
                 (setf extbips (clean-redundancy extbips (cadadr blob))))
                (blob (setf otherblobs (cdr otherblobs)))
                ((and bip 
                      (all-bipcompatible (cadr bip) extbips)
                      (allcompatible (cadr bip) extblobs))
                 (nconc-p extbips (list bip))
                 (setf otherbips (cdr otherbips)))
                (bip (setf otherbips (cdr otherbips)))))))
    
    (let ((tabres (make-array (list (length extblobs) 10)))      ;Vec Class Sup-blob Sup-cy 1stH Sup>4 Sup=4 2ndH  Sup>4   Sup=4
          ; 0    1      2       3      4     5     6    7     8       9
          (root nil))                                            ;After rooting            bestHyb           SupMaj SupMin SupBoth
      
      (dotimes (i (length extblobs))
        (let* ((b (nth i extblobs))
               (cycle (best-cycle b))
               (vec (cadr cycle))
               (class (numbclass vec))
               (supblob (car b))
               (supcycle (car cycle))
               (2-besthyb (2-best-hybrid vec allblobs)))                                 
          (setf (aref tabres i 0) vec)
          (setf (aref tabres i 1) class)
          (setf (aref tabres i 2) supblob)
          (setf (aref tabres i 3) supcycle)
          (dotimes (j 6)
            (setf (aref tabres i (+ 4 j)) (aref 2-besthyb j)))))
      
      (setf root (if *outgroup* *outgroup-code* (best-rooting tabres)))       
      (clean-tabres root tabres allblobs)                        ;select best compatible hybrid among 1st & 2nd, reroot the cycle, compute supp origins
      
      (printparsing)
      (print '*)
      (print '*-----BLOBS-----)
      (print '*)
      (dotimes (i (length extblobs))
        (nconc-p lcycs (list (aref tabres i 0)))
        (print '*****)
        (print (list 'blob-support (aref tabres i 2)))
        (print (list 'best-cycle 'support (aref tabres i 3) (aref tabres i 0)))
        (print (list 'hybrid-support (+ (aref tabres i 5)(aref tabres i 6))))
        (print '*****)))
    (print '*)
    (print '*--NON-REDUNDANT-BIPARTITIONS--)
    (print '*)
    (pp extbips)
    (print '*)
    (print '*--IRRESOLUTION--)   
    (print '*)
    (print (list 'irresolution '= (resolutionlevel extblobs extbips)))
    (print '*)
    (print '*--EASYNESS--)   
    (print '*)
    (let ((d (easyness)))
      (print (list 'EASYNESS d)))
    (print '*)
    (let ((fr (net-freq lcycs (lcadr extbips))))
      (print (list 'NET-FREQUENCY fr)))
    (print '**********)))

(defun lcadr (l)
  (let ((res nil))
    (dolist (x l res)
      (nconc-p res (list (cadr x))))))

(defun keepmajbips (lbips maj)
  (let ((res nil))
    (dolist (bip lbips res)
      (if (<= (car bip) maj)
          (return res)
        (nconc-p res (list bip))))))

(defun otherbips (lbips maj)
  (let ((res nil))
    (dolist (bip lbips res)
      (when (<= (car bip) maj)
        (nconc-p res (list bip))))))

(defun keepmajblobs (l maj)
  (let ((res nil))
    (dolist (x l res)
      (if (<= (car x) maj)(return res)
        (nconc-p res (list x))))))

(defun otherblobs (lblobs maj)
  (let ((res nil))
    (dolist (x lblobs res)
      (when (<= (car x) maj)        
        (nconc-p res (list x))))))

(defun resolutionlevel (lblobs lbips)
  (let ((level (- *tax* 3))
        (redundantbips nil))
    (dolist (blob lblobs)
      (setf level (- level (- (numbclass (cadadr blob)) 3)))
      (setf redundantbips (ajoutebip redundantbips (extbips (cadadr blob)))))
    (- level (length lbips) (length redundantbips))))

(defun ajoutebip (lb1 lb2)
  (let ((res (copy-list lb1)))
    (dolist (b lb2 res)
      (unless (member b res :test 'bip=)
        (nconc-p res (list b))))))

(defun extbips (blob)
  (let ((res nil)
        (num (numbclass blob)))
    (dotimes (i num res)
      (unless (trivialbip i blob)
        (let ((vec (make-array *tax* :initial-element 0)))
          (dotimes (j *tax*)
            (when (= i (aref blob j))
              (setf (aref vec j) 1)))
          (nconc-p res (list (reorder vec))))))))

(defun trivialbip (i blob)
  (let ((res 0))
    (dotimes (j *tax* T)
      (when (= i (aref blob j))
        (setf res (+ 1 res)))
      (when (> res 1) (return nil)))))

(defun themax (otherblobs otherbips)
  (cond ((and (spec-< (caar otherblobs) *stop-addition*)
              (spec-< (caar otherbips) *stop-addition*))
         nil)
        ((null otherblobs) 'bip)
        ((null otherbips) 'blob)
        ((or (>= (caar otherblobs) *trsh-blobs*)
             (>= (caar otherblobs)(caar otherbips))) 
         'blob)
        (T 'bip)))

(defun spec-< (x1 x2)
(if x1 (< x1 x2) T))

(defun allcompat-blobbip (blob lbips)
  (dolist (bip lbips T)
    (unless (compatible-blobbip blob (cadr bip))
      (return nil))))

(defun compatible-blobbip (b1 b2)                 ; b1 is a blob, b2 is a bipartition ; dont check for redundancy
  (let* ((nc1 (numbclass b1))                     
         (tab (make-array (list (+ 1 nc1) 3) :initial-element 0))
         (ind nil))
    (dotimes (x *tax*)
      (when (= 0 (aref tab (aref b1 x)(aref b2 x)))
        (setf (aref tab (aref b1 x)(aref b2 x)) 1)
        (setf (aref tab (aref b1 x) 2)
          (+ 1 (aref tab (aref b1 x) 2)))
        (setf (aref tab nc1 (aref b2 x))
          (+ 1 (aref tab nc1 (aref b2 x))))))
    (and (or (= 1 (aref tab nc1 0))
             (= 1 (aref tab nc1 1)))
         (dotimes (i nc1 T)
           (cond ((= 1 (aref tab i 2)) nil)
                 ((null ind) (setf ind T))
                 (T (return nil)))))))

(defun clean-redundancy (lbips blob)
  (let ((res nil))
    (dolist (bip lbips res)
      (unless (bipinblob (cadr bip) blob)
        (nconc-p res (list bip))))))

(defun bipinblob (bip blob)
  (let ((lbips (extbips blob)))
    (dolist (b lbips nil)
      (when (bip= bip b)
        (return T)))))

(defun bip= (v1 v2)
  (let ((x (reorder v1))
        (y (reorder v2)))
    (dotimes (i *tax* T)
      (unless (= (aref x i)(aref y i))
        (return nil)))))

(defun best-cycle (blob)
  (let ((res nil)(ind 0))
    (dolist (b (cdr blob) res)
      (let ((x (cadr b))(val 0))
        (dolist (c (cdr blob))
          (let ((y (cadr c)))
            (when (cycle-equality x y)
              (setf val (+ val (car c))))))
        (when (> val ind)
          (setf ind val)
          (setf res (list val x)))))))

(defun 2-best-hybrid (vec allblobs)
(let ((lres nil)
      (vres (make-array 6))
      (numcl (numbclass vec)))
  (dotimes (i numcl)
    (let ((x (transfhybrid vec i))               
(count>4 0)
(count=4 0))
(dolist (blob allblobs)
  (dolist (c (cdr blob))
    (let ((y (cadr c)))
      (when (hybrid-equality x y)
        (if (> (numbclass y) 4)
            (setf count>4 (+ count>4 (car c)))
          (setf count=4 (+ count=4 (car c))))))))
(nconc-p lres (list (list i count>4 count=4)))))
(setf lres (sort lres 'comparehyb))
(setf lres (append (car lres)(cadr lres)))
(dotimes (i 6 vres)
  (setf (aref vres i)(nth i lres)))))

(defun comparehyb (l1 l2)
  (let ((c1> (cadr l1))
        (c1= (caddr l1))
        (c2> (cadr l2))
        (c2= (caddr l2)))
    (> (+ c1> C1=)(+ C2> C2=))))

(defun sizehyb (vec)
  (let ((res 0))
    (dotimes (i *tax* res)
      (when (= 0 (aref vec i))(setf res (+ res 1))))))


(defun transfhybrid (b k) 
  (let ((v (copy-v b)))
    (dotimes (i *tax* v)
      (cond ((= (aref v i) 0)(setf (aref v i) k))
       ((= (aref v i) k)(setf (aref v i) 0))
       (T T)))))  

(defun best-rooting (tabres)
  (let ((max nil)(ind4+ -1)(ind5+ -1)(indsize>4 nil))
    (dotimes (tax *tax* max)
      (let ((S4+ 0)(S5+ 0))
        (dotimes (i (array-dimension tabres 0))
          (if (> (numbclass (aref tabres i 0)) 4) (setf indsize>4 T)
            (setf indsize>4 nil))
          
          (cond ((and indsize>4
                      (= (aref tabres i 4)
                         (aref (aref tabres i 0) tax)))
                 (setf S4+ (+ S4+ (aref tabres i 8)(aref tabres i 9)))
                 (setf S5+ (+ S5+ (aref tabres i 8)(aref tabres i 9))))
           
           ((= (aref tabres i 4)
               (aref (aref tabres i 0) tax))
            (setf S4+ (+ S4+ (aref tabres i 8)(aref tabres i 9))))
           
           (indsize>4
            (setf S4+ (+ S4+ (aref tabres i 5)(aref tabres i 6)))
            (setf S5+ (+ S5+ (aref tabres i 5)(aref tabres i 6))))
           
           (T (setf S4+ (+ S4+ (aref tabres i 5)(aref tabres i 6))))))
        
        (when (or (> s4+ ind4+)(and (= s4+ ind4+)(> s5+ ind5+)))
          (setf max tax)
          (setf ind4+ s4+)
          (setf ind5+ s5+))))))

(defparameter *equality* nil)

(defun clean-tabres (root tabres allblobs)
  (let ((nblobs (array-dimension tabres 0)))
    (dotimes (i nblobs)
      (let ((vec (aref tabres i 0)))                                               ;select the root-compatible hybrid among 2
        (when (= (aref vec root)(aref tabres i 4))
          (print '********************************************************************)
          (setf *EQUALITY* T)
          (setf (aref tabres i 4)(aref tabres i 7))
          (setf (aref tabres i 5)(aref tabres i 8))
          (setf (aref tabres i 6)(aref tabres i 9)))
        (setf (aref tabres i 0) (re-root vec (aref tabres i 4)))                   ;re-root cycle with selected hybrid 
        (let ((suporigins (supp-separate-origins (aref tabres i 0) allblobs)))
          (setf (aref tabres i 7) (car suporigins))                                ;support of hybrid origins
          (setf (aref tabres i 8) (cadr suporigins))
          (setf (aref tabres i 9) (caddr suporigins)))))))

(defun supp-separate-origins (cycle allblobs)
  (let ((count-both 0)
        (count-left 0)
        (count-right 0)
        (numc (numbclass cycle)))
    (dolist (blob allblobs 
                  (list (max count-right count-left)
                        (min count-right count-left)
                        count-both))
      (dolist (c (cdr blob))
        (let ((y (cadr c)))
          (cond ((not (hybrid-equality cycle y)) T)
                ((or (origins-equality cycle y)
                     (origins-equality cycle (reverseord y numc)))
(setf count-both (+ count-both (car c)))
(setf count-left (+ count-left (car c)))
(setf count-right (+ count-right (car c))))
((or (left-origin-equality cycle y)
     (left-origin-equality cycle (reverseord y numc)))
 (setf count-left (+ count-left (car c))))
((or (right-origin-equality cycle y)
     (right-origin-equality cycle (reverseord y numc)))
 (setf count-right (+ count-right (car c))))))))))

(defun re-root (cycle hybrid)
  (let* ((l (length cycle))
         (num (numbclass cycle))
         (vec (make-array l)))
    (dotimes (i l vec)
      (setf (aref vec i)(mod (- (aref cycle i) hybrid) num))))) 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; EQUALITY OF CYCLES, HYBRIDS ETC ;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; EQUALITY OF CIRCULAR ORDERS ;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; SAME VECTOR ENCODING AS BLOBS ;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; ORDER = BLOCK INTEGER CODE ;;;;;;;;;;;;;;;;;;;; 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defun order-equality (v1 v2)
  (dotimes (x *tax* T)
    (unless (= (aref v1 x)(aref v2 x))
      (return nil))))

(defun cycle-equality (v1 v2)
  (let ((nc (numbclass v1))
        (ncb (numbclass v2)))
    (if (= nc ncb)
        (or (cycle=aux v1 v2 nc)
            (cycle=aux v1 (reverseord v2 nc) nc))
      nil)))

(defun cycle=aux (v1 v2 nc)
  (dotimes (i nc nil)
    (when (order-equality v1 (rotate v2 i nc))
      (return T))))

(defun rotate (v shift nc)
  (let ((res (make-array *tax*)))
    (dotimes (i *tax* res)
      (setf (aref res i)
        (mod (+ (aref v i) shift) nc)))))

(defun reverseord (v nc)
  (let ((res (make-array *tax*)))
    (dotimes (i *tax* res)
(setf (aref res i)
  (mod (- (aref v i)) nc)))))

(defun rooted-cyc= (lx ly)
  (or (order-equality (cadr lx)(cadr ly))
   (order-equality (cadr lx)
                   (reverseord (cadr ly) (numbclass (cadr ly))))))

(defun hybrid-equality (v1 v2) 
  (dotimes (i *tax* T)
    (when (or (and (= 0 (aref v1 i))
               (not (= 0 (aref v2 i))))
           (and (= 0 (aref v2 i))
                (not (= 0 (aref v1 i)))))
      (return nil))))

(defun origins-equality (V1 v2)
  (let ((n1 (- (numbclass v1) 1))
        (n2 (- (numbclass v2) 1)))
    (dotimes (i *tax* T)
      (when (or (and (= 1 (aref v1 i))
                     (not (= 1 (aref v2 i))))
                (and (= 1 (aref v2 i))
                     (not (= 1 (aref v1 i))))
                (and (= n1 (aref v1 i))
                     (not (= n2 (aref v2 i))))
                (and (= n2 (aref v2 i))
                     (not (= n1 (aref v1 i)))))
        (return nil)))))

(defun left-origin-equality (V1 v2)
  (dotimes (i *tax* T)
    (when (or (and (= 1 (aref v1 i))
                   (not (= 1 (aref v2 i))))
              (and (= 1 (aref v2 i))
                   (not (= 1 (aref v1 i)))))
      (return nil))))

(defun right-origin-equality (V1 v2)
  (let ((n1 (- (numbclass v1) 1))
        (n2 (- (numbclass v2) 1)))
    (dotimes (i *tax* T)
      (when (or (and (= n1 (aref v1 i))
                     (not (= n2 (aref v2 i))))
                (and (= n2 (aref v2 i))
                     (not (= n1 (aref v1 i)))))
        (return nil)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; DATA SET EASYNESS ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; SEARCH FOR ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; INPUT NETWORK WITH HIGHEST FREQUENCY ;;;;;;;;;;
;;;;;;;;;;;;;;;; WHEN THIS FREQUENCY IS LARGER THAN 50% ;;;;;;;;
;;;;;;;;;;;;;;;; THE DATASET IS VERY EASY ;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; EASINESS = FREQUENCY ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;; WE ALSO COMPUTE CONSENSUS NET FREQUENCY ;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(defun net_equality (i j)
  (and (lcycles-equality (nth i *all_cycles*)(nth j *all_cycles*))
       (lbips-equality (nth i *all_bipartitions*)(nth j *all_bipartitions*))))

(defun lbips-equality (lb1 lb2)
  (if (= (length lb1)(length lb2))
      (dolist (b1 lb1 T)
        (unless (member b1 lb2 :test 'bip=)
          (return nil)))
    nil))

(defun lcycles-equality (lc1 lc2)
  (if (= (length lc1)(length lc2))
      (dolist (c1 lc1 T)
        (unless (member c1 lc2 :test 'rcyc=)
          (return nil)))
    nil))

(defun rcyc= (x y)
  (or (order-equality x y)
      (order-equality x
                      (reverseord y (numbclass y)))))

(defparameter *freq* nil)

(defun net-frequency (ithnet)
  (let ((l (length *all_cycles*))
        (count 1))
    (setf (aref *freq* ithnet 0) ithnet)
    (dotimes (i l count)
      (when (and (> i ithnet)
                 (net_equality i ithnet))
(setf (aref *freq* i 0) ithnet) 
(setf count (+ 1 count))))
(dotimes (i l count)
  (when (and (aref *freq* i 0)(= (aref *freq* i 0) ithnet))
    (setf (aref *freq* i 1) count)))))

(defun easyness ()
  (let ((l (length *all_cycles*))
        (max 0))
    (setf *freq* (make-array (list l 2) :initial-element nil))
    (dotimes (i l max)
      (unless (aref *freq* i 0)
        (let ((f (net-frequency i)))
          (setf max (max max f)))))))  

(defun net-freq (lcycs lbips)
  (let ((l (length *all_cycles*)))
    (dotimes (i l 0)
      (when (and (lcycles-equality (nth i *all_cycles*) lcycs)
                 (lbips-equality (nth i *all_bipartitions*) lbips))
        (return (aref *freq* i 1))))))




       