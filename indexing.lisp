;;; Aggressive constant-folding on vector indexing operations
;;;
;;; Only handles simple expression involving addition and multiplication
;;; by constants.

(declaim (inline ref (setf ref)))
(defun ref (array index)
  (aref array index))

(defun (setf ref) (value array index)
  (setf (aref array index) value))

(defvar *current-multiplier*)

(defun simplify-multiplied-expressions (term-table)
  (let ((values '()))
    (maphash (lambda (k v)
               (push `(* ,(if (rest v)
                              `(+ ,@v)
                              (first v))
                         ,k)
                     values))
             term-table)
    (let ((sorted (sort values (lambda (x y)
                                 (if (and (numberp x)
                                          (numberp y))
                                     (< x y)
                                     (numberp x)))
                        :key #'third)))
      (assert (not (eql 1 (third (first sorted)))))
      sorted)))

(defun split-noises (noises env)
  (let ((constants '())
        (variants  '()))
    (dolist (noise noises (values (nreverse constants)
                                  (nreverse variants)))
      (if (constantp noise env)
          (push noise constants)
          (push noise variants)))))

(defun build-index-expression (offset constant rest expressions)
  (let ((addends '()))
    (cond ((and (/= offset 0) constant)
           (push `(+ ,offset ,@constant) addends))
          ((/= offset 0)
           (push offset addends))
          ((rest constant)
           (push `(+ ,@constant) addends))
          (constant
           (push constant addends)))
    (let ((variant (append expressions rest)))
      (cond ((rest variant)
             (push `(+ ,@variant) addends))
            (variant
             (push (first variant) addends))))
    (cond ((rest addends)
           `(+ ,@addends))
          (addends
           (first addends))
          (t 0))))

(defun simplify-indexing-expression (index env)
  (let ((terms  (make-hash-table :test #'equal)) ; multiplier -> operands
        (noise '())
        (offset 0)
        (*current-multiplier* 1))
    (labels ((simplify-atom (term)
               (cond #+sbcl
                     ((constantp term env)
                      (incf offset (* *current-multiplier*
                                      (sb-int:constant-form-value term env))))
                     #-sbcl
                     ((integerp term)
                      (incf offset (* *current-multiplier* term)))
                     #-sbcl
                     ((constantp term env)
                      (push `(* ,term ,*current-multiplier*) noise))
                     ((= 1 *current-multiplier*)
                      (push term noise))
                     (t
                      (push term (gethash *current-multiplier* terms)))))
             (simplify+ (terms)
               (map nil #'simplify terms))
             (simplify* (terms)
               (let ((literal    *current-multiplier*)
                     (constants '())
                     (others    '()))
                 (dolist (term terms)
                   (cond #+sbcl
                         ((constantp term env)
                          (setf literal (* literal
                                           (sb-int:constant-form-value term env))))
                         #-sbcl
                         ((integerp term)
                          (setf literal (* literal term)))
                         #-sbcl
                         ((constantp term env)
                          (push term constants))
                         (t
                          (push term others))))
                 (when (zerop literal)
                   (return-from simplify*))
                 (setf constants (nreverse constants)
                       others    (nreverse others))
                 (cond ((rest others)
                        (push `(* ,@others (* ,literal ,@constants)) noise))
                       ((and constants others)
                        (push (first others)
                              (gethash `(* ,literal ,@constants) terms)))
                       (constants
                        (push `(* ,literal ,@constants) noise))
                       (t
                        (let ((*current-multiplier* literal))
                          (simplify (or (first others) 1)))))))
             (simplify  (term)
               (typecase term
                 ((cons (eql *))
                  (simplify* (rest term)))
                 ((cons (eql +))
                  (simplify+ (rest term)))
                 (t
                  (simplify-atom term)))))
      (simplify index))
    (multiple-value-call #'build-index-expression
      offset (split-noises noise env)
      (simplify-multiplied-expressions terms))))

(define-compiler-macro ref (array index &environment env)
  `(aref ,array ,(simplify-indexing-expression index env)))

(define-compiler-macro (setf ref) (value array index &environment env)
  (let ((_value (gensym "VALUE")))
    `(let ((,_value ,value))
       (setf (aref ,array ,(simplify-indexing-expression index env))
             ,_value))))
