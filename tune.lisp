(defun run-all-instances (function block-size instances)
  (declare (type simple-vector instances)
           (type function function)
           (type index block-size))
  (let* ((len   (length instances))
         (times (make-array (ceiling len block-size)))
         last-value)
    (loop for i below len by block-size
          for j upfrom 0
          for last of-type index = (min len (+ i block-size))
          do (multiple-value-bind (value time)
                 (sb-vm::with-cycle-counter
                   (let (dst)
                     (loop for i from i below last
                           do (setf dst (funcall function (aref instances j)))
                           finally (return dst))))
               (setf (aref times j) (/ (float time 1d0)
                                       (- last i))
                     last-value     value)))
    (sort times #'<)
    (values (aref times (truncate (length times) 20))
            last-value)))

(defun make-inputs (size)
  (let ((count (max 2 (truncate (ash 32 20) ;; at least 32 MB
                                (* size 16)))))
    (map-into (make-array count)
              (let ((count 0))
                (lambda ()
                  (fill (make-array size :element-type 'complex-sample)
                        (complex (incf count) 0d0)))))))

(defun time-inputs (makers n &optional (block-size 16))
  (when (atom makers)
    (setf makers (list makers)))
  (let ((funs (mapcar (lambda (maker)
                        (if (atom maker)
                            (funcall maker n)
                            (apply (first maker) n (rest maker))))
                      makers)))
    (flet ((run-it ()
             (let* ((outputs '())
                    (cycles
                      (loop with inputs = (make-inputs n)
                            for fun in funs
                            collect (multiple-value-bind (cycles dst)
                                        (run-all-instances fun
                                                           block-size
                                                           inputs)
                                      (push dst outputs)
                                      cycles))))
               (loop for maker in makers
                     for output in (rest outputs)
                     unless (every (lambda (x y)
                                     (< (/ (abs (- x y))
                                           (max x y 1d0))
                                        1d-5))
                                   (first outputs) output)
                       do (format t "Mismatch for ~A~%" maker))
               cycles)))
      (mapcar #'min (run-it) (run-it)))))

(defun make-bordeaux-fft (n)
  (let ((dst (make-array n :element-type 'complex-sample))
        (instance (bordeaux-fft::make-fft-instance n)))
    (lambda (src)
      (bordeaux-fft::fft-common instance src dst))))

(defun make-small-fft (n)
  (let* ((fft (compile nil `(lambda (dst src twiddle)
                              (declare (type complex-sample-array dst src
                                             twiddle)
                                       (ignorable twiddle)
                                       (optimize speed (safety 0)))
                              ,(gen-fft/small n :startd 0 :starts 0 :twiddle 'twiddle))))
         (twiddle (make-twiddle-factors n 1))
         (dst     (make-array n :element-type 'complex-sample)))
    (declare (type function fft))
    (lambda (src)
      (funcall fft dst src twiddle))))

(defun make-medium-fft (n)
  (let* ((fft (compile nil `(lambda (dst src tmp twiddle ck)
                              (declare (type complex-sample-array dst src
                                             tmp twiddle ck)
                                       (ignorable twiddle)
                                       (optimize speed (safety 0)))
                              ,(gen-fft/medium n :startd 0 :starts 0
                                                 :startt 0
                                                 :twiddle 'twiddle
                                                 :cooley-tukey 'ck))))
         (twiddle (make-twiddle-factors n 1)) ;; only need sqrt
         (ck      (make-all-factors (integer-length (1- n)) 1))
         (dst     (make-array n :element-type 'complex-sample))
         (tmp     (make-array n :element-type 'complex-sample)))
    (declare (type function fft))
    (lambda (src)
      (funcall fft dst src tmp twiddle ck))))

(defun make-large-fft (n &optional (lower 'gen-fft/medium))
  (let* ((fft (compile nil `(lambda (dst src tmp twiddle ck)
                              (declare (type complex-sample-array dst src
                                             tmp twiddle ck)
                                       (ignorable twiddle tmp)
                                       (optimize speed (safety 0)))
                              ,(gen-fft/large n
                                              :startd 0
                                              :starts 0
                                              :startt 0
                                              :twiddle 'twiddle
                                              :cooley-tukey 'ck
                                              :lower lower))))
         (twiddle (make-twiddle-factors n 1)) ;; only need sqrt
         (ck      (make-all-factors (integer-length (1- n)) 1))
         (dst     (make-array n :element-type 'complex-sample))
         (tmp     (make-array n :element-type 'complex-sample)))
    (declare (type function fft))
    (lambda (src)
      (funcall fft dst src tmp twiddle ck))))
