(defconstant +blocking-factor+ 1)

(defun medium-fft (small-fft twiddle size
                   vec dst tmp ck-factors)
  (declare (type function small-fft)
           (type size size)
           (type complex-sample-array vec dst tmp
                 twiddle ck-factors))
  (let ((half-size (ash 1 (truncate (integer-length (1- size))
                                    2))))
    (loop for i of-type index below half-size do
      (loop for j of-type index below half-size
            for k of-type index from i by half-size
            do (setf (ref tmp j) (ref vec k)))
      (funcall small-fft
               dst tmp twiddle
               (* i half-size)   0))
    (map-into dst #'* dst ck-factors)
    (loop for i of-type index below half-size do
      (loop for j of-type index below half-size
            for k of-type index from i by half-size
            do (setf (ref tmp j) (ref dst k)))
      (funcall small-fft
               tmp tmp twiddle
               half-size 0)
      (loop for j of-type index from half-size
            for k of-type index from i by half-size
            below size
            do (setf (ref dst k) (ref tmp j))))
    dst))

(let ((fun (compile nil `(lambda (dst src twiddle startd starts)
                           ,(gen-fft/n 16))))
      (twiddle (bordeaux-fft::make-twiddle-factors 16 1))
      (ck-factors (bordeaux-fft::make-cooley-tuckey-factors 16 16 1))
      (src *vec*)
      (dst (make-array 256 :element-type 'complex-sample))
      (tmp (make-array 256 :element-type 'complex-sample)))
  (medium-fft fun twiddle 256
              src dst tmp ck-factors))
