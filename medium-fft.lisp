(defconstant +blocking-factor+ 4)

(defun gen-simple-fft/medium (size
                              &key (dst 'dst)
                                (src 'src)
                                (tmp 'tmp)
                                (startd 'startd)
                                (starts 'starts)
                                (startt 'startt)
                                (strides 1)
                                (strided 1)
                                (twiddle 'twiddle)
                                (cooley-tukey 'cooley-tukey))
  (assert (evenp (integer-length (1- size))))
  (assert (= 1 strided))
  (let ((half-size (ash 1 (truncate (integer-length (1- size))
                                    2))))
    `(flet ((sub-fft (dst src twiddle
                          startd starts)
              (declare (type complex-sample-array dst src twiddle)
                       (type index startd starts))
              ,(gen-fft/small half-size)))
       (loop for i of-type index below ,half-size by ,+blocking-factor+
             for j of-type index from 0 by ,(* +blocking-factor+ half-size)
             do (loop for count of-type index from ,half-size above 0
                      for j from ,startt
                      for k of-type index from (+ i ,starts) by ,half-size
                      do (setf ,@(loop
                                   for block below +blocking-factor+
                                   append
                                   `((ref ,tmp (+ j ,(* block half-size)))
                                     (ref ,src (* (+ k ,block) ,strides))))))
             do (progn
                  ,@(loop
                      for block below +blocking-factor+
                      collect
                      `(sub-fft ,dst ,tmp ,twiddle
                                (+ j ,startd ,(* block half-size))
                                (+ ,startt ,(* block half-size))))
                  (loop for i of-type index from (+ j ,startd)
                        for idx of-type index from j
                        for count below ,(* half-size +blocking-factor+)
                        do (setf (ref ,dst i)
                                 (* (ref ,dst i)
                                    (ref ,cooley-tukey idx))))))
       (loop for i of-type index below ,half-size by +blocking-factor+ do
         (loop for count from ,half-size above 0
               for j of-type index from ,startt
               for k of-type index from (+ i ,startd) by ,half-size
               do (setf ,@(loop
                            for block below +blocking-factor+
                            append `((ref ,tmp (+ j ,(* block half-size)))
                                     (ref ,dst (+ k ,block))))))
         do
         (progn
           ,@(loop for block below +blocking-factor+
                   collect
                   `(sub-fft ,tmp ,tmp ,twiddle
                             (+ ,startt ,(* +blocking-factor+ half-size)
                                        ,(* block half-size))
                             (+ ,startt ,(* block half-size))))
           (loop for count from ,half-size above 0
                 for j of-type index from (+ ,startt ,(* +blocking-factor+ half-size))
                 for k of-type index from (+ i ,startd) by ,half-size
                 do (setf ,@(loop
                              for block below +blocking-factor+
                              append `((ref ,dst (+ k ,block))
                                       (ref ,tmp
                                            (+ j ,(* block half-size)))))))))
       dst)))

(defun gen-generic-fft/medium (size
                               &key (dst 'dst)
                                 (src 'src)
                                 (tmp 'tmp)
                                 (startd 'startd)
                                 (starts 'starts)
                                 (startt 'startt)
                                 (strides 1)
                                 (strided 1)
                                 (twiddle 'twiddle)
                                 (cooley-tukey 'cooley-tukey))
  (let* ((size1 (ash 1 (truncate (integer-length (1- size))
                                 2)))
         (size2 (/ size size1)))
    (assert (integerp size2))
    `(progn
       (flet ((rec (startd starts)
                ,(gen-fft/small size2
                                :dst dst :src tmp :twiddle twiddle
                                :startd 'startd
                                :starts 'starts
                                :strided strided)))
         (loop for i of-type index below ,size1 by +blocking-factor+
               for j of-type index from 0 by ,(* +blocking-factor+ size2)
               do (loop for count of-type index from ,size2 above 0
                        for j from ,startt
                        for k of-type index from (+ i ,starts) by ,size1
                        do (setf ,@(loop
                                     for block below +blocking-factor+
                                     append `((ref ,tmp (+ j ,(* block size2)))
                                              (ref ,src (* (+ k ,block) ,strides))))))
               do (progn
                    ,@(loop
                        for block below +blocking-factor+
                        collect `(rec (+ j ,(* block size2) ,startd)
                                      (+ ,startt ,(* block size2))))
                    (loop for i of-type index from (+ j ,startd) by ,strided
                          for idx of-type index from j
                          for count below ,(* size2 +blocking-factor+)
                          do (setf (ref ,dst i)
                                   (* (ref ,dst i)
                                      (ref ,cooley-tukey idx)))))))
       (flet ((rec (startd starts)
                ,(gen-fft/small size1
                                :dst tmp :src tmp :twiddle twiddle
                                :startd 'startd
                                :starts 'starts)))
         (loop for i of-type index below ,size2 by +blocking-factor+ do
           (loop for count from ,size1 above 0
                 for j of-type index from ,startt
                 for k of-type index from (+ i ,startd) by (* ,size2 ,strided)
                 do (setf ,@(loop
                              for block below +blocking-factor+
                              append `((ref ,tmp (+ j ,(* block size1)))
                                       (ref ,dst (+ k ,block))))))
           do (progn
                ,@(loop
                    for block below +blocking-factor+
                    collect `(rec (+ ,startt ,(* +blocking-factor+ size1)
                                     ,(* block size1))
                                  (+ ,startt ,(* block size1))))
                (loop for count from ,size1 above 0
                      for j of-type index
                        from (+ ,startt ,(* +blocking-factor+ size1))
                      for k of-type index from (+ i ,startd) by (* ,size2 ,strided)
                      do (setf
                          ,@(loop for block below +blocking-factor+
                                  append `((ref ,dst (+ k ,block))
                                           (ref ,tmp (+ j ,(* block size1))))))))))
       dst)))

(defun gen-fft/medium (size &rest args
                       &key (dst 'dst)
                         (src 'src)
                         (tmp 'tmp)
                         (startd 'startd)
                         (starts 'starts)
                         (startt 'startt)
                         (strides 1)
                         (strided 1)
                         (twiddle 'twiddle)
                         (cooley-tukey 'cooley-tukey))
  (declare (ignore dst src tmp
                   startd starts startt
                   strides
                   twiddle cooley-tukey))
  (if (and (evenp (integer-length (1- size)))
           (eql strided 1))
      (apply 'gen-square-fft/medium size args)
      (apply 'gen-generic-fft/medium size args)))

#+nil
(let ((fun (compile nil `(lambda (dst src twiddle startd starts)
                           ,(gen-fft/small 16))))
      (twiddle (bordeaux-fft::make-twiddle-factors 16 1))
      (ck-factors (bordeaux-fft::make-cooley-tuckey-factors 16 16 1))
      (src *vec*)
      (dst (make-array 256 :element-type 'complex-sample))
      (tmp (make-array 256 :element-type 'complex-sample)))
  (medium-fft fun twiddle 256
              src dst tmp ck-factors))

#+nil
(defun medium-fft (small-fft twiddle size
                   vec dst tmp ck-factors)
  (declare (type function small-fft)
           (type size size)
           (type complex-sample-array vec dst tmp
                 twiddle ck-factors)
           (optimize speed (safety 0)))
  (let ((half-size (ash 1 (truncate (integer-length (1- size))
                                    2))))
    (loop for i of-type index below half-size
          for dsts of-type index by half-size
          do (loop for j of-type index below half-size
                   for k of-type index from i by half-size
                   do (setf (ref tmp j) (ref vec k)))
             (funcall small-fft
                      dst  tmp twiddle
                      dsts 0)
             (loop for i from dsts
                   for count from half-size above 0
                   do (setf (ref dst i)
                            (* (ref dst i)
                               (ref ck-factors i)))))
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

