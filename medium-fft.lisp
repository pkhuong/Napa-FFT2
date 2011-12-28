(defconstant +default-blocking-factor+ 4)

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
                                (cooley-tukey 'cooley-tukey)
                                (blocking-factor +default-blocking-factor+))
  (assert (evenp (integer-length (1- size))))
  (assert (= 1 strided))
  (let ((half-size (ash 1 (truncate (integer-length (1- size))
                                    2))))
    (assert (<= blocking-factor half-size))
    `(flet ((sub-fft (dst src twiddle
                          startd starts)
              (declare (type complex-sample-array dst src twiddle)
                       (type index startd starts))
              ,(gen-fft/small half-size)))
       ;; copy columns from src to tmp
       (loop for i of-type index below ,half-size by ,blocking-factor
             for j of-type index from 0 by ,(* blocking-factor half-size)
             do (loop for count of-type index from ,half-size above 0
                      for j of-type index from ,startt
                      for k of-type index from (+ (* ,strides i) ,starts)
                        by (* ,strides ,half-size)
                      do (setf ,@(loop
                                   for block below blocking-factor
                                   append
                                   `((ref ,tmp (+ j ,(* block half-size)))
                                     (ref ,src (+ k (* ,block ,strides)))))))
             do (progn
                  ;; FFT columns in scratch space
                  ;; write result to dst, in rows
                  ,@(loop
                      for block below blocking-factor
                      collect
                      `(sub-fft ,dst ,tmp ,twiddle
                                (+ j ,startd ,(* block half-size))
                                (+ ,startt ,(* block half-size))))
                  (loop for i of-type index from (+ j ,startd)
                        for idx of-type index from j
                        for count below ,(* half-size blocking-factor)
                        do (setf (ref ,dst i)
                                 (* (ref ,dst i)
                                    (ref ,cooley-tukey idx))))))
       (loop for i of-type half-index below ,half-size by blocking-factor do
         (loop for count of-type half-index from ,half-size above 0
               for j of-type index from ,startt
               for k of-type index from (+ i ,startd) by ,half-size
               do (setf ,@(loop
                            for block below blocking-factor
                            append `((ref ,tmp (+ j ,(* block half-size)))
                                     (ref ,dst (+ k ,block))))))
         do
         (progn
           ,@(loop for block below blocking-factor
                   collect
                   `(sub-fft ,tmp ,tmp ,twiddle
                             (+ ,startt ,(* blocking-factor half-size)
                                        ,(* block half-size))
                             (+ ,startt ,(* block half-size))))
           (loop for count of-type half-index from ,half-size above 0
                 for j of-type index from (+ ,startt ,(* blocking-factor half-size))
                 for k of-type index from (+ i ,startd) by ,half-size
                 do (setf ,@(loop
                              for block below blocking-factor
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
                                 (cooley-tukey 'cooley-tukey)
                                 (blocking-factor +default-blocking-factor+))
  (let* ((size1 (ash 1 (truncate (integer-length (1- size))
                                 2)))
         (size2 (/ size size1)))
    (assert (integerp size2))
    (assert (<= blocking-factor size1))
    (assert (<= blocking-factor size2))
    `(progn
       (flet ((rec (dst src startd starts twiddle)
                ,(gen-fft/small size2
                                :dst 'dst :src 'src :twiddle 'twiddle
                                :startd 'startd
                                :starts 'starts
                                :strided strided)))
         ;; copy columns to scratch
         (loop for i of-type index below ,size1 by blocking-factor
               for j of-type index from 0 by ,(* blocking-factor size2)
               do (loop for count of-type index from ,size2 above 0
                        for j of-type index from ,startt
                        for k of-type index from (+ (* i ,strides) ,starts)
                          by (* ,size1 ,strides)
                        do (setf ,@(loop
                                     for block below blocking-factor
                                     append `((ref ,tmp (+ j ,(* block size2)))
                                              (ref ,src (+ k (* ,block ,strides)))))))
               do (progn
                    ,@(loop
                        for block below blocking-factor
                        collect `(rec ,dst ,tmp
                                      (+ (* ,strided (+ j ,(* block size2)))
                                         ,startd)
                                      (+ ,startt ,(* block size2))
                                      ,twiddle))
                    (loop for i of-type index from (+ (* j ,strided) ,startd) by ,strided
                          for idx of-type index from j
                          for count of-type index below ,(* size2 blocking-factor)
                          do (setf (ref ,dst i)
                                   (* (ref ,dst i)
                                      (ref ,cooley-tukey idx)))))))
       (flet ((rec (vec startd starts twiddle)
                ,(gen-fft/small size1
                                :dst 'vec :src 'vec :twiddle 'twiddle
                                :startd 'startd
                                :starts 'starts)))
         (loop for i of-type index below ,size2 by blocking-factor do
           (loop for count of-type half-index from ,size1 above 0
                 for j of-type index from ,startt
                 for k of-type index from (+ (* i ,strided) ,startd)
                   by (* ,size2 ,strided)
                 do (setf ,@(loop
                              for block below blocking-factor
                              append
                              `((ref ,tmp (+ j ,(* block size1)))
                                (ref ,dst (+ k (* ,block ,strided)))))))
           do (progn
                ,@(loop
                    for block below blocking-factor
                    collect `(rec ,tmp
                                  (+ ,startt ,(* blocking-factor size1)
                                     ,(* block size1))
                                  (+ ,startt ,(* block size1))))
                (loop for count of-type half-index from ,size1 above 0
                      for j of-type index
                        from (+ ,startt ,(* blocking-factor size1))
                      for k of-type index from (+ (* i ,strided) ,startd)
                        by (* ,size2 ,strided)
                      do (setf
                          ,@(loop
                              for block below blocking-factor
                              append
                              `((ref ,dst (+ k ,block))
                                (ref ,tmp (+ j (* ,strided
                                                  ,(* block size1)))))))))))
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
                         (cooley-tukey 'cooley-tukey)
                         (blocking-factor +default-blocking-factor+))
  (declare (ignore dst src tmp
                   startd starts startt
                   strides
                   twiddle cooley-tukey
                   blocking-factor))
  (if (and (evenp (integer-length (1- size)))
           (eql strided 1))
      (apply 'gen-simple-fft/medium size args)
      (apply 'gen-generic-fft/medium size args)))

