(defun gen-square-fft/large (half-size
                             &key (dst 'dst)
                               (src 'src)
                               (tmp 'tmp)
                               (startd 'startd)
                               (starts 'starts)
                               (startt 'startt)
                               (scale  1d0)
                               (twiddle 'twiddle)
                               (cooley-tukey-large 'cooley-tukey-large)
                               (cooley-tukey-size1 'cooley-tukey-size1)
                               (cooley-tukey-size2 'cooley-tukey-size2)
                               &aux (scalep (not (onep scale))))
  (declare (ignore cooley-tukey-size2))
  `(flet ((rec (dst src tmp startd starts startt
                twiddle ck)
            (declare (type complex-sample-array dst src tmp)
                     (type index startd starts startt)
                     (type complex-sample-array twiddle ck))
            ,(gen-fft/medium half-size :dst 'tmp
                                       :src 'src
                                       :tmp 'tmp
                                       :startd 'startd
                                       :starts 'starts
                                       :startt 'startt
                                       :twiddle 'twiddle
                                       :cooley-tukey 'ck)))
     (loop for i of-type index below size by half-size
           do (rec ,dst ,src ,tmp
                   (+ ,startd i)
                   (+ ,starts i)
                   ,startt
                   ,twiddle
                   ,cooley-tukey-size1))
     ,(generate-transpose half-size half-size nil
                          :vec dst :tmp tmp
                          :vecs startd :tmps startt
                          :twiddle cooley-tukey-large
                          :twiddle-start 0)
     (loop for i of-type index below size by half-size
           do (let ((start-dst (+ ,startd i))
                    (start-tmp (+ ,startt ,half-size))
                    ,@(and scalep '((.scale2.    0d0))))
                (declare (type index start-dst start-tmp)
                         ,@(and scalep '((type double-float scale2))))
                ,(and scalep `(setf .scale2. ,scale))
                (rec ,tmp ,dst ,tmp
                     ,startt
                     start-dst
                     start-tmp
                     ,twiddle
                     ,cooley-tukey-size1)
                ,(generate-blit half-size
                                :dst dst :src tmp
                                :startd 'start-dst :starts 'start-tmp
                                :scale (if scalep '.scale2. 1d0))))))

(defun gen-rect-fft/large (size1 size2
                             &key (dst 'dst)
                               (src 'src)
                               (tmp 'tmp)
                               (startd 'startd)
                               (starts 'starts)
                               (startt 'startt)
                               (scale  1d0)
                               (twiddle 'twiddle)
                               (cooley-tukey-large 'cooley-tukey-large)
                               (cooley-tukey-size1 'cooley-tukey-size1)
                               (cooley-tukey-size2 'cooley-tukey-size2)
                               &aux (scalep (not (onep scale))))
  `(progn
     (flet ((rec (dst src tmp startd starts startt
                      twiddle ck)
              (declare (type complex-sample-array dst src tmp)
                       (type index startd starts startt)
                       (type complex-sample-array twiddle ck))
              ,(gen-fft/medium size2 :dst 'tmp
                                     :src 'src
                                     :tmp 'tmp
                                     :startd 'startd
                                     :starts 'starts
                                     :startt 'startt
                                     :twiddle 'twiddle
                                     :cooley-tukey 'ck)))
       (loop for i of-type index below size by size2
             do (rec ,dst ,src ,tmp
                     (+ ,startd i)
                     (+ ,starts i)
                     ,startt
                     ,twiddle
                     ,cooley-tukey-size2)))
     ,(generate-transpose size2 size1 nil
                          :vec dst :tmp tmp
                          :vecs startd :tmps startt
                          :twiddle cooley-tukey-large
                          :twiddle-start 0)
     (flet ((rec (dst src tmp startd starts startt
                      twiddle ck)
              (declare (type complex-sample-array dst src tmp)
                       (type index startd starts startt)
                       (type complex-sample-array twiddle ck))
              ,(gen-fft/medium size1 :dst 'tmp
                                     :src 'src
                                     :tmp 'tmp
                                     :startd 'startd
                                     :starts 'starts
                                     :startt 'startt
                                     :twiddle 'twiddle
                                     :cooley-tukey 'ck)))
       (loop for i of-type index below size by size1
             do (let ((start-dst (+ ,startd i))
                      (start-tmp (+ ,startt ,size1))
                      ,@(and scalep '((.scale2.    0d0))))
                  (declare (type index start-dst start-tmp)
                           ,@(and scalep '((type double-float scale2))))
                  ,(and scalep `(setf .scale2. ,scale))
                  (rec ,tmp ,dst ,tmp
                       ,startt
                       start-dst
                       start-tmp
                       ,twiddle
                       ,cooley-tukey-size1)
                  ,(generate-blit size1
                                  :dst dst :src tmp
                                  :startd 'start-dst :starts 'start-tmp
                                  :scale (if scalep '.scale2. 1d0)))))))

(defun gen-fft/large (size &rest args
                      &key (dst 'dst)
                        (src 'src)
                        (tmp 'tmp)
                        (startd 'startd)
                        (starts 'starts)
                        (startt 'startt)
                        (scale  1d0)
                        (twiddle 'twiddle)
                        (cooley-tukey-large 'cooley-tukey-large)
                        (cooley-tukey-size1 'cooley-tukey-size1)
                        (cooley-tukey-size2 'cooley-tukey-size2))
  (declare (ignore dst src tmp startd starts startt
                   scale twiddle
                   cooley-tukey-large cooley-tukey-size1
                   cooley-tukey-size2))
  (let* ((size1 (ash 1 (truncate (integer-length (1- size))
                                 2)))
         (size2 (/ size size1)))
    (assert (integerp size2))
    (if (= size1 size2)
        (apply 'gen-square-fft/large size1 args)
        (apply 'gen-rect-fft/large size1 size2 args))))
