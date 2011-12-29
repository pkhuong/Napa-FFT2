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
                               &aux (scalep (not (onep scale)))
                                 (size (* half-size half-size)))
  (declare (ignore cooley-tukey-size2))
  `(flet ((rec (dst src tmp startd starts startt
                twiddle ck)
            (declare (type complex-sample-array dst src tmp)
                     (type index startd starts startt)
                     (type complex-sample-array twiddle ck))
            ,(gen-fft/medium half-size :dst 'dst
                                       :src 'src
                                       :tmp 'tmp
                                       :startd 'startd
                                       :starts 'starts
                                       :startt 'startt
                                       :twiddle 'twiddle
                                       :cooley-tukey 'ck)))
     (loop for i of-type index below ,size by ,half-size
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
     (loop for i of-type index below ,size by ,half-size
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
                                :startd 'start-dst :starts startt
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
                               &aux (scalep (not (constantp scale)))
                                 (size (* size1 size2)))
  `(progn
     (flet ((rec (dst src tmp startd starts startt
                      twiddle ck)
              (declare (type complex-sample-array dst src tmp)
                       (type index startd starts startt)
                       (type complex-sample-array twiddle ck))
              ,(gen-fft/medium size1 :dst 'dst
                                     :src 'src
                                     :tmp 'tmp
                                     :startd 'startd
                                     :starts 'starts
                                     :startt 'startt
                                     :twiddle 'twiddle
                                     :cooley-tukey 'ck)))
       (loop for i of-type index below ,size by ,size1
             do (rec ,dst ,src ,tmp
                     (+ ,startd i)
                     (+ ,starts i)
                     ,startt
                     ,twiddle
                     ,cooley-tukey-size1)))
     ,(generate-transpose size1 size2 nil
                          :vec dst :tmp tmp
                          :vecs startd :tmps startt
                          :twiddle cooley-tukey-large
                          :twiddle-start 0)
     (flet ((rec (dst src tmp startd starts startt
                      twiddle ck)
              (declare (type complex-sample-array dst src tmp)
                       (type index startd starts startt)
                       (type complex-sample-array twiddle ck))
              ,(gen-fft/medium size2 :dst 'dst
                                     :src 'src
                                     :tmp 'tmp
                                     :startd 'startd
                                     :starts 'starts
                                     :startt 'startt
                                     :twiddle 'twiddle
                                     :cooley-tukey 'ck)))
       (loop for i of-type index below ,size by ,size2
             do (let ((start-dst (+ ,startd i))
                      (start-tmp (+ ,startt ,size2)))
                  (declare (type index start-dst start-tmp))
                  (rec ,tmp ,dst ,tmp
                       ,startt
                       start-dst
                       start-tmp
                       ,twiddle
                       ,cooley-tukey-size2)
                  (let ,(and scalep '((scale2 0d0)))
                    ,@(and scalep `((declare (type double-float scale2))
                                    (setf scale2 ,scale)))
                    ,(generate-blit size2
                                    :dst dst :src tmp
                                    :startd 'start-dst :starts startt
                                    :scale (if scalep 'scale2 scale))))))))

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

(defun split-int (x)
  (let ((len (integer-length (1- x))))
    (values (ash 1 (truncate len 2))
            (ash 1 (truncate (1+ len) 2)))))

(defun test-fft (size input)
  (let* ((size1 (ash 1 (truncate (integer-length (1- size))
                                 2)))
         (size2 (/ size size1))
         (twiddle (bordeaux-fft::make-twiddle-factors size2 1))
         (ck      (make-cooley-tukey-factors size1 size2 1))
         (ck1     (multiple-value-call #'make-cooley-tukey-factors
                    (split-int size1) 1))
         (ck2     (multiple-value-call #'make-cooley-tukey-factors
                    (split-int size2) 1))
         (tmp     (make-array size :element-type 'complex-sample))
         (dst     (make-array size :element-type 'complex-sample)))
    (let ((transpose (compile nil `(lambda (vec tmp)
                                     (declare (type complex-sample-array vec tmp)
                                              (optimize speed (safety 0)))
                                     ,(generate-transpose size2 size1 nil
                                                          :vecs 0 :tmps 0))))
          (fft       (compile nil `(lambda (src dst tmp
                                            twiddle ck ck1 ck2)
                                     (declare (type complex-sample-array src tmp dst
                                                    twiddle ck ck1 ck2)
                                              (optimize speed (safety 0))
                                              (ignorable twiddle ck1 ck2))
                                     ,(gen-fft/large size :startd 0 :starts 0 :startt 0
                                                          :twiddle 'twiddle
                                                          :cooley-tukey-large 'ck
                                                          :cooley-tukey-size1 'ck1
                                                          :cooley-tukey-size2 'ck2)))))
      (time
       (progn
         (funcall transpose input tmp)
         (funcall fft input dst tmp
                  twiddle ck ck1 ck2 .5d0)
         (funcall transpose dst tmp)))
      dst)))


#+nil
(reduce #'max (map '(simple-array double-float 1)
                   (lambda (x y)
                     (abs (- x y)))
                   *x* *y*))