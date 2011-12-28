;; large in-place transpose (with a scratch buffer for non-square
;; matrices)


(defconstant +transpose-base-size+ (ash 1 5))

(declaim (inline %transpose! %transpose-into))

(defun %transpose! (vec size start stride)
  (declare (type complex-sample-array vec)
           (type half-size size stride)
           (type index start)
           (optimize speed (safety 0)))
  (labels ((in-place (size start)
             (declare (type half-size size)
                      (type index start))
             (cond ((<= size +transpose-base-size+)
                    (loop
                      for i of-type half-index below size by 2
                      for start1 of-type index from start by (* 2 stride)
                      for start2 of-type index from start by 2
                      do (loop for start1 of-type index from start1
                               for start2 of-type index from start2 by stride
                               for j of-type half-index from i above 0
                               do (rotatef (ref vec start1)
                                           (ref vec start2))
                                  (rotatef (ref vec (+ start1 stride))
                                           (ref vec (+ start2 1)))
                               finally (rotatef (ref vec (+ start1 1))
                                                (ref vec (+ start2 stride)))))
                    vec)
                   (t
                    (let* ((size/2      (truncate size 2))
                           (long-stride (* size/2 stride)))
                      (in-place size/2
                                start)
                      (swap size/2
                            (+ start size/2)
                            (+ start long-stride))
                      (in-place size/2
                                (+ start long-stride size/2))))))
           (swap (size start1 start2)
             (declare (type half-size size)
                      (type index start1 start2))
             (cond ((<= size +transpose-base-size+)
                    (loop
                      for i from size above 0 by 4
                      for start1 of-type index from start1 by (* 4 stride)
                      for start2 of-type index from start2 by 4
                      do (loop for j from size above 0
                               for start1a of-type index from start1
                               for start1b of-type index from (+ start1 (* 2 stride))
                               for start2 of-type index from start2 by stride
                               do (rotatef (aref vec start1a)
                                           (aref vec start2))
                                  (rotatef (aref vec (+ start1a stride))
                                           (aref vec (+ start2 1)))
                               (rotatef (aref vec start1b)
                                        (aref vec (+ start2 2)))
                               (rotatef (aref vec (+ start1b stride))
                                        (aref vec (+ start2 3)))))
                    vec)
                   (t
                    (let* ((size/2 (truncate size 2))
                           (long-stride (* size/2 stride)))
                      (swap size/2
                            start1
                            start2)
                      (swap size/2
                            (+ start1 size/2)
                            (+ start2 long-stride))
                      (swap size/2
                            (+ start1 long-stride)
                            (+ start2 size/2))
                      (swap size/2
                            (+ start1 long-stride size/2)
                            (+ start1 long-stride size/2)))))))
    (in-place size start)))

(defun %transpose-into (dst src size startd starts strided strides)
  (declare (type complex-sample-array dst src)
           (type size size)
           (type index starts startd)
           (type half-size strided strides)
           (optimize speed (safety 0)))
  (labels ((rec (size startd starts)
             (declare (type half-size size)
                      (type index startd starts))
             (cond ((<= size +transpose-base-size+)
                    (loop
                      for i from size above 0 by 4
                      for startd of-type index from startd by (* 4 strided)
                      for starts of-type index from starts by 4
                      do (loop for j from size above 0
                               for startd1 of-type index from startd
                               for startd2 of-type index from (+ startd (* 2 strided))
                               for starts of-type index from starts by strides
                               do (setf (aref dst startd1) (aref src starts)
                                        (aref dst (+ startd1 strided)) (aref src (+ starts 1))
                                        (aref dst startd2)  (aref src (+ starts 2))
                                        (aref dst (+ startd2 strided)) (aref src (+ starts 3)))))
                    dst)
                   (t
                    (let* ((size/2 (truncate size 2))
                           (long-strided (* size/2 strided))
                           (long-strides (* size/2 strides)))
                      (rec size/2
                           startd
                           starts)
                      (rec size/2
                           (+ startd size/2)
                           (+ starts long-strides))
                      (rec size/2
                           (+ startd long-strided)
                           (+ starts size/2))
                      (rec size/2
                           (+ startd long-strided size/2)
                           (+ starts long-strides size/2)))))))
    (rec size startd starts)))

(defun transpose (vec tmp size1 size2 total vecs tmps)
  (declare (type complex-sample-array vec tmp)
           (type half-size size1 size2)
           (type size total)
           (type index vecs tmps))
  (flet ((%transpose-into (dst src size startd starts strided strides)
           (%transpose-into dst src size startd starts strided strides)))
    (cond ((= size1 size2)
           (return-from transpose (%transpose! vec size1 vecs size1)))
          ((< size1 size2)
           (let* ((size  size1)
                  (block (* size size)))
             (%transpose-into tmp vec size
                              0   0
                              size2 size1)
             (%transpose-into tmp vec size
                              size1 block
                              size2 size1)))
          (t
           (let* ((size  size2)
                  (block (* size size)))
             (%transpose-into tmp vec size
                              0   0
                              size2 size1)
             (%transpose-into tmp vec size
                              block size2
                              size2 size1)))))
  (loop for count from total above 0 by 4
        for vecs from vecs by 4
        for tmps from tmps by 4
        do (setf (ref vec (+ vecs 0)) (ref tmp (+ tmps 0))
                 (ref vec (+ vecs 1)) (ref tmp (+ tmps 1))
                 (ref vec (+ vecs 2)) (ref tmp (+ tmps 2))
                 (ref vec (+ vecs 3)) (ref tmp (+ tmps 3))))
  vec)

(declaim (notinline %transpose! %transpose-into))
