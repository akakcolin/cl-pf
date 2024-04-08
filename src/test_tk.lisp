(require 'asdf)
(require 'ltk)                                              ; 1
(use-package :ltk)                                          ; 2



(defun hello-1()
  (with-ltk ()
   (let ((b (make-instance 'button
                           :master nil
                           :text "Press Me"
                           :command (lambda ()
                                      (format t "Hello World!~&")))))
     (pack b))))

(hello-1)

(defun hello-2()
  (with-ltk ()
   (let* ((f (make-instance 'frame))
          (b1 (make-instance 'button
                             :master f
                             :text "Button 1"
                             :command (lambda () (format t "Button1~&"))))
          (b2 (make-instance 'button
                             :master f
                             :text "Button 2"
                             :command (lambda () (format t "Button2~&")))))
     (pack f)
     (pack b1 :side :left)
     (pack b2 :side :left)
     (configure f :borderwidth 3)
     (configure f :relief :sunken)
     )))

(hello-2)

(with-ltk ()
  (let ((button (make-instance 'button
                               :text "hello"
                               :command (lambda ()
                                          (format t "clicked")))))
    (grid button 0 0)))
