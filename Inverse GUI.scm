(define (make-new-rpvar name default type)
	(if (not (rp-var-object name))
		(rp-var-define name default type #f)
	)
)

; 定义2个用于交互的变量
(make-new-rpvar 'myudf/sigma 0.01 'real)
(make-new-rpvar 'myudf/co 100 'real)

; 定义一个函数，用于处理GUI界面
(define gui-dialog-box
	(let(
			(dialog-box #f)
			(table)
			(sigma)
			(co)
			(btn1)
			(btn2)
		)

		; 定义回调函数，加载后自动调用
		(define (update-cb . args)
			(cx-set-real-entry sigma 0.01)
			(cx-set-real-entry co 100.0)
		)

		; 定义回调函数，点击按钮OK后调用
		(define (apply-cb . args)
			(display "The dialog closed!")
		)  
		
		; 定义回调函数，点击按钮Click后调用
		(define (btn-cb1 . args)
			(rpsetvar 'myudf/sigma (cx-show-real-entry sigma)) ; 将GUI中的控件值赋值给变量
			(rpsetvar 'myudf/co (cx-show-real-entry co))  ; 将GUI中的控件值赋值给变量
			(%run-udf-apply 7) ; 调用UDF宏，这里的7为UDF宏中的mode参数
		)
		
		(define (btn-cb2 . args)
			(display "\n")
			(display "******************\n")
			(display "*Have a good day!*\n")
			(display "******************\n")
		)
		
       ; 创建GUI元素
		(lambda args
			(if (not dialog-box)
				(let ()
					(set! dialog-box (cx-create-panel "Source Inverse" apply-cb update-cb))
					(set! table (cx-create-table dialog-box "Some useful parms" 'border #t 'below 0 'right-of 0))
			   
					(set! sigma (cx-create-real-entry table "For sigma:" 'row 1 'col 1))
					(set! co (cx-create-real-entry table "For Co:" 'row 1 'col 2)) 
					(set! btn1 (cx-create-button table "Apply" 'activate-callback btn-cb1 'row 2 'col 1))
					(set! btn2 (cx-create-button table "ZHHL"'activate-callback btn-cb2 'row 2 'col 2))					
				) 
			)         
			(cx-show-panel dialog-box)
		) 
    ) 
) 
 
(gui-dialog-box)
