; Author: WHU ZHHL
; email: zhhl_email@qq.com
; data: 2022-04-01
; 用于辅助伴随概率法寻源得计算

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

		; 定义回调函数，加载后自动调用，用于初始化
		(define (update-cb . args)
			(cx-set-real-entry sigma 0.01)
			(cx-set-real-entry co 100.0)
		)

		; 定义回调函数，点击按钮OK后调用
		;(define (apply-cb . args)
		;	(display "\nThe dialog closed!\n")
		;) 
		(define apply-cb #f)
		
	  	; 定义回调函数，点击按钮Apply and Execute UDF后调用
		(define (btn-cb1 . args)
			(rpsetvar 'myudf/sigma (cx-show-real-entry sigma)) ; 将GUI中的控件值赋值给变量
			(rpsetvar 'myudf/co (cx-show-real-entry co))  ; 将GUI中的控件值赋值给变量
			(%run-udf-apply 7) ; 调用UDF宏，这里的7为UDF宏中的mode参数
			(ti-menu-load-string "define/user-defined/execute-on-demand \"adjoint_probability::libudf\"") ;执行Define_on_Demand宏				
		)
		
		; 定义回调函数，点击按钮ZHHL后调用
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
					(set! dialog-box (cx-create-panel "Source Inverse: 2022-04-02" apply-cb update-cb))
					(set! table (cx-create-table dialog-box "Some Useful Parms:" 'border #t 'below 0 'right-of 0))
			   
					(set! sigma (cx-create-real-entry table "For sigma:" 'row 1 'col 1))
					(set! co (cx-create-real-entry table "For Co:" 'row 1 'col 2)) 
					(set! btn1 (cx-create-button table "Execute UDF"'activate-callback btn-cb1 'row 2 'col 1))	
					(set! btn2 (cx-create-button table "ZHHL" 'activate-callback btn-cb2 'row 2 'col 2))					
				) 
			)         
			(cx-show-panel dialog-box)
		) 
    ) 
) 
 
;调用函数 
(gui-dialog-box)
