# -*- coding: utf-8 -*-
"""
author： Hongliang Zhang - WHU
date：   2022-08-14
log: 1. rewrite adapted to pyFluent
     2. couple with 'socket' to change (x,y) location in UDF
     3. I can't finish the external compilation of udf，and won't use this code untill UDF VC++ supporting for 2022R2
"""

import ansys.fluent.core as pyfluent
import numpy as np
import random
import time
import matplotlib.pyplot as plt
#%% 利用pyFluent包连接fluent,该包的使用仅限于Fluent 2022R2以后版本

# 定义工作目录
import_filename = r'F:\ZHHL\TE_Doctor\CASES\case220823\fluent16-uniformmesh-0814'
# UDF_Path = r'F:\ZHHL\TE_Doctor\CASES\case220626\python_fluent\python\udf_source.c'

session = pyfluent.launch_fluent(version = "2d", precision='double',processor_count = 8, show_gui=False)

# 用工作站远程计算，连接现有的窗口
# session = pyfluent.launch_fluent(ip='192.168.31.230', port=63993, start_instance=False)

tui = session.solver.tui
root = session.solver.root

#%% 一些函数的声明

# 计算似然函数值
def PROB(point_1_value, point_2_value, point_3_value, point_1_true_value, point_2_true_value, point_3_true_value, sigma):
    sqrt_2pi=np.power(2*np.pi,0.5)
    coef=1/(sqrt_2pi*sigma)
    powercoef=-1/(2*np.power(sigma,2))
    mypow_1=powercoef*(np.power((point_1_value - point_1_true_value),2))
    mypow_2=powercoef*(np.power((point_2_value - point_2_true_value),2))
    mypow_3=powercoef*(np.power((point_3_value - point_3_true_value),2))
    mypow=mypow_1 + mypow_2 +  mypow_3
    return coef*(np.exp(mypow))

#%% 用于声明和设定一些需要用到的参数

# 存储逐时污染物浓度的数组
point_1_values  = []
point_2_values  = []
point_3_values  = []

# 设定调整区间 ΔX 和 ΔY
delta = 0.5

# 设定标准差
sigma = 0.1

# 接受率计数器
count_accept = 0
count_reject = 0

# 统计（X，Y）随迭代的变化
history_X = []
history_Y = []

#%% 生成监测点数据

# 读入case和data文件
root.file.read(file_type="case", file_name=import_filename)
root.file.read(file_type="data", file_name=import_filename)

# 首先加载污染源UDF
# 不能在此进行编译，因UDF采用C++，所以需要进行外部边界，在此只进行加载
# tui.define.user_defined.compiled_functions('compile', 'libudf', 'yes', 'udf_source.cpp')
tui.define.user_defined.compiled_functions('load' , 'libudf')

# 设定真实污染源的位置 
X_true = 2.3
Y_true = 2.6

#------------------------------------------------------------------------------
#下面进行UDF的数据操作...
with open('X_Y_value.txt', 'w') as file_write:
    file_write.write(str(X_true) + str(Y_true))
# 执行Define_on_demand宏
tui.define.user_defined.execute_on_demand ('"python_udf_socket::libudf"')
#UDF数据操作结束...
# -----------------------------------------------------------------------------

# 初始化patch
tui.solve.patch('air', [], 'uds-0', 'no', '0')

# fluent计算
tui.solve.set.equations('flow', 'no')
tui.solve.set.equations('ke', 'no')
tui.solve.set.equations('uds-0', 'yes')
tui.solve.iterate(10)

# 污染物云图显示
# tui.display.contour(' uds-0-scalar', '0', '0.5')

# 统计真实污染源时检测点处污染物的浓度
point_1_result = root.solution.report_definitions.compute(report_defs=['report-def-1'])
point_1_true_value = point_1_result['report-def-1'][0]
point_2_result = root.solution.report_definitions.compute(report_defs=['report-def-2'])
point_2_true_value = point_2_result['report-def-2'][0]
point_3_result = root.solution.report_definitions.compute(report_defs=['report-def-3'])
point_3_true_value = point_3_result['report-def-3'][0]

# 打印真实污染源时检测点处污染物的浓度
print(f"point_1_true_value = {point_1_true_value}")
print(f"point_2_true_value = {point_2_true_value}")
print(f"point_3_true_value = {point_3_true_value}")

#%% 利用MCMC反演污染源的位置

#随机产生初始点
X = random.uniform(0.0, 9.0)
X = round(X, 1)
Y = random.uniform(0.0, 3.0)
Y = round(Y, 1)

######## 迭代开始，计时开始 ##########
start = time.perf_counter() 
       
for i in range(100):
    
    # 初始化patch
    tui.solve.patch('air', [], 'uds-0', 'no', '0')


    #--------------------------------------------------------------------------
    #下面进行UDF的数据操作...
    with open('X_Y_value.txt', 'w') as file_write:
        file_write.write(str(X) + str(Y))
    # 执行Define_on_demand宏
    tui.define.user_defined.execute_on_demand('"python_udf_socket::libudf"')
    #UDF数据操作结束...
    # -------------------------------------------------------------------------

    # fluent计算
    tui.solve.set.equations('flow', 'no')
    tui.solve.set.equations('ke', 'no')
    tui.solve.set.equations('uds-0', 'yes')
    tui.solve.iterate(10)
    
    # 统计污染物的浓度  
    point_1_result = root.solution.report_definitions.compute(report_defs=['report-def-1'])
    point_1_value  = point_1_result['report-def-1'][0]
    point_1_values.append(point_1_value)    
    point_2_result = root.solution.report_definitions.compute(report_defs=['report-def-2'])
    point_2_value  = point_2_result['report-def-2'][0]
    point_2_values.append(point_2_value)  
    point_3_result = root.solution.report_definitions.compute(report_defs=['report-def-3'])
    point_3_value  = point_3_value  = point_3_result['report-def-3'][0]
    point_3_values.append(point_3_value)

    # 计算似然函数
    Likehood       = PROB(point_1_value, point_2_value, point_3_value, point_1_true_value, point_2_true_value, point_3_true_value, sigma)
    
    # 根据建议分布确定新的污染物的位置(X_new,Y_new)（假设建议分布为均匀分布）
    X_new = random.uniform(X - delta, X + delta)
    X_new = round(X_new, 1)
    Y_new = random.uniform(Y - delta, Y + delta)
    Y_new = round(Y_new, 1)
    # 超X计算边界罚回
    while X_new < 0 or X_new > 9:
        X_new = random.uniform(X - delta, X + delta)
        X_new = round(X_new, 1)
    # 超Y计算边界罚回    
    while Y_new < 0 or Y_new > 3:
        Y_new = random.uniform(Y - delta, Y + delta)
        Y_new = round(Y_new, 1)  

    #--------------------------------------------------------------------------
    #下面进行UDF的数据操作...
    with open('X_Y_value.txt', 'w') as file_write:
        file_write.write(str(X_new) + str(X_new))
    # 执行Define_on_demand宏
    tui.define.user_defined.execute_on_demand('"python_udf_socket::libudf"')
    #UDF数据操作结束...
    # -------------------------------------------------------------------------

    # fluent计算
    tui.solve.set.equations('flow', 'no')
    tui.solve.set.equations('ke', 'no')
    tui.solve.set.equations('uds-0', 'yes')
    tui.solve.iterate(10)
    
    # 统计污染物的浓度  
    point_1_result = root.solution.report_definitions.compute(report_defs=['report-def-1'])
    point_1_value  = point_1_result['report-def-1'][0]
    point_1_values.append(point_1_value)  
    point_2_result = root.solution.report_definitions.compute(report_defs=['report-def-2'])
    point_2_value  = point_2_result['report-def-2'][0]
    point_2_values.append(point_2_value)   
    point_3_result = root.solution.report_definitions.compute(report_defs=['report-def-3'])
    point_3_value  = point_3_value  = point_3_result['report-def-3'][0]
    point_3_values.append(point_3_value)

    # 计算新坐标下的似然函数    
    Likehood_new   = PROB(point_1_value, point_2_value, point_3_value, point_1_true_value, point_2_true_value, point_3_true_value, sigma)

    # 确定是否接受
    u = random.uniform(0, 1)
    temp = min(1, Likehood_new/Likehood)
    # 接受
    if u < temp:
        X = X_new
        Y = Y_new
        count_accept = count_accept + 1
   # 拒绝
    else:
        X = X
        Y = Y
        count_reject = count_reject + 1
        
    history_X.append(X)
    history_Y.append(Y)

    # 打印当面MCMC迭代步数
    print(f"current i = {i}")
    
######## 迭代结束，计时结束 ##########    
end = time.perf_counter()

# 打印运行时间
print("运行时间为", round((end-start)/60, 1), 'mins')

# %% 展示结果

# 打印接受率和反演得到的污染源位置
print(f"接受率 = {count_accept/(count_accept + count_reject)}") 
print(f"X_find = {X_new}")
print(f"Y_find = {Y_new}")
print(f"X_true = {X_true}")
print(f"Y_true = {Y_true}")

# 直方图显示
plt.hist(history_X, bins=30, rwidth=500000, density=True)
plt.hist(history_Y, bins=30, rwidth=500000, density=True)
