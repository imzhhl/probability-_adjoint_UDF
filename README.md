# FLUENT_UDF

各类个人用FLUENT UDF代码

博士毕业后公开

---

# ✨Source Identification Using Adjoint Probability Method

## 单个污染源正向模拟

**Step 1. 设置第一个污染源的位置**
  ``` c
  #define direct_no1_x 162.6
  #define direct_no1_y 0.195
  #define direct_no1_z 166.5
  ```
**Step 2. 设置污染源的大小**
  ``` c
  #define	 di 0.6
  ```
**Step 3. FLUENT预先设置**

    预开8个UDS（根据实际需求决定）
  
    UDS-0----------正向模拟结果（单个污染源）
    UDS-1----------测点1的逆向模拟结果（单个污染源）
    UDS-2----------测点1的UDS-1的概率计算结果
    UDS-3----------测点2的逆向模拟结果（单个污染源）
    UDS-4----------测点2的UDS-3的概率计算结果
    UDS-5----------测点3的逆向模拟结果（单个污染源）
    UDS-6----------测点3的UDS-3的概率计算结果
    UDS-7----------测点1、测点2、测点3的联合概率

**Step 4. 编译加载**

    挂载source源项到UDS-0：direct_source_no1
    
    设置UDS diffusivity
    



## 逆向溯源

