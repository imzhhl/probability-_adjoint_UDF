# FLUENT_UDF

各类个人用FLUENT UDF代码

博士毕业后公开

---

# ✨Source Identification Using Adjoint Probability Method

## 单个污染源正向模拟

**Step 1. 设置第一个污染源的位置**
  ``` c
  /*第一个污染源n1的坐标*/
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

    挂载source源项到UDS-0: direct_source_no1
    
    设置UDS diffusivity: diff
    
    只计算UDS方程（因为单向耦合）


## 单个污染源逆向溯源

**Step 1. 设置三个监测点的位置**

```c
/*第一个检测点n1的坐标*/
#define inverse_no1_x 800
#define inverse_no1_y 620
#define inverse_no1_z 1.0

/*第二个检测点n2的坐标*/
#define inverse_no2_x 670
#define inverse_no2_y 320
#define inverse_no2_z 1.0

/*第三个检测点n3的坐标*/
#define inverse_no3_x 825
#define inverse_no3_y 325
#define inverse_no3_z 1.0
```
**Step 2. 编译加载**
    
    挂载source源项到UDS-1: inverse_source_no1
    挂载source源项到UDS-3: inverse_source_no2
    挂载source源项到UDS-5: inverse_source_no3
    
    设置UDS diffusivity: diff
    
    设置UDS-1的Flux Function: adjoint_flux
    设置UDS-3的Flux Function: adjoint_flux  
    设置UDS-5的Flux Function: adjoint_flux
    

