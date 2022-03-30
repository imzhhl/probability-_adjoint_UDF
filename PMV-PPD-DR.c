/* 使用时将此文档的扩展名改为.c 然后编译或解释到fluent当中 每个参数都有说明 用户可以根据需要调整 调用UDM需要首先设置Number of User-Defined Memory Locations(注意不是存储在Node节点) 设置为3,此处存储了三组数据，分别为PMV,PPD,DR */
/*user define memory 0变量为PMV*/
/*user define memory 1变量为PPD*/
/*user define memory 2变量为DR*/
#include "udf.h"

/* Author:wanghao SOUTH CHINA UNIVERSITY of TECHNOLOGY BEEL 2014*/
/* Author:zhanghongliang DALIAN UNIVERSITY of TECHNOLOGY FIE 2019*/
/*增加计算DR值(吹风感系数)和PPD值代码 2019*/

DEFINE_ON_DEMAND(PMV_PPD_DR_calc)
{
    /*定义real数据类型fluent会转译为real数据类型*/
	#if !RP_HOST
	
    real t_wh;                                  /*空气温度,℃*/
    real speed_u_wh;                            /* u方向速度分量,m/s */
    real speed_v_wh;                            /* v方向速度分量,m/s*/ 
    real speed_w_wh;                            /* w方向速度分量,m/s*/ 
    real speed_wh;                              /* 三个方向的和速度,m/s*/ 
    real re_zhhl;                               /* 雷诺数*/ 
    real rho_zhhl;                              /* 空气密度*/ 
    real mu_zhhl;                               /*动力粘度系数*/
    real diam_zhhl=3.0;                         /*房间特征长度,置换通风取房间高度,根据实际情况修改*/
    real tu_zhhl;                               /* 湍动强度,10%~60%之间,混合通风默认40%,置换通风默认20%*/ 
    real RH_wh=40.0;                            /* 输入相对湿度 水蒸气饱和时输入100 */
    real pa_wh;                                 /*水蒸气分压力 后面通过相对湿度计算水蒸气分压力,pa*/
    real icl_wh=0.155;                          /*服装热阻 在计算时单位不是clo 此处将0.3clo换算成0.05㎡K/W*/
    real fcl_wh;                                /*穿衣服人体外表面与裸体人体表面积之比*/
    real tcl_wh;                                /*服装外表面温度,℃*/
    real tcl1_wh;
    real tcl2_wh;
    real temperary1_wh;
    real temperary2_wh;
    real hc_wh;                                 /*表面换热系数,w/m2K*/
    real tr_wh;                                 /*平均辐射温度 可以后期自己定义,℃*/ 
    real M_wh=69.78;                            /*新陈代谢量,W*/ 
    real W_wh=0.0;                              /*机械做功,W*/ 
    real a_wh;                                  /*heat loss by radiation refer to<EN ISO 7730:2005 page24> */ 
    real b_wh;                                  /*heat loss by convection refer to<EN ISO 7730:2005 page24> */ 
    real c_wh;                                  /*heat loss diff.through skin refer to<EN ISO 7730:2005 page24> */ 
    real d_wh;                                  /*heat loss by sweating(comfort) refer to<EN ISO 7730:2005 page24> */ 
    real e_wh;                                  /*latent respiration heat loss refer to<EN ISO 7730:2005 page24> */ 
    real f_wh;                                  /*dry respiration heat loss refer to<EN ISO 7730:2005 page24> */         
    real L_wh;
    real PMV_wh;                                /*计算所得PMV值*/
    real PPD_zhhl;                              /*计算所得PPD值*/
    real DR_zhhl;                               /*计算所得DR(吹风感系数）值*/
	
	Domain *d;                                  /* declare domain pointer since it is not passed as an argument to the DEFINE macro  */
    Thread *t;                                  /* 声明Thread指针 */
    cell_t c;                                   /* 声明cell_t变量，用于循环宏 */
  
    d = Get_Domain(1);                          /*!!!!!这个要置于所有定义的变量之后 Get the domain using Fluent utility,单相流取值为1*/

    thread_loop_c(t,d)                          /* Loop over all cell threads in the domain 因为loop是针对每个cell的所以针对cell的计算都要包含到loop里面*/
    {
        begin_c_loop(c,t)
        {
            /*计算PMV值--------------------------------------------------------------------------------------------*/
            t_wh=C_T(c,t)-273.15;               /*获取空气温度 ℃*/    
            speed_u_wh=C_U(c,t);                /*获取u方向速度分量 m/s */
            speed_v_wh=C_V(c,t);                /*获取v方向速度分量 m/s */
            speed_w_wh=C_W(c,t);                /*获取w方向速度分量 m/s */
            tr_wh=t_wh;                         /*此处把平均辐射温度近似处理为空气温度 ℃ 可以后期自己定义*/

            pa_wh=RH_wh*10.0*exp(16.6536-4030.183/(t_wh+235.0));                                /*计算水蒸气分压力*/
            speed_wh=sqrt(pow(speed_u_wh,2.0)+pow(speed_v_wh,2.0)+pow(speed_w_wh,2.0));         /*计算绝对空气相对速度*/

            if (icl_wh<0.078)
                fcl_wh=1.0+1.290*icl_wh;
            else
                fcl_wh=1.05+0.645*icl_wh;                                                       /*判断穿衣服人体外表面与裸体人体表面积之比*/                                 
 
            tcl1_wh=t_wh+273.0+(35.5-t_wh)/(3.5*icl_wh+0.1);                                    /*计算服装表面温度的初始值*/
            tcl2_wh=(308.7-0.028*(M_wh-W_wh)+icl_wh*fcl_wh*hc_wh*(t_wh+273.0)+icl_wh*0.0000000396*fcl_wh*pow((tr_wh+273.0),4.0)-icl_wh*0.0000000396*fcl_wh*pow(tcl1_wh,4.0))/(1.0+icl_wh*fcl_wh*hc_wh);/*计算中服装表面温度 K */
  
            temperary1_wh=2.38*pow(fabs(tcl2_wh-t_wh-273.0),0.25);                              /*自然换热系数*/ 
            temperary2_wh=12.1*pow(speed_wh,0.5);                                               /*强迫换热系数*/     
            if (temperary1_wh<temperary2_wh)
                hc_wh=temperary2_wh;
            else
                hc_wh=temperary1_wh;                                                            /*判断表面换热系数 并赋值*/  
 


            while (fabs(tcl1_wh-tcl2_wh)>0.001)                                                 /*通过判断差的绝对值进行循环逼近真实值 fabs表示取浮点数据类型的绝对值*/                 
            {
                tcl1_wh=tcl2_wh;                                                                /*循环过程中把上次的计算结果赋给tcl1_wh做下一次计算*/
  
                temperary1_wh=2.38*pow(fabs(tcl1_wh-t_wh-273.0),0.25);                            /*自然换热系数*/ 
                temperary2_wh=12.1*pow(speed_wh,0.5);                                           /*强迫换热系数*/     
                if (temperary1_wh<temperary2_wh)
                    hc_wh=temperary2_wh;
                else
                    hc_wh=temperary1_wh;                                                        /*判断表面换热系数 并赋值*/  
                // tcl2_wh=35.7-0.028*(M_wh-W_wh)-icl_wh*(0.0000000396*fcl_wh*(pow((tcl1_wh+273),4)-pow((tr_wh+273),4))+fcl_wh*hc_wh*(tcl1_wh-t_wh));
                tcl2_wh=(308.7-0.028*(M_wh-0.0)+icl_wh*fcl_wh*hc_wh*(t_wh+273.0)+icl_wh*0.0000000396*fcl_wh*pow((tr_wh+273.0),4.0)-icl_wh*0.0000000396*fcl_wh*pow(tcl1_wh,4.0))/(1.0+icl_wh*fcl_wh*hc_wh);/*计算中服装表面温度 K */
            }

            tcl_wh=tcl2_wh-273;                                                               /*最终的服装表面温度 ℃ */

            a_wh=0.0000000396*fcl_wh*(pow((tcl_wh+273.0),4.0)-pow((tr_wh+273),4.0));
            b_wh=fcl_wh*hc_wh*(tcl_wh-t_wh);
            c_wh=0.00305*(5733.0-6.99*(M_wh-W_wh)-pa_wh);
            if((M_wh-W_wh-58.15)>0.0) d_wh=0.42*(M_wh-W_wh-58.15); else d_wh=0.0;               /*heat loss by sweating(comfort) refer to<EN ISO 7730:2005 page24> */
            e_wh=0.000017*M_wh*(5867.0-pa_wh);
            f_wh=0.0014*M_wh*(34.0-t_wh);            
            L_wh=M_wh-W_wh-(a_wh+b_wh+c_wh+d_wh+e_wh+f_wh);

            PMV_wh=(0.303*exp(-0.036*M_wh)+0.028)*L_wh;
            
            /*计算PPD值--------------------------------------------------------------------------------------------*/
            PPD_zhhl=100.0-95.0*exp(-0.03353*pow(PMV_wh,4.0)-0.2179*pow(PMV_wh,2.0));

            /*计算DR值---------------------------------------------------------------------------------------------*/
            if(speed_wh<0.05) speed_wh=0.05;
            rho_zhhl=C_R(c,t);
            mu_zhhl=C_MU_L(c,t);
            re_zhhl=rho_zhhl*speed_wh*diam_zhhl/mu_zhhl;                                        /*计算雷诺数*/
            tu_zhhl=100.0*0.16*pow(re_zhhl,-0.125);
            //tu_zhhl=20.0;
            DR_zhhl=(34.0-t_wh)*pow((speed_wh-0.05),0.62)*(0.37*speed_wh*tu_zhhl+3.14);
            if(DR_zhhl>100.0) DR_zhhl=100.0;


            /*UDM设置---------------------------------------------------------------------------------------------*/
            /*调用UDM需要首先设置Number of User-Defined Memory Locations，该值设置为3*/

	        C_UDMI(c,t,0)=PMV_wh;
	        C_UDMI(c,t,1)=PPD_zhhl;  
	        C_UDMI(c,t,2)=DR_zhhl;

        }
        end_c_loop(c,t)
    }
	#endif
}
