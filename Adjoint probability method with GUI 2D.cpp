/************************************************************************************/
/*
User-Defined Function for source identification using Adjoint probability method
Author: Hongliang Zhang_WHU                                                                   
Date:   2021-09-07
Type:	2D
*/

/*UDS allocate structure---------------------------------------------------------------------*/
/*
UDS-0----------正向模拟结果（单个污染源）
UDS-1----------测点1的逆向模拟结果（单个污染源）
UDS-2----------测点1的UDS-1的概率计算结果
UDS-3----------测点2的逆向模拟结果（单个污染源）
UDS-4----------测点2的UDS-2的概率计算结果
UDS-5----------测点2的逆向模拟结果（单个污染源）
UDS-6----------测点2的UDS-3的概率计算结果
UDS-7----------测点1，测点2，测点3的联合概率
*/

/*LOG:---------------------------------------------------------------------------------------------*/
/*
0.   2021.09.07-将原先的伴随概率法用UDF实现，代码用VC++ Fluent UDF辅助编写（http://VcUdfStudio.bitbucket.io）
1.   2021.09.07-正向计算结果存于UDS-0，no1逆向计算的结果存于UDS-1，概率计算结果存于UDS-2
2.   2021.09.07-Only used for Fluent >=14.0
3.   2021.09.07-查找多个点的cell
4.   2021.09.08-no2逆向计算的结果存于UDS-3，概率计算结果存于UDS-4
5.   2021.09.08-差不多可以用了
6.   2021.12.07-重新测试
7.   2021.12.08-改成三个测点反演
8.   2021-12.09-全部完成
9.   2022-01-21-改为三维. V2.0
10. 2022-04-02-增加GUI控件
/************************************************************************************/
#include "udf.h"

extern "C"
{
#if RampantReleaseMajor>=13
#include "cxndsearch.h" /*需要的头文件*/
#endif
};

/*扩散系数*/
#define DIFF 10

/*定义污染源和监测点的大小范围*/
#define	 di 5

/*第一个污染源n1的坐标*/
#define direct_no1_x 250
#define direct_no1_y 425

/*第二个污染源n2的坐标*/
// #define direct_no2_x 0
// #define direct_no2_y 0

// /*第三个污染源n3的坐标*/
// #define direct_no3_x 0
// #define direct_no3_y 0

/*第一个检测点n1的坐标*/
#define inverse_no1_x 350
#define inverse_no1_y 420

/*第二个检测点n2的坐标*/
#define inverse_no2_x 670
#define inverse_no2_y 320

/*第三个检测点n3的坐标*/
#define inverse_no3_x 730
#define inverse_no3_y 520

/*随机误差，数据分布的标准差*/
float sigma = 0.01;

/*放大系数*/
float co = 100.0;

/*----------------------------------------------------------------------------------------------------------------------------------------*/
/*当点击GUI界面中的apply按钮时运行此宏*/
DEFINE_EXECUTE_FROM_GUI(check, libudf, mode)
{
	/*这里的mode对应于scm中的回调函数中的相应变量*/
	if(mode == 7)
	{
		Message0("\n-------------------------------------\n");
		sigma = RP_Get_Real("myudf/sigma");
		Message0("Sigma is '%f'\n", sigma);

		co = RP_Get_Real("myudf/co");
		Message0("Co is '%f'\n", co);

		Message0("The UDF has been called successfully!\n");
		Message0("\n"); 
	}

	else
	{
		Message0("Error!\n");
	}
}

/*源项扩散系数*/
DEFINE_DIFFUSIVITY(diff,c,t,i)
{
	return DIFF;
}

/*真实污染源释放source宏1*/
DEFINE_SOURCE(direct_source_no1, c, t, dS, eqn)
{
	real x[ND_ND];
	real source;
	C_CENTROID(x,c,t);

	if(x[0]>direct_no1_x - di&&x[0]<direct_no1_x + di&&x[1]>direct_no1_y - di &&x[1]<direct_no1_y + di )
	{
		source=100;
	}
	else
	{
		source=0;
	}
	dS[eqn] = 0;
	return source;
}

/*----------------------------------------------------------------------------------------------------------------------------------------*/

/*名义污染源释放source宏1*/
DEFINE_SOURCE(inverse_source_no1, c, t, dS, eqn)
{
	real x[ND_ND];
	real source;
	C_CENTROID(x,c,t);

	if(x[0]>inverse_no1_x - di&&x[0]<inverse_no1_x + di&&x[1]>inverse_no1_y - di &&x[1]<inverse_no1_y + di)
	{
		source=1;
	}
	else
	{
		source=0;
	}
	dS[eqn] = 0;
	return source;
}

/*名义污染源释放source宏2*/
DEFINE_SOURCE(inverse_source_no2, c, t, dS, eqn)
{
	real x[ND_ND];
	real source;
	C_CENTROID(x,c,t);

	if(x[0]>inverse_no2_x - di&&x[0]<inverse_no2_x + di&&x[1]>inverse_no2_y - di &&x[1]<inverse_no2_y + di)
	{
		source=1;
	}
	else
	{
		source=0;
	}
	dS[eqn] = 0;
	return source;
}

/*名义污染源释放source宏3*/
DEFINE_SOURCE(inverse_source_no3, c, t, dS, eqn)
{
	real x[ND_ND];
	real source;
	C_CENTROID(x,c,t);

	if(x[0]>inverse_no3_x - di&&x[0]<inverse_no3_x + di&&x[1]>inverse_no3_y - di &&x[1]<inverse_no3_y + di)
	{
		source=1;
	}
	else
	{
		source=0;
	}
	dS[eqn] = 0;
	return source;
}

/*----------------------------------------------------------对流项取反-----------------------------------------------*/
DEFINE_UDS_FLUX(Adjoint_flux,f,t,i)
{
	cell_t c0, c1 = -1;
	Thread *t0, *t1 = NULL;

	real NV_VEC(psi_vec), NV_VEC(A), flux = 0.0;

	c0 = F_C0(f,t);
	t0 = F_C0_THREAD(f,t);
	F_AREA(A, f, t);

	if (BOUNDARY_FACE_THREAD_P(t))
	{
		real dens;

		if (NNULLP(THREAD_STORAGE(t,SV_DENSITY)))
			dens = F_R(f,t);
		else
			dens = C_R(c0,t0);

		NV_DS(psi_vec,  =, F_U(f,t), F_V(f,t), F_W(f,t), *, dens);

		flux = -1 * NV_DOT(psi_vec, A);
	}
	else
	{
		c1 = F_C1(f,t);
		t1 = F_C1_THREAD(f,t); 

		NV_DS(psi_vec,  =, C_U(c0,t0),C_V(c0,t0),C_W(c0,t0),*,C_R(c0,t0));
		NV_DS(psi_vec, +=, C_U(c1,t1),C_V(c1,t1),C_W(c1,t1),*,C_R(c1,t1));

		flux = -1 * NV_DOT(psi_vec, A)/2.0;
	}

	return flux;

}
/*--------------------------------------------------------------------------------------------------------------------------------*/

DEFINE_ON_DEMAND(adjoint_probability)
{
	/* Make sure there are enough user-defined scalars. */
	if (n_uds < 5)
		Internal_Error("NOT enough UDS allocated!!!");

	/*-----------------------------------------根据坐标获取no1、no2坐标的cell------begin----------------------------------------------*/
	real yp1, yp2, yp3; /*不同位置的正向测量值*/
	double PI = 3.1415926;
	real ratio = co; /*放大系数*/
	double sigma_temp = sigma;
	double TEMP_0_a = 0.0, TEMP_2_a = 0.0, TEMP_4_a = 0.0, TEMP_6_a = 0.0;
	double TEMP_0_b = 0.0, TEMP_2_b = 0.0, TEMP_4_b = 0.0, TEMP_6_b = 0.0;
	real cross_UDSI_max = 0.0;

	cell_t c, c1, c2, c3; 
	Thread *t, *t1, *t2, *t3; 
	real coord_Cell1[ND_ND];					/*找到的cell1的中心坐标 */
	real coord_Cell2[ND_ND];					/*找到的cell2的中心坐标 */
	real coord_Cell3[ND_ND];					/*找到的cell3的中心坐标 */
	real coord_Source[ND_ND];				/*找到的污染源的中心坐标 */
	CX_Cell_Id* cx_cell1 = NULL; 
	CX_Cell_Id* cx_cell2 = NULL; 
	CX_Cell_Id* cx_cell3 = NULL; 
	real Pt_to_find1[ND_ND]={inverse_no1_x, inverse_no1_y};		/*指定寻找点1的坐标*/
	real Pt_to_find2[ND_ND]={inverse_no2_x, inverse_no2_y};		/*指定寻找点2的坐标*/
	real Pt_to_find3[ND_ND]={inverse_no3_x, inverse_no3_y};		/*指定寻找点3的坐标*/
	double Pt1[ND_ND], Pt2[ND_ND],  Pt3[ND_ND]; 
	ND_Search *domain_table = NULL; 
	domain_table = CX_Start_ND_Point_Search(domain_table,TRUE,-1); /* 准备开始查找点 */
	NV_V(Pt1, =, Pt_to_find1); 
	NV_V(Pt2, =, Pt_to_find2); 
	NV_V(Pt3, =, Pt_to_find3); 
	cx_cell1 = CX_Find_Cell_With_Point(domain_table, Pt1, 0); 
	cx_cell2 = CX_Find_Cell_With_Point(domain_table, Pt2, 0); 
	cx_cell3 = CX_Find_Cell_With_Point(domain_table, Pt3, 0); 

	Domain *domain = Get_Domain(1);

	if (NULL != cx_cell1 || NULL != cx_cell2 ) 
	{ 
		c1=RP_CELL(cx_cell1); /*找到的cell1序号 */
		c2=RP_CELL(cx_cell2); /*找到的cell2序号 */
		c3=RP_CELL(cx_cell3); /*找到的cell3序号 */
		t1 = RP_THREAD(cx_cell1); /*找到的cell1线索*/
		t2 = RP_THREAD(cx_cell2); /*找到的cell2线索*/
		t3 = RP_THREAD(cx_cell3); /*找到的cell3线索*/
		C_CENTROID(coord_Cell1,c1,t1); /*获取cell1中心坐标 */
		C_CENTROID(coord_Cell2,c2,t2); /*获取cell2中心坐标 */
		C_CENTROID(coord_Cell3,c3,t3); /*获取cell3中心坐标 */
		yp1 = C_UDSI(c1, t1, 0); /*no1处的正向测量值*/
		yp2 = C_UDSI(c2, t2, 0); /*no2处的正向测量值*/
		yp3 = C_UDSI(c3, t3, 0); /*no3处的正向测量值*/

		Message0("\n"); 
		Message0("\n"); 
		Message0("Coordinate of the specified point: (%g, %g, %g)\n",Pt_to_find1[0],Pt_to_find1[1]); 
		Message0("Coordinate of the cell found: (%g, %g, %g)\n",coord_Cell1[0],coord_Cell1[1]);
		Message0("Value of NO1 = %g\n\n",yp1);

		Message0("Coordinate of the specified point: (%g, %g, %g)\n",Pt_to_find2[0],Pt_to_find2[1]); 
		Message0("Coordinate of the cell found: (%g, %g, %g)\n",coord_Cell2[0],coord_Cell2[1]);
		Message0("Value of NO2 = %g\n\n",yp2);

		Message0("Coordinate of the specified point: (%g, %g, %g)\n",Pt_to_find3[0],Pt_to_find3[1]); 
		Message0("Coordinate of the cell found: (%g, %g, %g)\n",coord_Cell3[0],coord_Cell3[1]);
		Message0("Value of NO3 = %g\n\n",yp3);
	}
	else
	{
		Message("Could not find cell at [%g,%g,%g]!\n",Pt_to_find1[0],Pt_to_find1[1]);
		Message("Could not find cell at [%g,%g,%g]!\n",Pt_to_find2[0],Pt_to_find2[1]);
		Message("Could not find cell at [%g,%g,%g]!\n",Pt_to_find3[0],Pt_to_find3[1]);
	}
	domain_table = CX_End_ND_Point_Search(domain_table); //结束查找 
	/*-----------------------------------------根据坐标获取no1、no2坐标的cell------end------------------------------------------------*/

	//for(sigma_temp = 1.0; sigma_temp <= 10.0; sigma_temp++)
	{
		//		for (ratio = 1.0; ratio <= 500.0;)
		{

			thread_loop_c(t, domain)
			{
				begin_c_loop(c, t)
				{	
					TEMP_2_b = TEMP_2_b + (C_UDSI(c, t, 1) * 1.0 / (sigma * sqrt(2.0 * PI))) * exp(-(( ratio * C_UDSI(c, t, 1) - yp1 ) * ( ratio * C_UDSI(c, t, 1) - yp1 )) / (2.0 * sigma * sigma));
					TEMP_4_b = TEMP_4_b + (C_UDSI(c, t, 3) * 1.0 / (sigma * sqrt(2.0 * PI))) * exp(-(( ratio * C_UDSI(c, t, 3) - yp2 ) * ( ratio * C_UDSI(c, t, 3) - yp2 )) / (2.0 * sigma * sigma));
					TEMP_6_b = TEMP_6_b + (C_UDSI(c, t, 5) * 1.0 / (sigma * sqrt(2.0 * PI))) * exp(-(( ratio * C_UDSI(c, t, 5) - yp3 ) * ( ratio * C_UDSI(c, t, 5) - yp3 )) / (2.0 * sigma * sigma));
				}
				end_c_loop(c,t)
			}

			thread_loop_c(t, domain)
			{
				begin_c_loop(c, t)
				{	

					TEMP_2_a = (C_UDSI(c, t, 1) * 1.0 / (sigma * sqrt(2.0 * PI))) * exp(-(( ratio * C_UDSI(c, t, 1) - yp1 ) * ( ratio * C_UDSI(c, t, 1) - yp1 )) / (2.0 * sigma * sigma));
					C_UDSI(c, t, 2) = 100.0 * TEMP_2_a / TEMP_2_b; /*概率大小, % */

					TEMP_4_a = (C_UDSI(c, t, 3) * 1.0 / (sigma * sqrt(2.0 * PI))) * exp(-(( ratio * C_UDSI(c, t, 3) - yp2 ) * ( ratio * C_UDSI(c, t, 3) - yp2 )) / (2.0 * sigma * sigma));
					C_UDSI(c, t, 4) = 100.0 * TEMP_4_a / TEMP_4_b; /*概率大小, % */

					TEMP_6_a = (C_UDSI(c, t, 5) * 1.0 / (sigma * sqrt(2.0 * PI))) * exp(-(( ratio * C_UDSI(c, t, 5) - yp3 ) * ( ratio * C_UDSI(c, t, 5) - yp3 )) / (2.0 * sigma * sigma));
					C_UDSI(c, t, 6) = 100.0 * TEMP_6_a / TEMP_6_b; /*概率大小, % */

					C_UDSI(c, t, 7) = C_UDSI(c, t, 2) + C_UDSI(c, t, 4) + C_UDSI(c, t, 6);

					if(C_UDSI(c, t, 2) * C_UDSI(c, t, 4) * C_UDSI(c, t, 6)>= 0.1 || C_UDSI(c, t, 2) * C_UDSI(c, t, 4)>= 0.1 || C_UDSI(c, t, 4) * C_UDSI(c, t, 6)>= 1 || C_UDSI(c, t, 2) * C_UDSI(c, t, 6)>= 0.1 )
					{
						C_CENTROID(coord_Source, c, t); /*获取cell中心坐标 */
						Message0("Coordinate of the origial source: (%d, %d)\n",direct_no1_x, direct_no1_y); /*TUI输出污染源坐标 */
						Message0("Coordinate of the find source: (%g, %g)\n",coord_Source[0], coord_Source[1]); /*TUI输出反演源坐标 */
						Message0("Successful! \n");
						goto here; /*跳出所有循环 */
					}

				}
				end_c_loop(c, t)
			}			
			ratio = ratio+10.0;
		} 
	}
here:;
	Message0("\nDone! \n");

}
