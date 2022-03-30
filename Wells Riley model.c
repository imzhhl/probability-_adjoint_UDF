/*******************************************************/
/* Wells-Riley model used Euler method same to QianHua**/
/* Author:  ZhangHongliang_WHU**************************/
/* Date:    2021/05/25**********************************/
/* 2021/07/03 debug*************************************/
/* 2021/07/26 and W_R model code************************/
/* 2021/07/27 debug*************************************/
/* 2022/02/23 许航改************************************/
/*******************************************************/
/* NOTE：***********************************************/
/* 1. 编译运行******************************************/
/* 2. fluent当中每个参数都有说明 用户可以根据需要调整****/
/* 3. 设置Number of User-Defined Memory Locations为1****/
/* 4. 此处存储了1组数据，为感染概率**********************/
/* 4. user define memory 0变量为Wells_Riley_cal*********/


#include "udf.h"
#define GX 9.81
#define DP 0.0000025
#define PP 800 
#define PA 1.197
#define NI 0.0000182/*NI is the Dynamic viscosity of airborne,T=22*/
#define SC 0.9
#define TT 298.15
#define MU 0.0000182
#define BV 110.4
#define K 0.1 	/*k is the constant death rate of airborne*/
#define e_t 1 /*exposure time, hour*/


/*version tag*/
DEFINE_PROFILE(version_02_22,thread,i)
	{
		face_t f;
		begin_f_loop(f,thread)
			{
		
			}
		end_f_loop(f,thread)
	}


/*Define the location and strength of the release source*/
DEFINE_SOURCE(souce_begin, c, t, dS, eqn)
	{
		real x[ND_ND];
		real source;
		C_CENTROID(x, c, t);
		if(x[0]>7.3&&x[0]<8.15&&x[1]>1.1&&x[1]<1.3&&x[2]>15&&x[2]<16)
			source=78; 
		else
			source=0.0;
		dS[eqn]=0;
		return source;
	}


/*Define the diffusivity of biological particles in room air*/
DEFINE_DIFFUSIVITY(diffu_coef,c,t,i)
	{
		return C_R(c,t) * 2.88e-05 + C_MU_EFF(c,t) / 0.7;
	}	


/*Define the drift flux and settling velocity of convection term*/
DEFINE_UDS_FLUX(Adjoint_flux,f,t,i)
	{
    	cell_t c0, c1 = -1;
	
    	Thread *t0, *t1 = NULL;

    	real US, MV_1, MV_2, MV_3, mu_g, NV_VEC(psi_vec), NV_VEC(A), flux = 0.0;
		real gravity_vector[]={0,0,-1},area[ND_ND];
	
    	c0 = F_C0(f,t);        /*returns the index of a face’s neighboring c0 cell*/
    	t0 = F_C0_THREAD(f,t); /* Get the Thread id of the c0 cell adjacent to the face*/
    	F_AREA(A, f, t);
	
		MV_1 = TT/288.15;
		MV_2 = pow(MV_1, 1.5);
		MV_3 = (288.15+BV)/(TT+BV);
		mu_g = MV_2*MV_3;

		US=0.545*DP*DP*(PP-PA)/mu_g;

  		if (BOUNDARY_FACE_THREAD_P(t)) /*returns TRUE if Thread *t is a boundary face thread*/
    		{
    			real dens;

      			if (NNULLP(THREAD_STORAGE(t,SV_DENSITY)))/*test whether the memory for a particular field variable has already been allocated on a given Thread or not.*/
        			dens = F_R(f,t);
      			else
        			dens = C_R(c0,t0);

      			NV_DS(psi_vec,  =, F_U(f,t), F_V(f,t), F_W(f,t), *, dens);
      			flux = NV_DOT(psi_vec, A);
   			}
 	 	else
    		{
     			  c1 = F_C1(f,t); /*returns the index of a face’s neighboring c1 cell*/
      			t1 = F_C1_THREAD(f,t); /* Get the Thread id of the c1 cell adjacent to the face*/

      			NV_DS(psi_vec,  =, C_U(c0,t0),C_V(c0,t0),C_W(c0,t0),*,C_R(c0,t0));
      			NV_DS(psi_vec, +=, C_U(c1,t1),C_V(c1,t1),C_W(c1,t1),*,C_R(c1,t1));

      			flux = NV_DOT(psi_vec, A)/2.0;
    		}
	
		return flux+C_R(c0,t0)*US*NV_DOT(gravity_vector,A);
    
	}

/*Define the deposition rate of biological particles on the upward surface*/
DEFINE_PROFILE(upward_boundrycondition,tf,i)
	{
		face_t f;
		cell_t c0;
		real x[ND_ND];
		Thread *t0 = THREAD_T0(tf);
		real bm, bn, farea, NV_VEC(A);
		real taow, ufr; /* friction velocity */
		real rplus;
		real aa, bna, bnb, ba, bb, aresult, bresult;
		real Iresult;
		real vdu;
	
		begin_f_loop(f,tf)
			{
				c0 = F_C0(f,tf);
				F_CENTROID(x,f,tf);
				taow = NV_MAG(F_STORAGE_R_N3V(f,tf, SV_WALL_SHEAR));
				ufr = sqrt (taow/PA);
				rplus = (DP/2)*ufr/((PP-PA)*GX*DP*DP/(18*NI));/* Stokes  gravitational setting velocity */
				aa = log((10.92*1.035744+4.3)*(10.92*1.035744+4.3)*(10.92*1.035744+4.3)/(1.111111+0.0609));
				aresult = 0.5*aa-13.64334;
				bna = 11.31+rplus;
				bnb = 0.0007669*rplus*rplus*rplus;
				ba = log(bna*bna*bna/(11.31+bnb));
				bb = atan((2*rplus-11.31)/19.59001);
				bresult = 0.5*ba + 1.732*bb;
				Iresult = 3.3931*(aresult-bresult)+39;
				vdu = (PP-PA)*GX*DP*DP/(18*NI)/(1-exp((-1)*(PP-PA)*GX*DP*DP/(18*NI)*Iresult/ufr));
				bm = C_UDSI(c0,t0,0);
				farea = F_AREA(A,f,tf);
				bn = vdu*farea*bm;
				F_PROFILE(f,tf,i) = bm - bn;
			}
		end_f_loop(f,tf)
	}


/*Define the deposition rate of biological particles on the downward surface*/	
DEFINE_PROFILE(downward_boundrycondition,tf,i)
	{
		face_t f;
		cell_t c0;
		real x[ND_ND];
		real A[ND_ND];
		real taoarea, area;
		Thread *t, *c_thread;
		Thread *t0 = THREAD_T0(tf);/*cell thread pointer for cell c0*/
		real bm, bn, farea;
		real taow, ufr; /* friction velocity */
		real rplus;
		real aa, bna, bnb, ba, bb, btan, aresult, bresult;
		real Iresult;
		real vdd;
	
	
		begin_f_loop(f,tf)
			{
				c0 = F_C0(f,tf);/* Get the cell id of the cell adjacent to the face*/
				F_CENTROID(x,f,tf);
				F_AREA(A,f,tf);
				taow = NV_MAG(F_STORAGE_R_N3V(f,tf, SV_WALL_SHEAR));
				area = NV_MAG(A);
				taoarea = taow/area;
				ufr = sqrt (taoarea/PA);
				rplus = (DP/2)*ufr/((PP-PA)*GX*DP*DP/(18*NI));/* Stokes  gravitational setting velocity */
				aa = log((10.92*1.035744+4.3)*(10.92*1.035744+4.3)*(10.92*1.035744+4.3)/(1.111111+0.0609));
				aresult = 0.5*aa-13.64334;
				bna = 11.31+rplus;
				bnb = 0.0007669*rplus*rplus*rplus;
				ba = log(bna*bna*bna/(11.31+bnb));
				bb = atan((2*rplus-11.31)/19.59001);
				bresult = 0.5*ba + 1.732*bb;
				Iresult = 3.3931*(aresult-bresult)+39;
				vdd = (PP-PA)*GX*DP*DP/(18*NI)/(exp((-1)*(PP-PA)*GX*DP*DP/(18*NI)*Iresult/ufr)-1);
				bm = C_UDSI(c0,t0,0);
				farea = F_AREA(A,f,tf);
				bn = vdd*bm*farea;
				F_PROFILE(f,tf,i) = bn;
			}
		end_f_loop(f,tf)
	}

/*Define the deposition rate of biological particles on the vertical surface*/
DEFINE_PROFILE(vertical_boundrycondition,tf,i)
	{
		face_t f;
		cell_t c0;
		real x[ND_ND];
		Thread *t0 = THREAD_T0(tf);
		real bm, bn, farea, NV_VEC(A);
		real taow, ufr; /* friction velocity */
		real rplus;
		real aa, bna, bnb, ba, bb, aresult, bresult;
		real Iresult;
		real vdv;
	
		begin_f_loop(f,tf)
			{
				c0 = F_C0(f,tf);
				F_CENTROID(x,f,tf);
				taow = NV_MAG(F_STORAGE_R_N3V(f,tf, SV_WALL_SHEAR));
				ufr = sqrt (taow/PA);
				rplus = (DP/2)*ufr/((PP-PA)*GX*DP*DP/(18*NI));/* Stokes  gravitational setting velocity */
				aa = log((10.92*1.035744+4.3)*(10.92*1.035744+4.3)*(10.92*1.035744+4.3)/(1.111111+0.0609));
				aresult = 0.5*aa-13.64334;
				bna = 11.31+rplus;
				bnb = 0.0007669*rplus*rplus*rplus;
				ba = log(bna*bna*bna/(11.31+bnb));
				bb = atan((2*rplus-11.31)/19.59001);
				bresult = 0.5*ba + 1.732*bb;
				Iresult = 3.3931*(aresult-bresult)+39;
				vdv = ufr/Iresult;
				bm = C_UDSI(c0,t0,0);
				farea = F_AREA(A,f,tf);
				bn = vdv*farea*bm;
				F_PROFILE(f,tf,i) = bm - bn;
			}
		end_f_loop(f,tf)
	}


/* Define W_R calculation code*/
DEFINE_ON_DEMAND(Wells_Riley_cal)
	{
		real P;  				/*The infection risk of one susceptible*/
		real PE; 				/*The escape probability*/
		real N;  				/*the number of microorganisms or concentration of the quanta*/
		real p=0.36;  	/*呼吸通风量6L/min=0.36m3/h*/
			
		Domain *d;               /* declare domain pointer since it is not passed as an argument to the DEFINE macro  */
    Thread *t;               /* 声明Thread指针 */
    cell_t c;                /* 声明cell_t变量，用于循环宏 */
    d = Get_Domain(1);       /* 这个要置于所有定义的变量之后 Get the domain using Fluent utility,单相流取值为1*/	
		thread_loop_c(t, d)
			{
				begin_c_loop(c, t)
					{
						N=(1-K)*C_UDSI(c,t,0);
						PE=exp(-p*N*e_t);
						P =1-PE;
						/*UDM设置----------------------------------------------------------------*/
						/*调用UDM需要首先设置Number of User-Defined Memory Locations，该值设置为1*/
						C_UDMI(c,t,0)=P;												
					}			
				end_c_loop(c, t)
			}
	}
