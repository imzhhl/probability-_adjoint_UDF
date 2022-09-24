
/*==============================================*/
/* 1. Implementation of  mixed boundary (3st) condition for a UDS				*/
/* 2. Implementation of inverse flux function for UDS										*/
/* 3. Implementation of diffusivity coefficient for UDS										*/
/* 4. the UDS	source term is applied by splitting Mesh in FLUENT					*/
/*==============================================*/

#include "udf.h"
#include "sg.h"

enum
{
	PHI0,	/*PHI0 = 0*/
	PHI1,	/*PHI1 = 1*/
	PHI2,	/*PHI2 = 2*/
	PHI3,	/*PHI3 = 3*/
	PHI4,	/*PHI4 = 4*/
	PHI5,	/*PHI5 = 5*/
	PHI6,	/*PHI6 = 6*/
	PHI7,	/*PHI7 = 7*/
	PHI8,	/*PHI8 = 8*/
	PHI9,	/*PHI9 = 9*/
	N_REQUIRED_UDS	/*N_REQUIRED_UDS = 10*/
};

/* Names of the user-defined scalar to be used */

/*扩散系数*/
#define DIFF 50

/*定义污染源和监测点的大小范围*/
#define di 10

/*定义污染源强度*/
#define intensity 100

/*第一个污染源n1的坐标*/
float direct_no1_x  = 100;
float direct_no1_y =  0;

/*第一个检测点n1的坐标*/
#define inverse_no1_x 200
#define inverse_no1_y 150

/*第二个检测点n2的坐标*/
#define inverse_no2_x 500
#define inverse_no2_y 150

/*第三个检测点n3的坐标*/
#define inverse_no3_x 800
#define inverse_no3_y 150

/*第四个检测点n4的坐标*/
#define inverse_no4_x 200
#define inverse_no4_y 0

/*第五个检测点n5的坐标*/
#define inverse_no5_x 500
#define inverse_no5_y 0

/*第六个检测点n6的坐标*/
#define inverse_no6_x 800
#define inverse_no6_y 0

/*第七个检测点n7的坐标*/
#define inverse_no7_x 200
#define inverse_no7_y -150

/*第八个检测点n8的坐标*/
#define inverse_no8_x 500
#define inverse_no8_y -150

/*第九个检测点n9的坐标*/
#define inverse_no9_x 800
#define inverse_no9_y -150

/*源项扩散系数*/
DEFINE_DIFFUSIVITY(diff,c,t,i)
{
	return DIFF;
}

/*源项扩散系数*/
DEFINE_DIFFUSIVITY(diff_2,c,t,i)
{
	return 1.0 * 2.88e-05 + C_MU_EFF(c,t) / 0.7;
}

DEFINE_DIFFUSIVITY(mean_age_diff, c, t, i)
{
	real Schmidt = 1.0;
	return C_MU_L(c, t)/Schmidt + C_MU_T(c, t)/Schmidt;
}

DEFINE_DIFFUSIVITY(Vf_diffusivity, c, t, i)
{
	return  C_MU_T(c,t)/0.7; 
}



DEFINE_ADJUST(n_required_UDS_adjust, domain)
{
	/* Make sure there are enough user defined-scalars. */
	if (n_uds < N_REQUIRED_UDS)
		Internal_Error("not enough user-defined scalars allocated");
}

DEFINE_PROFILE(inverse_1_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI1);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI1)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI1 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI1),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI1);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_2_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI2);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI2)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI2 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI2),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI2);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_3_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI3);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI3)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI3 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI3),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI3);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_4_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI4);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI4)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI4 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI4),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI4);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_5_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI5);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI5)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI5 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI5),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI5);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_6_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI6);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI6)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI6 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI6),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI6);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_7_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI7);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI7)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI7 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI7),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI7);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_8_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI8);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI8)))) 
	/*NULLP returns TRUE if storage is not allocated,
	checks whether the storage of the gradient of the user defined scalar 
	with index PHI8 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI8),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI8);
	if (temp3<0) temp3=0;
	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

DEFINE_PROFILE(inverse_9_profile, thread, nv)
{
  /* constants must be specified correctly for the mixed BC in the preamble section */
  real Un;

  /* ==================== */

  face_t f;
  
  real A[ND_ND], dG[ND_ND], dr0[ND_ND], es[ND_ND], dr, A_by_es;
  real NV_VEC(psi_vec); /*velocity vector*/
  real A_unit[ND_ND];
  real Af;
  real beta0, gamma;
  real temp1, temp2, temp3;
  
  Thread *t0=thread->t0; 
  
  begin_f_loop(f, thread)
    {

     /* identify the cell thread adjacent to the face thread f */
	cell_t c0 = F_C0(f, thread); 
 
	BOUNDARY_FACE_GEOMETRY(f, thread, A, dr, es, A_by_es, dr0);
	Af=NV_MAG(A); //面法向向量的模
	gamma=C_UDSI_DIFF(c0, t0, PHI9);//扩散系数

	NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
	//NV_D(psi_vec,  =, F_U(f,thread), F_V(f,thread), F_W(f,thread));	
	F_AREA(A_unit, f, thread);		
	NV_S(A_unit ,/=, NV_MAG(A_unit)); /*获得面A的单位矢量*/
	Un = NV_DOT(psi_vec, A_unit);

	if (NULLP(T_STORAGE_R_NV(t0, SV_UDSI_G(PHI9)))) 
	/*NULLP returns TRUE if storage is not allocated,
	Checking whether the storage of the gradient of the user defined scalar 
	with index PHI9 has been allocated*/
		
		beta0=0;  /* if gradient is not allocated and stored yet, bypass
						 the following macro (it happens when saved
						 case/data files are being read during read-case-data */
	else
		BOUNDARY_SECONDARY_GRADIENT_SOURCE(beta0, SV_UDSI_G(PHI9),dG, es, A_by_es, gamma);

	/* temporary variables used in the profile expression */
	temp1=gamma*A_by_es/dr;
	temp2=Un*Af;
    temp3=C_UDSI(c0, t0, PHI9);
	//(temp3<0) temp3=0;

	F_PROFILE(f, thread, nv) = (temp1*temp3- beta0)/(ABS(temp2) + temp1);
    }
  end_f_loop(f, thread)
  
}

/*----------------------------------------------------------对流项取正/反----------------------------------------------------------*/

DEFINE_UDS_FLUX(Direct_flux,f,t,i)
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

		flux =  NV_DOT(psi_vec, A);
	}
	else
	{
		c1 = F_C1(f,t);
		t1 = F_C1_THREAD(f,t); 

		NV_DS(psi_vec,  =, C_U(c0,t0),C_V(c0,t0),C_W(c0,t0),*,C_R(c0,t0));
		NV_DS(psi_vec, +=, C_U(c1,t1),C_V(c1,t1),C_W(c1,t1),*,C_R(c1,t1));

		flux = NV_DOT(psi_vec, A)/2.0;
	}

	return flux;
}

DEFINE_UDS_FLUX(Inverse_flux,f,t,i)
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

		flux =  NV_DOT(psi_vec, A);
	}
	else
	{
		c1 = F_C1(f,t);
		t1 = F_C1_THREAD(f,t); 

		NV_DS(psi_vec,  =, C_U(c0,t0),C_V(c0,t0),C_W(c0,t0),*,C_R(c0,t0));
		NV_DS(psi_vec, +=, C_U(c1,t1),C_V(c1,t1),C_W(c1,t1),*,C_R(c1,t1));

		flux = NV_DOT(psi_vec, A)/2.0;
	}

	return (-1.0)*flux;
}

DEFINE_UDS_FLUX(adjoint_flux,f,t,i)
{
	if(i == 0 || i == 10)
		return F_FLUX(f, t);
	else
		return -1.0*F_FLUX(f, t);
}

DEFINE_ON_DEMAND(value)
{
	real NV_VEC(psi_vec), NV_VEC(A);
	real NV_VEC(i);
	real NV_VEC(j);
	real NV_VEC(k);
	face_t f;
	int ID = 12;
	cell_t c0;
	Thread *t0;
	Domain *domain = Get_Domain(1);
	Thread *f_thread = Lookup_Thread(domain, ID);

	begin_f_loop(f, f_thread)
	{
		c0 = F_C0(f,f_thread);
		t0 = F_C0_THREAD(f,f_thread);
		//		Message0("\nC_U = %0.6f\n", C_U(c0,t0));
		//		Message0("C_V = %0.6f\n", C_V(c0,t0));
		//		Message0("F_U = %0.6f\n", F_U(f,f_thread));
		//		Message0("F_V = %0.6f\n", F_V(f,f_thread));
		//		NV_D(psi_vec,  =, F_U(f,f_thread), F_V(f,f_thread), F_W(f,f_thread));
		NV_D(psi_vec,  =, C_U(c0,t0), C_V(c0,t0), C_W(c0,t0));
		NV_D(i ,=, 1, 0 ,0);
		NV_D(j ,=, 0, 1 ,0);
		NV_D(k ,=, 0, 0 ,1);
		F_AREA(A,f,f_thread);
		NV_S(A ,/=, NV_MAG(A)); //获得面A的单位矢量
		Message0("\nNV_DOT(i, A) = %0.6f\n",NV_DOT(i, A));
		Message0("NV_DOT(j, A) = %0.6f\n",NV_DOT(j, A));
		Message0("NV_DOT(k, A) = %0.6f\n",NV_DOT(k, A));
		//		Message0("\nu_i = %0.6f\n",NV_DOT(psi_vec ,A));
		//		Message0("u_j = %0.6f\n",NV_DOT(psi_vec,A));
		//		Message0("u_k = %0.6f\n",NV_DOT(psi_vec ,A));
		//		Message0("\nC_UDSI -F_UDSI= %0.6f\n",C_UDSI(c0,t0,0)-F_UDSI(f,f_thread,0));

	}
	end_f_loop(f, f_thread)
}
