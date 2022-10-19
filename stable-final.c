/* **************************************************************
** �ȶ���� **
**************************************************************
Fluent UDFs ����ģ������ABL����

1. ͨ������Ĳ�������ABL�����Ŀ���
2. ȷ���������expertģʽ�½��м���
3. UDFͨ�����뷽ʽ����
4. ZHHL�޸ĵĵط��Ѿ������˱��(ZHHL)

-------------------------------------------
Owner: Hongliang Zhang <zhhl_email@qq.com>
Check Date: 2022-10-18*/

#include "udf.h"
#include "math.h"

/* ģ�Ͳ��� */
#define Cmu 0.09
#define vonKarman 0.4
#define Ce1 1.21
#define Ce2 1.92
#define sigma_k 1.0
#define sigma_e 1.3
#define sigma_theta 1.0
#define PrTurb 0.85

/* ���ٲ��� */
#define z0 0.0128 /*����ֲڶ�, m*/
#define Cs 0.5 /* �ֲڶȳ��� */
#define uStar 0.320 /* Ħ���ٶ�, uStar = (vonKarman*uRef)/log(zRef/z0); */
#define Lin 330 /* Ī������, ��UDF��Lin�������0 */
#define densOper 1.0 /* �����ܶ�. ԭ��Ϊ1.0919. ���Ϊ1.0 */  /*ZHHL*/


/* ********************** ����ٶ� ********************** */
DEFINE_PROFILE(inlet_V_Stable,t,i)
{
	double x[ND_ND];
	double z;
	double phiM;
	face_t f;

	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0; /*ZHHL*/

		phiM = 1.0 + 5.0*(z/Lin);
		F_PROFILE(f, t, i) = (uStar/vonKarman)*(log(z/z0) + phiM -1.0);
	}
	end_f_loop(f,t)
}

/* ********************** ��� k ********************** */
DEFINE_PROFILE(inlet_k_Stable,t,i)
{
	double x[ND_ND];
	double z;
	double phiE,phiM;
	face_t f;

	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0;

		phiM = 1.0 + 5.0*(z/Lin);
		phiE = phiM-z/Lin;
		F_PROFILE(f,t,i) = (pow((double)uStar,(double)2.0)/sqrt(Cmu))*pow((double)phiE/phiM,(double)0.5);
	}
	end_f_loop(f,t)
}

/* ********************** ��� epsilon ********************** */
DEFINE_PROFILE(inlet_e_Stable,t,i)
{
	double x[ND_ND];
	double z;
	double phiE, phiM;
	face_t f;

	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0;

		phiM = 1.0 + 5.0*(z/Lin);
		phiE = phiM-z/Lin;
		F_PROFILE(f,t,i) = phiE*pow(uStar,3.0)/(vonKarman*z);
	}
	end_f_loop(f,t)
}

/* ************************ k Դ��**********************/
/*��ģ�Ͳ���Ҫ�������ģ��*/

DEFINE_SOURCE(k_source_DTU_Stable,c,t,dS,eqn)
{
	double fSt, phiM, phiE, phiH, CkD, source, Gb, Sk, uStarLocal;
	double x[ND_ND];
	double z, L;
	C_CENTROID(x,c,t);
	z=x[2]+z0; /*ZHHL*/
	L = Lin;

	if (N_ITER > 5) {
		phiM = 1.0 + 5.0*(z/L);
		phiE = phiM-z/L;
		phiH = 1.0 + 5.0*(z/L);
		uStarLocal = pow((double)C_K(c,t),(double)0.5)*pow((double)Cmu,(double)0.25)*pow((double)phiM,(double)0.25)*pow((double)phiE,(double)-0.25);

		fSt = 2.0-(z/L) - 10.0*(z/L)*(1.0-2.0*(z/L) + 10.0*(z/L));
		CkD = pow((double)vonKarman,(double)2.0)/(sigma_k*sqrt(Cmu));
		Gb = -C_MU_T(c,t)*pow((double)sqrt(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2]),(double)2.0)*((z/L)/(sigma_theta))*(phiH/pow((double)phiM,(double)2.0)); /* DTU Formulation */
		Sk = pow((double)uStarLocal,(double)3.0)/(vonKarman*L)*(1.0 - (phiH)/(sigma_theta*phiM) - 0.25*CkD*pow((double)phiM,(double)-3.5)*pow((double)phiE,(double)-1.5)*fSt);
		source = -densOper*Sk + Gb;
	}
	else {
		source = 0.0;
	}

	dS[eqn] = 0.0;
	return source;
}


/* ************************ Epsilon Դ��**********************
Epsilon���ݶȵĺ���, Ϊ�����ݶ���Ҫ��expertģʽ�½��м���(tui��������: solve/set/expert)
DTU����ͨ���ı�Ce3���޸�ģ��, ����ѡ������������Gb=0����ͨ��Դ�����������Gb*/

/* DTU Epsilon ���ڱ��ʽ�޸�Ce3�����Դ�� */
DEFINE_SOURCE(epsilon_source_DTU_Stable,c,t,dS,eqn)
{
	double x[ND_ND];
	double z, L;
	double Gb, C3e, source;
	double phiM, phiH, phiE, fe;
	C_CENTROID(x,c,t);

	z = x[2] +z0;
	L = Lin;

	if (N_ITER > 5) {
		phiM = 1.0 + 5.0*(z/L);
		phiE = phiM-z/L;
		phiH = 1.0 + 5.0*(z/L);
		fe = pow((double)phiM,(double)-2.5)*(2.0*phiM-1.0);
		
		Gb = -C_MU_T(c,t)*pow((double)sqrt(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2]),(double)2.0)*((z/L)/(sigma_theta))*(phiH/pow((double)phiM,(double)2.0)); /* DTU Formulation */
		C3e = sigma_theta/(z/L)*(phiM/phiH)*(Ce1*phiM-Ce2*phiE+(Ce2-Ce1)*pow((double)phiE,(double)-0.5)*fe);
		source = Ce1*C_D(c,t)/C_K(c,t)*C3e*Gb;
	}
	else {
		C3e = 0.0;
		Gb = 0.0;
		source = 0.0;
	}

	dS[eqn] = 0.0;
	return source;
}

/* ************************* ��ʼ��, ���ڰ������������ͷ�ֹ����ݶ�ʱ��ɢP42, �����������Ͳ����� ************************** */

DEFINE_INIT(initStable,d)
{
	cell_t c;
	Thread *t;
	double x[ND_ND];
	double phiM, phiE, phiH, pressure, potenTemp, z, zAMSL, L ;
	thread_loop_c(t,d)
	{
		begin_c_loop_all(c,t)
		{
			C_CENTROID(x,c,t);
			z = x[2] + z0;
			L = Lin;

				phiM = 1.0 + 5.0*(z/L);
				phiE = phiM-z/L;
				phiH = 1.0 + 5.0*(z/L);
				C_U(c,t) = (uStar/vonKarman)*(log(z/z0) +phiM -1.0); /* x velocity */ /*ZHHL*/
				C_V(c,t) = 0.0; /*y velocity */
				C_W(c,t) = 0.0; /* z velocity */
				C_K(c,t) = (pow((double)uStar,(double)2.0)/sqrt(Cmu))*pow((double)phiE/phiM,(double)0.5); /* k */
				C_D(c,t) = phiE*pow((double)uStar,(double)3.0)/(vonKarman*z); /* epsilon */
				C_P(c,t) = 0.0; /*Pressure*/
		}
		end_c_loop_all(c,t)
	}
}


/* ********************** ����ֲڶ� ********************** */

/* ���ʹ���Զ����ABL�������溯��, ����Ҫʹ�øôֲڶ� */
DEFINE_PROFILE(wallRoughness,t,i)
{
	double x[ND_ND];
	face_t f;
	
	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		F_PROFILE(f,t,i) = z0;
	}
	end_f_loop(f,t)
}

/* ���δʹ���Զ����ABL�������溯��, ���ø�������ı���ֲڶ� */
DEFINE_PROFILE(wallRoughnessModified,t,i)
{
	double x[ND_ND];
	face_t f;
	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		F_PROFILE(f,t,i) = 9.793*z0/Cs;
	}
	end_f_loop(f,t)
}

/* ************************* ���溯�� ************************** */

/* Designed around u/uStar = 1/K*log(z/z0) �ο�����: Improved k-e model and wall function formulation for the RANS simulation of ABL flows, Parente et al
�øñ��溯��������Ҫz0*9.73/Cs ���, ABLģ���еĴֲڶȳ���ֱ�ӵ��ڵ�һ����Ԫ�߶�= 2*z0*/

DEFINE_WALL_FUNCTIONS(ABL_logLaw, f, t, c0, t0, wf_ret, yPlus, Emod)
{
	double ustar_ground, E_prime, yPlus_prime, zp, dx_mag, wf_value;
	double mu=C_MU_L(c0,t0);
	double xf[ND_ND];
	double xc[ND_ND];
	double dx[ND_ND];

	F_CENTROID(xf, f, t);
	C_CENTROID(xc, c0,t0);

	dx[0] = xc[0] - xf[0];
	dx[1] = xc[1] - xf[1];
	dx[2] = xc[2] - xf[2];
	dx_mag = NV_MAG(dx);
	zp = dx_mag;

	ustar_ground = pow((double)C_K(c0,t0),(double)0.5)*pow((double)Cmu,(double) 0.25);
	E_prime = (mu/densOper)/(z0*ustar_ground);
	yPlus_prime = (zp+z0)*ustar_ground/(mu/densOper);

	switch (wf_ret)
	{
	case UPLUS_LAM:
		wf_value = yPlus_prime;
		break;
	case UPLUS_TRB:
		wf_value = log(E_prime*yPlus_prime)/vonKarman;
		/*wf_value = log(Emod*yPlus)/vonKarman; Standard Fluent*/
		break;
	case DUPLUS_LAM:
		wf_value = 1.0;
		break;
	case DUPLUS_TRB:
		wf_value = 1.0/(vonKarman*yPlus_prime);
		break;
	case D2UPLUS_TRB:
		wf_value = -1.0/(vonKarman*yPlus_prime*yPlus_prime);
		break;
	default:
		printf("Wall function return value unavailable\n");
	}
	return wf_value;
}