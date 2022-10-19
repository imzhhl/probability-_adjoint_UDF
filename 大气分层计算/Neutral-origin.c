#include "udf.h"
#include "math.h"

/* **************************************************************
** Neutral **
**************************************************************
Fluent UDFs for simulating neutral ABL flow

Control via the defined parameters
Ensure the solver is in expert mode
Use compiled UDF method

Model Axis. xz = inlet/outlet plane, yz = sides, z = AGL, origin at the inlet, positive in direction of flow and AGL

C_UDMI - 3 User memory slots, 1User scalar slot
Wall Distance
Cor x
Cor y

C_UDSI
wallPhi - See description in define cell wall distance

-------------------------------------------
Owner: Hendri Breedt <u10028422@tuks.co.za>
Date: 09/11/2017
Version: 00 - Public release */

/* Model Constants - DTU */
#define Cmu 0.03
#define vonKarman 0.4
#define Ce1 1.21
#define Ce2 1.92
#define sigma_k 1.0
#define sigma_e 1.3
#define sigma_theta 1.0
#define PrTurb 0.85

/* Model Constants AM */
#define CmuAM 0.033
#define vonKarmanAM 0.42
#define Ce1AM 1.176
/* The rest are the same as the DTU model */

/* Wind speed relations */
#define z0 0.03 /*m*/
#define Cs 0.5 /* Roughness Constant */
#define uStar 0.1439 /* uStar = (vonKarman*uRef)/log(zRef/z0); */
#define ablHeight 1000.0 /* Height of ABL, this is the height for fixed values of all profiles and sources */

/* Site */
#define globalLat -33.0 /* Latitude of the origin in degrees - This is a dummy value for confidentiality*/
#define siteElevation 0.0 /* Altitude of site AMSL - If you specify the operating pressure from site data then DO NOT change this value. */
#define earthRot 0.000072921159 /* Earth rotational speed */
#define offset 477.0 /* Use to control the z value, this is deducted from the mesh z coordinate. This is the height AGL of the inlet location of the mesh */
#define offsetY -3000.0 /* This is deducted from the local lattitude in the corliolis calculation */

/* General */
#define pi 3.141592
#define g -9.80665
#define R 8.3144598 /* Universal Gas Constant - Dry Air */
#define M 0.0289644 /* Molar mass of Earth's air */
#define Lb -0.0065 /* Standard temperature lapse rate */

/* Operating Conditions - Material Air */
#define presOper 101325 /* Operating Pressure Pa - Internal Solver Pressure. This is the pressure specified at 0m and for this you can use lowest mast pressure reading */
#define tempOper 288.16 /* Operating Temperature - Internal Solver Standard Tempearture. This is the temperature based from the lowest measurement height on the mast. But can be left as the standard value */
#define densOper 1.0919 /* Problem density */
#define Cp 1006.43
#define beta 0.032
#define viscosity 1.7894e-05

/* Initilization */
/* Due to HAGL variations and Fluent not being able to compute cell distance before initialiazing we have to manually set the initialiaze values. These are used for z values lower than maxZInit, afterwards it returns to the inlt profile values */
#define maxZInit 1000.0 /* Height before using init values from inlet profiles */
#define initVelocity 10.0 /* y velocity */
#define initK 2.0 /* k */
#define initEpsilon 2.0 /* epsilon */

/* ********************** Inlet Velocity ********************** */
DEFINE_PROFILE(inletVelocityNeutral, t, i)
{
real x[ND_ND];
real z;
face_t f;

begin_f_loop(f, t)
{
F_CENTROID(x,f,t);
z = x[2] + z0 - offset;
if (z > ablHeight){
z = ablHeight;
}
F_PROFILE(f, t, i) = (uStar/vonKarman)*log(z/z0);
}
end_f_loop(f, t)
}

/* ********************** Inlet k ********************** */
DEFINE_PROFILE(inlet_k_Neutral, t, i)
{
real x[ND_ND];
face_t f;

begin_f_loop(f, t)
{
F_CENTROID(x,f,t);
F_PROFILE(f, t, i) = pow(uStar,2.0)/sqrt(Cmu);
}
end_f_loop(f, t)
}

/* ********************** Inlet epsilon ********************** */
DEFINE_PROFILE(inlet_e_Neutral, t, i)
{
real x[ND_ND];
real z;
face_t f;

begin_f_loop(f, t)
{
F_CENTROID(x,f,t);
z = x[2] + z0 - offset;
if (z > ablHeight){
z = ablHeight;
}
F_PROFILE(f, t, i) = pow(uStar,3.0)/(vonKarman*z);
}
end_f_loop(f, t)
}


/* ********************** Wall Roughness ********************** */

/* Use this if you are using the ABL log law wall function */
DEFINE_PROFILE(wallRoughness,t,i)
{
real x[ND_ND];
face_t f;
begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
F_PROFILE(f,t,i) = z0; /* Use this if you are using the ABL log law wall function */
}
end_f_loop(f,t)
}

/* Modified wall roughness */
DEFINE_PROFILE(wallRoughnessModified,t,i)
{
real x[ND_ND];
face_t f;
begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
F_PROFILE(f,t,i) = 9.793*z0/Cs;
}
end_f_loop(f,t)
}

/* ************************* Cell Wall Distance ************************** */
/* To Use: Define a UDS with Flux Function = none, no Inlet Diffusion
Add Material Property "UDS Diffusivity"; defined-per-uds: constant, Coefficient = 1 [kg/ms]
Add Source Terms for User Scalars in the cell zone: Source Term = 1
Set Boundary Conditions for User Scalar: Specified Value = 0 on all boundaries to which the distance should be computed (boundary lower in the attached sample case); Specified Flux = 0 on all other boundaries.
Define a User-Defined Memory Location in which the UDF stores the computed distance
Hook to Fluent */

DEFINE_EXECUTE_AT_END(computeSelectedWallDistance)
{
Domain *d=Get_Domain(1);
Thread *t;
cell_t c;
real wallPhi, gradWallPhi, wallDistance;

/* Check if UDM and UDS exist */
if (N_UDM < 3 || N_UDS < 1) {
Message0("\n Error: No UDM or no UDS defined! Abort UDF execution.\n");
return;
}

/* Loop over all threads and cells to compute the wall distance */
thread_loop_c(t,d)
{
begin_c_loop(c,t)
{
/* Retrieve wallPhi from UDS-0 */
wallPhi = C_UDSI(c,t,0);
/* Compute magnitude of gradient of wallPhi */
gradWallPhi = NV_MAG(C_UDSI_G(c,t,0));
/* Compute local wall distance */
wallDistance = -gradWallPhi + sqrt(MAX(gradWallPhi*gradWallPhi + 2*wallPhi, 0));

/* Store local wall distance in UDM-0 */
C_UDMI(c,t,0) = wallDistance; /* Call C_UDMI(c,t,0) to retrieve the wall distance */
}
end_c_loop(c,t)
}
}

/* ****************** Sources *********************
**************** Corliolis Force *************** */
DEFINE_SOURCE(Coriolis_X_source,c,t,dS,eqn)
{
real x[ND_ND];
real source;
real Lat, density;

C_CENTROID(x,c,t);

Lat = globalLat + (x[1] - offsetY)*9.0066*1e-6; /* Add the local lattitude change converted from m to degrees */
density = C_R(c,t);

source = 2.0*earthRot*sin(Lat * 3.1459/180)*density*C_V(c,t);
dS[eqn] = 0.0;
C_UDMI(c,t,1) = source;
return source;
}

DEFINE_SOURCE(Coriolis_Y_source,c,t,dS,eqn)
{
real x[ND_ND];
real source;
real Lat, density;

C_CENTROID(x,c,t);

Lat = globalLat + (x[1] - offsetY)*9.0066*1e-6; /* Add the local lattitude change converted from m to degrees */
density = C_R(c,t);

source = -2.0*earthRot*sin(Lat * 3.1459/180)*density*C_U(c,t);
dS[eqn] = 0.0;
C_UDMI(c,t,2) = source;
return source;
}

/* ************************* Initilization ************************** */

DEFINE_INIT(initNeutral,d)
{
cell_t c;
Thread *t;
real x[ND_ND];
real z;
/* loop over all cell threads in the domain */
thread_loop_c(t,d)
{
/* loop over all cells */
begin_c_loop_all(c,t)
{
C_CENTROID(x,c,t);
z = x[2] + z0;
if (z > ablHeight){
z = ablHeight;
}

if (z > maxZInit){
C_U(c,t) = 0.0; /*x velocity */
C_V(c,t) = (uStar/vonKarman)*log(z/z0); /* y velocity */
C_W(c,t) = 0.0; /* z velocity */
C_K(c,t) = pow(uStar,2.0)/sqrt(Cmu); /* k */
C_D(c,t) = pow(uStar,3.0)/(vonKarman*z); /* epsilon */
C_P(c,t) = 0.0; /*Pressure*/
}
else{
C_U(c,t) = 0.0;
C_V(c,t) = initVelocity;
C_W(c,t) = 0.0;
C_K(c,t) = pow(uStar,2.0)/sqrt(Cmu);
/* C_K(c,t) = initK; */
C_D(c,t) = initEpsilon;
C_P(c,t) = 0.0;
}
}
end_c_loop_all(c,t)
}
}


/* ************************* Wall Functions ************************** */

/* Designed around u/uStar = 1/K*log(z/z0) ref: Improved k-e model and wall function formulation for the RANS simulation of ABL flows, Parente et al
Removes the need for multiplying z0 by 9.73/Cs and can thus use roughness lengths directly from ABL modelling with first cell height = 2*z0*/


DEFINE_WALL_FUNCTIONS(ABL_logLaw, f, t, c0, t0, wf_ret, yPlus, Emod)
{
real ustar_ground, E_prime, yPlus_prime, zp, dx_mag, wf_value;
real mu=C_MU_L(c0,t0);
real xf[ND_ND];
real xc[ND_ND];
real dx[ND_ND];

F_CENTROID(xf, f, t);
C_CENTROID(xc, c0,t0);

dx[0] = xc[0] - xf[0];
dx[1] = xc[1] - xf[1];
dx[2] = xc[2] - xf[2];
dx_mag = NV_MAG(dx);
zp = dx_mag;

ustar_ground = pow(C_K(c0,t0),0.5)*pow(Cmu, 0.25);
E_prime = (mu/densOper)/(z0*ustar_ground);
yPlus_prime = (zp+z0)*ustar_ground/(mu/densOper);

switch (wf_ret)
{
case UPLUS_LAM:
wf_value = yPlus;
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