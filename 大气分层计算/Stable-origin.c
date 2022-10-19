#include "udf.h"
#include "math.h"

/* **************************************************************
** Stable **
**************************************************************
Fluent UDFs for simulating stable ABL flow

Control via the defined parameters
Ensure the solver is in expert mode
Use compiled UDF method

C_UDMI - 12 User memory slots, 1 User scalar slot
Wall Distance
Cor x
Cor y
k DTU
k Dtu Norm
epsilon Fluent
epsilon AM
epsilon AM Ce3
epsilon AM Gb
epsilon DTU
epsilon DTU - Ce3
DTU Gb

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
#define uStar 0.1407
#define Lin 124.7334 /* L - Inlet L must be > 0 to use this UDF set!!! */
#define Lmast 222.1774 /* L at the mast position Interpolation is performed from the inlety to the mast for the L values so that at the inlet the value is Lin and at the mast the value is L mast */
#define T0 288.0
#define Tstar 0.0232
#define ablHeight 600.0 /* Height of ABL, this is the height for fixed values of all profiles */

/* Site */
#define globalLat -33.0 /* Latitude of the origin in degrees - This is a dummy value for confidentiality*/
#define siteElevation 0.0 /* Altitude of site AMSL - If you specify the operating pressure from site data then DO NOT change this value. */
#define earthRot 0.000072921159 /* Earth rotational speed */
#define offset 477.0 /* Use to control the z value, this is deducted from the mesh z coordinate. This is the height AGL of the inlet location of the mesh */
#define offsetY -3000.0 /* This is deducted from the local lattitude in the corliolis calculation */
#define mastLocation 8687.0


/* General */
#define pi 3.141592
#define g -9.80665
#define R 8.3144598 /* Universal Gas Constant - Dry Air */
#define M 0.0289644 /* Molar mass of Earth's air */
#define Lb -0.0065 /* Standard temperature lapse rate */

/* Operating Conditions - Material Air */
#define presOper 101325 /* Operating Pressure Pa - Internal Solver Pressure. This is the pressure specified at 0m and for this you can use lowest mast pressure reading */
#define tempOper 288.16 /* Operating Temperature - Internal Solver Standard Tempearture. This is the temperature based from the lowest measurement height on the mast. But can be left as the standard value */
#define densOper 1.0800 /* Problem density */
#define Cp 1006.43
#define beta 0.0032
#define viscosity 1.7894e-05

/* Initilization */
/* Due to HAGL variations and Fluent not being able to compute cell distance before initialiazing we have to manually set the initialiaze values. These are used for z values lower than maxZInit, afterwards it returns to the inlt profile values */
#define maxZInit 1000.0 /* Height before using init values from inlet profiles */
#define initVelocity 10.0 /* y velocity */
#define initK 2.0 /* k */
#define initEpsilon 2.0 /* epsilon */

double linearInterpolation(double y);
/* ********************* Profiles ****************************** */

/* ********************** Inlet Velocity ********************** */
DEFINE_PROFILE(inletVelocityStable,t,i)
{
real x[ND_ND];
real z;
real phiM;
face_t f;

begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
z = x[2] + z0 - offset;
if (z > ablHeight){
z = ablHeight;
}
phiM = 1.0 + 5.0*(z/Lin);
F_PROFILE(f, t, i) = (uStar/vonKarman)*(log(z/z0) +phiM -1.0);
}
end_f_loop(f,t)
}


/* ********************** Inlet Temperature **********************
The site values for temperature are in potential temperature. This is converted back to standard temperature via the operating pressure of Fluent. */
DEFINE_PROFILE(inletTemperatureStable,t,i)
{
real x[ND_ND];
real z;
real phiM, potenTemp, pressure, zAMSL;
face_t f;

begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
z = x[2] + z0 - offset;
if (z > ablHeight){
z = ablHeight;
}
zAMSL = z + siteElevation;
phiM = 1.0 + 5.0*(z/Lin);
potenTemp = T0 + (Tstar/vonKarman)*(log(z/z0) +phiM -1.0);
pressure = presOper*pow(tempOper/(tempOper+Lb*zAMSL),(-g*M)/(R*Lb));
F_PROFILE(f,t,i) = potenTemp/(pow(presOper/pressure,0.286));
}
end_f_loop(f,t)
}


/* ********************** Inlet k ********************** */
DEFINE_PROFILE(inlet_k_Stable,t,i)
{
real x[ND_ND];
real z;
real phiE,phiM;
face_t f;

begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
z = x[2] + z0 - offset;
if (z > ablHeight){
z = ablHeight;
}
phiM = 1.0 + 5.0*(z/Lin);
phiE = phiM-z/Lin;
F_PROFILE(f,t,i) = (pow(uStar,2.0)/sqrt(Cmu))*pow(phiE/phiM,0.5);
}
end_f_loop(f,t)
}

/* ********************** Inlet epsilon ********************** */
DEFINE_PROFILE(inlet_e_Stable,t,i)
{
real x[ND_ND];
real z;
real phiE, phiM;
face_t f;

begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
z = x[2] + z0 - offset;
if (z > ablHeight){
z = ablHeight;
}
phiM = 1.0 + 5.0*(z/Lin);
phiE = phiM-z/Lin;
F_PROFILE(f,t,i) = phiE*pow(uStar,3.0)/(vonKarman*z);
}
end_f_loop(f,t)
}


/* *********************** Walls *********************** */

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

/* ********************** Wall Temperature ********************** */
DEFINE_PROFILE(wallTemperatureStable,t,i)
{
real x[ND_ND];
face_t f;
begin_f_loop(f,t)
{
F_CENTROID(x,f,t);
F_PROFILE(f,t,i) = T0;
}
end_f_loop(f,t)
}

/* ************************* Cell Wall Distance ************************** */
/* To Use: Define a UDS with Flux Function = none, no Inlet Diffusion
Add Material Property "UDS Diffusivity"; defined-per-uds: constant, Coefficient = 1 [kg/ms]
Add Source Terms for User Scalars in the cell zone: Source Term = 1
Set Boundary Conditions for User Scalar: Specified Value = 0 on all boundaries to which the distance should be computed (boundary lower in the attached sample case); Specified Flux = 0 on all other boundaries.
Define a User-Defined Memory Location in which the UDF stores the computed distance
Hook to define_excecute_at_end */

DEFINE_EXECUTE_AT_END(computeSelectedWallDistance)
{
Domain *d=Get_Domain(1);
Thread *t;
cell_t c;
real wallPhi, gradWallPhi, wallDistance;

/* Check if UDM and UDS exist */
if (N_UDM < 12 || N_UDS < 1) {
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

/* ******************** Sources ********************* */

/* **************** Corliolis Force *************** */
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

/* ************************ k **************************** */
/* No energy eqaution is solved with this model*/
DEFINE_SOURCE(k_source_DTU_Stable,c,t,dS,eqn)
{
real fSt, phiM, phiE, phiH, CkD, source, Gb, Sk, uStarLocal;
real x[ND_ND];
real z, L;
C_CENTROID(x,c,t);
z = C_UDMI(c,t,0) + z0;
L = linearInterpolation(x[1]);
if (z > ablHeight){
z = ablHeight;
}

if (N_ITER > 5) {
phiM = 1.0 + 5.0*(z/L);
phiE = phiM-z/L;
phiH = 1.0 + 5.0*(z/L);
uStarLocal = pow(C_K(c,t),0.5)*pow(Cmu,0.25)*pow(phiM,0.25)*pow(phiE,-0.25);

fSt = 2.0-(z/L) - 10.0*(z/L)*(1.0-2.0*(z/L) + 10.0*(z/L));
CkD = pow(vonKarman,2.0)/(sigma_k*sqrt(Cmu));
Gb = -C_MU_T(c,t)*pow(sqrt(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2]),2.0)*((z/L)/(sigma_theta))*(phiH/pow(phiM,2.0)); /* DTU Formulation */
Sk = pow(uStarLocal,3.0)/(vonKarman*L)*(1.0 - (phiH)/(sigma_theta*phiM) - 0.25*CkD*pow(phiM,-3.5)*pow(phiE,-1.5)*fSt);
source = -densOper*Sk + Gb;
}
else {
source = 0.0;
}

dS[eqn] = 0.0;
C_UDMI(c,t,3) = Sk;
C_UDMI(c,t,4) = Sk*vonKarman*z/pow(uStar,3.0);
return source;
}


/* ************************ Epsilon **********************
Epsilon is a function of the gradients and to save these the solver needs to be in expert mode
Issue: 'solve/set/expert' in the FLUENT window, and answer YES when it asks if you want to free temporary memory

Standard Fluent buoyancy treatment for epsilon
Checking advanced buoyancy treatmnent in the viscous model box adds in the formulation below
Changes in the model is made by changing Ce3 according to the AM or DTU method
Not checking the box sets Gb = 0, this term is then re added in by the sources below. Do not check the box in the viscous box! */
DEFINE_SOURCE(epsilon_source_Fluent_Stable,c,t,dS,eqn)
{
real Gb, C3e, source;

if (N_ITER > 5) {
Gb = beta*g*C_MU_T(c,t)/PrTurb*C_T_G(c,t)[2]; /* Standard Fluent Gb formulation, C_MU_T = Turbulent Viscosity, PrTurb = Turbulent Prandtl number, C_T_G = [partial_T/partial_xi] */
C3e = tanh(fabs(C_V(c,t)/C_U(c,t))); /* Standard Fluent C3e formulation, C_V = v velocity, C_U = x velocity */
source = Ce1*C_D(c,t)/C_K(c,t)*C3e*Gb; /* C_D = epsilon, C_K = k */
}
else {
source = 0.0; /* Only run this source after 15 iterations. The gradients can cuase divergence with an illposed initilization */
}
dS[eqn] = 0.0;
C_UDMI(c,t,5) = source;
return source;
}


/* ALot & Masson
/* Epsilon source treatment based on an anylytical expression for Ce3 */
/* Only valid of -2.3 < z/L < 2 and also highly sensitive*/
DEFINE_SOURCE(epsilon_source_AM_Stable,c,t,dS,eqn)
{
real x[ND_ND];
real z, L;
real Gb, C3e, source;
real a0, a1, a2, a3, a4, a5;
C_CENTROID(x,c,t);
z = C_UDMI(c,t,0) + z0;
L = linearInterpolation(x[1]);
if (z > ablHeight){
z = ablHeight;
}

if (N_ITER > 5 && z/L < 2.0) {
if (z/L < 0.33) {
a0 = 4.181;
a1 = 33.994;
a2 = -442.398;
a3 = 2368.12;
a4 = -6043.544;
a5 = 5970.776;
}
else {
a0 = 5.225;
a1 = -5.269;
a2 = 5.115;
a3 = -2.406;
a4 = 0.435;
a5 = 0;
}

Gb = beta*g*C_MU_T(c,t)/PrTurb*C_T_G(c,t)[2];
C3e = a0*pow((z/L),0) + a1*pow((z/L),1.0) + a2*pow((z/L),2.0) + a3*pow((z/L),3.0) + a4*pow((z/L),4.0) + a5*pow((z/L),5.0); /* AM C3e formulation */
}
else if (N_ITER > 5 && z/L > 2.0){
Gb = beta*g*C_MU_T(c,t)/PrTurb*C_T_G(c,t)[2];
C3e = 2.858999999999999;
}
else{
Gb = 0.0;
C3e = 0.0;
}

dS[eqn] = 0;
source = Ce1AM*C_D(c,t)/C_K(c,t)*C3e*Gb;
C_UDMI(c,t,6) = source;
C_UDMI(c,t,7) = C3e;
C_UDMI(c,t,8) = Gb*vonKarmanAM*z/pow(uStar,3.0);
return source;
}

/* 2 - This uses the DTU Gb formulation and is run without a temperature eqaution */
DEFINE_SOURCE(epsilon_source_AM_Stable_2,c,t,dS,eqn)
{
real x[ND_ND];
real z, L;
real Gb, C3e, phiM, phiH, source;
real a0, a1, a2, a3, a4, a5;
C_CENTROID(x,c,t);
z = C_UDMI(c,t,0) + z0;
L = linearInterpolation(x[1]);

if (z > ablHeight){
z = ablHeight;
}
phiM = 1.0 + 5.0*(z/L);
phiH = 1.0 + 5.0*(z/L);

if (N_ITER > 5 && z/L < 2.0) {
if (z/L < 0.33) {
a0 = 4.181;
a1 = 33.994;
a2 = -442.398;
a3 = 2368.12;
a4 = -6043.544;
a5 = 5970.776;
}
else {
a0 = 5.225;
a1 = -5.269;
a2 = 5.115;
a3 = -2.406;
a4 = 0.435;
a5 = 0;
}

Gb = C_MU_T(c,t)*pow(sqrt(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2]),2.0)*((z/L)/(sigma_theta))*(phiH/pow(phiM,2.0));
C3e = a0*pow((z/L),0) + a1*pow((z/L),1.0) + a2*pow((z/L),2.0) + a3*pow((z/L),3.0) + a4*pow((z/L),4.0) + a5*pow((z/L),5.0); /* AM C3e formulation */
}
else if (N_ITER > 5 && z/L >= 2.0){
Gb = C_MU_T(c,t)*pow(sqrt(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2]),2.0)*((z/L)/(sigma_theta))*(phiH/pow(phiM,2.0));
C3e = 2.858999999999999;
}
else{
Gb = 0.0;
C3e = 0.0;
}

dS[eqn] = 0;
source = Ce1AM*C_D(c,t)/C_K(c,t)*C3e*Gb;
C_UDMI(c,t,6) = source;
C_UDMI(c,t,7) = C3e;
C_UDMI(c,t,8) = Gb*vonKarmanAM*z/pow(uStar,3.0);
return source;
}

/* DTU
Epsilon source treatment based on an anylytical expression for Ce3 */
DEFINE_SOURCE(epsilon_source_DTU_Stable,c,t,dS,eqn)
{
real x[ND_ND];
real z, L;
real Gb, C3e, source;
real phiM, phiH, phiE, fe;
C_CENTROID(x,c,t);
z = C_UDMI(c,t,0) + z0;
L = linearInterpolation(x[1]);

if (z > ablHeight){
z = ablHeight;
}

if (N_ITER > 5) {
phiM = 1.0 + 5.0*(z/L);
phiE = phiM-z/L;
phiH = 1.0 + 5.0*(z/L);
fe = pow(phiM,-2.5)*(2.0*phiM-1.0);

Gb = -C_MU_T(c,t)*pow(sqrt(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2]),2.0)*((z/L)/(sigma_theta))*(phiH/pow(phiM,2.0)); /* DTU Formulation */
C3e = sigma_theta/(z/L)*(phiM/phiH)*(Ce1*phiM-Ce2*phiE+(Ce2-Ce1)*pow(phiE,-0.5)*fe);
source = Ce1*C_D(c,t)/C_K(c,t)*C3e*Gb;
}
else {
C3e = 0.0;
Gb = 0.0;
source = 0.0;
}

C_UDMI(c,t,9) = source;
C_UDMI(c,t,10) = C3e;
C_UDMI(c,t,11) = Gb*vonKarman*z/pow(uStar,3.0);
dS[eqn] = 0.0;
return source;
}

/* ************************* Initilization ************************** */

DEFINE_INIT(initStable,d)
{
cell_t c;
Thread *t;
real x[ND_ND];
real phiM, phiE, phiH, pressure, potenTemp, z, zAMSL, L ;
/* loop over all cell threads in the domain */
thread_loop_c(t,d)
{
/* loop over all cells */
begin_c_loop_all(c,t)
{
C_CENTROID(x,c,t);
z = x[2] + z0;
L = linearInterpolation(x[1]);
if (z > ablHeight){
z = ablHeight;
}

if (z > maxZInit){
phiM = 1.0 + 5.0*(z/L);
phiE = phiM-z/L;
phiH = 1.0 + 5.0*(z/L);
C_U(c,t) = 0.0; /*x velocity */
C_V(c,t) = (uStar/vonKarman)*(log(z/z0) +phiM -1.0); /* y velocity */
C_W(c,t) = 0.0; /* z velocity */
/* C_T(c,t) = potenTemp/(pow(presOper/pressure,0.286)); */ /* Temperature */
C_K(c,t) = (pow(uStar,2.0)/sqrt(Cmu))*pow(phiE/phiM,0.5); /* k */
C_D(c,t) = phiE*pow(uStar,3.0)/(vonKarman*z); /* epsilon */
C_P(c,t) = 0.0; /*Pressure*/
}
else{
C_U(c,t) = 0.0;
C_V(c,t) = initVelocity;
C_W(c,t) = 0.0;
C_K(c,t) = initK;
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
wf_value =
break;
case D2UPLUS_TRB:
wf_value = -1.0/(vonKarman*yPlus_prime*yPlus_prime);
break;
default:
printf("Wall function return value unavailable\n");
}
return wf_value;
}


/* ************************* Interpolation ************************** */
/* Currently does linear interpolation, Must be run with 180degree inlet location. This function can be expanded in future to bilinear (or more) to include more mast/WRF locations */
double linearInterpolation(double y)
{
double L;
if (y > mastLocation){
L = Lmast;
}
else{
L = (Lin*(mastLocation - y) + Lmast*(y - offsetY))/(mastLocation - offsetY); /* Local L */
}

return L;
}