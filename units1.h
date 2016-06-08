/* all physical constants go in here */
/* all physical constants IN CGS UNITS PLEASE!!! */
/* the conversion factors assume that all internal physics is done 
 * in cgs, and convert FROM CGS TO USER-UNITS. They are set in the 
 * .ctl file. ADD */

extern float massCF;
extern float lengthCF;
extern float timeCF;

#ifndef CONSTS
#define CONSTS
double GRAV_C =6.67428e-8;
double A_COEFF =5.760400e-05;
double C_LIGHT =2.99792458e+10;
double KES_COEFF =6.6524586e-25;
double KFF_COEFF =0.640000e+23;
double K_BOLTZ =1.3806503e-18;
double MH =1.67262158e-24;
#else
extern double GRAV_C;
extern double A_COEFF;
extern double C_LIGHT;
extern double KES_COEFF;
extern double KFF_COEFF;
extern double K_BOLTZ;
extern double MH;
#endif
