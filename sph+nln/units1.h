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
#define GRAV_C (6.67428e-8)
#define A_COEFF (5.760400e-05)
#define C_LIGHT (2.99792458e+10)
#define KES_COEFF (6.6524586e-25)
#define KFF_COEFF (0.640000e+23)
#define K_BOLTZ (1.3806503e-18)
#define MH (1.67262158e-24)
#endif
