/*purpose: take in any necessary values from calling routine, 
and calculate the cooling, and return the cooling value to 
the calling routine-note
-CE*/
#include <math.h>
#include <stdlib.h>
//#include "physics_sph.h"
//#include "vop.h"
#include "fastflpt.h"
//#include "timers.h"
//#include "error.h"
#include "cool.h"

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
#define NKERNEL_TABLE 80000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

//Counter_t SPHCnt, SPHrej, nbrMACCnt;

/*do I need any of these? -CE*/
//static float dvtable; /* == 0.0001 ... */
//static float invdvtable; /* == 10000.0 ... */
//static float cnormk;
//static float wij[MAX_INDEX];
//static float grwij[MAX_INDEX];
//static float Gamma = (float)(5.0/3.0);
//static float alpha = (float)1.0;
//static float beta = (float)2.5;
//static float epsil = (float)1e-2;
//static float heatf1 = (float)1.0;
//static int ndim;
//static int Nobj;
//static int add_offset;
//static float offset[NDIM];
//static float voffset[NDIM];
//static void (*bodyfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);
//static void (*cellfunc)(SinkSPH *sink, hcell **src_vec, int *res, int n);

/*
extern int do_diffusion;
extern int do_cooling;
*/


//void calc_lcool(double temp, double lcool)

double calc_lcool(double temp)
{
	double lcool;

	/*eventually replace this by a table look-up? -CE*/
	/* From Chris's email; fit to Dalgarno and McCray (ARA&A
	   1972, 10, 375) and Sutherland and Dopita (ApJS, 88,
	   253) */
	if (temp < 1.0e4)
		lcool = 1.0e-27 * exp(-1.0e2/temp) * sqrtf_fast(temp);
	else if (temp < 3.0e5)
		lcool=1.0e-21;
	else
		lcool=1.0e-21/(3.0*(log10(temp)-5.5)+1.0);

	//return;
	return lcool;
}
