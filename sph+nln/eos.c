#include <stdio.h>
#include <math.h>
#include <Msgs.h>
#include "singlio.h"
#include "units.h"

#define MAX_ITER 200

/* Need to change printf's to Msgfs or whatever */

float newtraph(double xl, double xr, double prec, double (*f)(double x), 
	     double (*df)(double x)) {
    double xguess = 2.0*xl*xr/(xl+xr);
    int i;
   
    if ((*f)(xl)*(*f)(xr) > 0.0) {
	singlPrintf("Bisect: %1.1e and %1.1e do not bracket a root\n", xl, xr);
        xguess = -99.0;
        return (float)xguess;
    }

    for(i = 0; i <= MAX_ITER; i++) {
	xguess = xguess - (*f)(xguess)/(*df)(xguess); /* Try Newton-Raphson */

	if (xguess < xl || xguess > xr)
	    xguess = (xl + xr) / 2.0; /* Fall back to bisection */
    
	if ( fabs((*f)(xguess)) < prec )
	    break;
	else if ((*f)(xl)*(*f)(xguess) < 0.0)
	    xr = xguess;
	else
	    xl = xguess;
    }

    if (i == MAX_ITER) {
	singlPrintf("Bisect: max iterations exceeded\n");
        xguess = -99.0;
    }

    return (float)xguess;
}


/* changed: these now need to be in cgs */
extern double eos_n;
extern double eos_u;

#include "physics_sph.h"

double uvst(double t) {
    return 1.5*eos_n*(K_BOLTZ)*t + t*t*t*t*(A_RAD) - eos_u;
}


double duvst(double t) {
    return 1.5*eos_n*(K_BOLTZ) + 4.0*t*t*t*(A_RAD);
}

/* pressure from different eos's */
/* eta = rho/rho_0 */
/* Schaefer et al. 2016, A&A */
double liquid_eos (double k_bulk, double eta) {
	return k_bulk * (eta - 1.0);
}

/* Schaefer et al. 2016, A&A */
double murnaghan_eos(double k_bulk, double n_M, double eta) {
	return k_bulk / n_M * (pow (eta, n_M) - 1.0);
}

/* Schaefer et al. 2016, A&A */
double tillotson_eos(double A_T, double B_T, double E_0, double a_T, double b_T, double alpha_T, double beta_T, double eta, double u, double rho) {
	/* first ~7 args are 'material constants' */
	double pressure;

	if (eta > 1.0) { /* compression */
		pressure = (a_T + b_T / (1 + u / (E_0 * eta * eta))) * rho * u;
		pressure += A_T * (eta - 1.0) + B_T * (eta - 1.0) * (eta - 1.0);
	} else { /* expansion/vaporization? */
		pressure = a_T * rho * u;
		pressure += (b_T * rho * u / (1.0 + u / (E_0 * eta * eta)) + 
					A_T * (eta - 1.0) * exp (-beta_T * (1./eta - 1.0)));
		pressure *= exp (-alpha_T * (1./eta - 1.0) * (1./eta - 1.0));
	}

	return pressure;
}

/* Wikipedia on Anton-Schmidt equation of state */
double anton_schmidt_eos(double k_bulk, double power_n, double eta) {
	/* here, technically eta = Vol/Vol_0, but mass should be 
	 * constant, so close enough */
	return -k_bulk * pow (eta, power_n) * log (eta);
}

