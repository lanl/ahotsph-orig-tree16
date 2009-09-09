#include <stdio.h>
#include <math.h>
#include <Msgs.h>
#include "units.h"

#define MAX_ITER 100

/* Need to change printf's to Msgfs or whatever */

float newtraph(double xl, double xr, double prec, double (*f)(double x), 
	     double (*df)(double x)) {
    double xguess = 2.0*xl*xr/(xl+xr);
    int i;
   
    Msgf(("in eos, K_BOLTZ=%g\n",K_BOLTZ));

    if ((*f)(xl)*(*f)(xr) > 0.0) {
	Msgf(("Bisect: %1.1e and %1.1e do not bracket a root\n", xl, xr));
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

    if (i == MAX_ITER) 
	Msgf(("Bisect: max iterations exceeded\n"));

    return (float)xguess;
}


#include "physics_sph.h"

extern double eos_n;
extern double eos_u;

double uvst(double t) {
    /*convert from cgs to user-units*/
    double kB=K_BOLTZ /((double)massCF /lengthCF/lengthCF *timeCF*timeCF);
    double acoef=A_COEFF/((double)massCF*lengthCF*timeCF*timeCF);
    return 1.5*eos_n*(kB)*t + t*t*t*t*(acoef) - eos_u;
}


double duvst(double t) {
    /*convert from cgs to user-units*/
    double kB=K_BOLTZ /((double)massCF /lengthCF/lengthCF *timeCF*timeCF);
    double acoef=A_COEFF/((double)massCF*lengthCF*timeCF*timeCF);
    return 1.5*eos_n*(kB) + 4.0*t*t*t*(acoef);
}
