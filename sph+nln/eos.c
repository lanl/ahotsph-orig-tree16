/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <Msgs.h>
#include <math.h>
#include <stdio.h>

#include "singlio.h"
#include "units.h"

#define MAX_ITER 200

/* Need to change printf's to Msgfs or whatever */

float newtraph(double xl, double xr, double prec, double (*f)(double x), double (*df)(double x)) {
    double xguess = 2.0 * xl * xr / (xl + xr);
    int i;

    if ((*f)(xl) * (*f)(xr) > 0.0) {
        singlPrintf("Bisect: %1.1e and %1.1e do not bracket a root\n", xl, xr);
        xguess = -99.0;
        return (float)xguess;
    }

    for (i = 0; i <= MAX_ITER; i++) {
        xguess = xguess - (*f)(xguess) / (*df)(xguess); /* Try Newton-Raphson */

        if (xguess < xl || xguess > xr)
            xguess = (xl + xr) / 2.0; /* Fall back to bisection */

        if (fabs((*f)(xguess)) < prec)
            break;
        else if ((*f)(xl) * (*f)(xguess) < 0.0)
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

double uvst(double t) { return 1.5 * eos_n * (K_BOLTZ)*t + t * t * t * t * (A_RAD)-eos_u; }


double duvst(double t) { return 1.5 * eos_n * (K_BOLTZ) + 4.0 * t * t * t * (A_RAD); }


/* pressure from different eos's */
/* eta = rho/rho_0 */
/* Schaefer et al. 2016, A&A */
double liquid_eos(double k_bulk, double eta) { return k_bulk * (eta - 1.0); }

/* Schaefer et al. 2016, A&A */
double murnaghan_eos(double k_bulk, double n_M, double eta) {
    return k_bulk / n_M * (pow(eta, n_M) - 1.0);
}

/* possibly basalt? */
// void setconst1(Material_t *m) {
//	m->rho0 = 2.86;
//    m->A = 1.80e11;
//	m->B = 1.80e11;
//    m->a = 0.5;
//    m->b = 1.3;
//    m->alpha = 5.;
//    m->beta = 5.;
//    m->u0 = 1.60e11;
//    m->Eiv = 3.50e10; /* E. of incipient vaporization */
//    m->Ecv = 1.80e11; /* E. of complete vaporizaton */
//    m->mu = 2.5e11;
//    m->umelt = 3.0e10;
//    m->yield = 3.5e10;
//    m->pweib = 6.2;
//    m->cweib = 1.e27;
//}

/* possibly iron? */
//  void setconst2(Material_t *m)
//{
//    //m->rho0 = 7.86;
//    m->rho0 = 8.05;
//    m->A = 1.28e12;
//    m->B = 1.05e12;
//    m->a = 0.5; /* like (gamma - 1.) term ?? */
//    m->b = 1.5;
//    m->alpha = 5.;
//    m->beta = 5.;
//    //m->u0 = 9.50e10;
//    m->u0 = 1.50e09;
//    m->Eiv = 1.42e10;
//    m->Ecv = 8.45e10;
//    m->mu = 0.0;
//    m->umelt = 1.0e10;
//    m->yield = 6.0e9;
//    m->pweib = 9.0;
//    m->cweib = 0.0;
//}

void tillotson_eos(double rho, double u, Material_t *m, double *pressure, double *cs) {
    double PC = 0.;  /* pressure, current */
    double csC = 0.; /* sound speed, current */
    double rho0m1 = 1. / m->rho0;
    double rhom1 = 1. / rho;
    double eta = rho * rho0m1;
    double mu = eta - 1.;
    double csmin = 0.25 * m->A * rho0m1;
    double Pmin = 0.; /* set floor on pressure */
    double c1 = u / (m->u0 * eta * eta);
    double c2 = 1. / (c1 + 1.);

    /* A. Brundage, 2013, 12th Hypervelocity Impact Symposium,
     * Procedia Engineering, 58, 461-470 describes an implementation of this eos - CIE */
    if (u > m->Eiv && eta < 1.) {
        /* eq. 3, but there u> u_cv */
        double d1 = m->rho0 * rhom1; /* = 1/eta ? */
        double d2 = d1 - 1.;
        double ex1 = (m->beta * d2 < 60) ? exp(-m->beta * d2) : 0.0;
        double ex2 = (m->alpha * d2 * d2 < 60) ? exp(-m->alpha * d2 * d2) : 0.0;
        *pressure = m->a * rho * u;
        *pressure += ex2 * (m->b * rho * u * c2 + m->A * mu * ex1);

        /* sound speed? */
        *cs = m->b * u * (3. * c1 + 1.) * c2 * c2 + 2. * m->alpha * d2 * m->b * d1 * u * c2;
        *cs += m->A * ex1 * ((2. * m->alpha * d2 + m->beta) * mu * d1 * rhom1 + rho0m1);
        *cs = *cs * ex2 + m->a * u;
        *cs += *pressure * rhom1 * (m->a + m->b * c2 * c2 * ex2);
        if (*cs < 0.)
            *cs = 0.;
    }

    if (u < m->Ecv || eta >= 1.) {
        /* eq. 2 */
        PC = (m->a + m->b * c2) * rho * u + m->A * mu + m->B * mu * mu;
        csC = m->a * u + rho0m1 * (m->A + 2. * m->B * mu) + m->b * u * (3. * c1 + 1.) * c2 * c2;
        csC += PC * rhom1 * (m->a + m->b * c2 * c2);

        if (u > m->Eiv && u < m->Ecv && eta < 1) {
            double e1 = m->Ecv - u;
            double e2 = u - m->Eiv;
            double e3 = 1. / (m->Ecv - m->Eiv);
            /* eq. 5 */
            *pressure = (e2 * *pressure + e1 * PC) * e3;
            *cs = (e2 * *cs + e1 * csC) * e3;
        } else {
            *pressure = PC;
            *cs = csC;
        }
    }
    if (*cs < csmin)
        *cs = csmin;
    if (*pressure < Pmin)
        *pressure = Pmin;
    *cs = sqrt(*cs);
}

/* Wikipedia on Anton-Schmidt equation of state */
double anton_schmidt_eos(double k_bulk, double power_n, double eta) {
    /* here, technically eta = Vol/Vol_0, but mass should be
     * constant, so close enough */
    return -k_bulk * pow(eta, power_n) * log(eta);
}
