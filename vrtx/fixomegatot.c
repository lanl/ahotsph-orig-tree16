/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

/* fixomegatot.c */


/*
October 1995: Fix the particle strengths so as to make sure that the sum of all
particle strengths is zero.
*/

#include "physics_vrtx.h"
#include "vop.h"

extern double omega_tot[3];

void FixOmegaTot(bodyptr btab, int nobj, int gnobj) {
    int i;
    bodyptr bp;
    double term[3];


    term[0] = omega_tot[0] / gnobj;
    term[1] = omega_tot[1] / gnobj;
    term[2] = omega_tot[2] / gnobj;

    for (i = 0; i < nobj; i++) {
        bp = btab + i;

        VV(Strength(bp), -= term);
    }
}
