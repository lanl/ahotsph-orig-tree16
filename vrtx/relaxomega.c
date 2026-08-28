/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

/* relaxomega.c */


/*
May 1995: particle strength update in Winckelmans' relaxation  scheme
(see thesis and JCP paper).
*/


#include "physics_vrtx.h"
#include "vop.h"

void RelaxOmega(bodyptr btab, int n, float relaxw) {
    int i;
    bodyptr bp;
    float omega[3];

    for (i = 0; i < n; i++) {
        bp = btab + i;

        /* divergence free vorticity field that the new particle strengths
           will ``best" approximate after all iterations of the relaxation
           scheme.  */

        omega[0] = Gradvel(bp)[2][1] - Gradvel(bp)[1][2];
        omega[1] = Gradvel(bp)[0][2] - Gradvel(bp)[2][0];
        omega[2] = Gradvel(bp)[1][0] - Gradvel(bp)[0][1];

        /* incremental change of particle strength for this iteration */

        Strength(bp)[0] += relaxw * Vol(bp) * (omega[0] - Omegat(bp)[0]);
        Strength(bp)[1] += relaxw * Vol(bp) * (omega[1] - Omegat(bp)[1]);
        Strength(bp)[2] += relaxw * Vol(bp) * (omega[2] - Omegat(bp)[2]);
    }
}
