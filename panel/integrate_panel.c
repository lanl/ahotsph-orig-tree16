/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

/* Panel code */


#include <math.h>

#include "mpmy.h"
#include "physics_panel.h"
#include "vop.h"


void Update(bodyptr btab, int n, float relax, float *residual) {
    int i;
    bodyptr bp;
    double sumsq, sumarea;
    float Veltot[3], uzl;
    MPMY_Comm_request req;

    sumsq = 0.;
    sumarea = 0.;
    for (i = 0; i < n; i++) {
        bp = btab + i;


        /* GSW Apr16 93: So far, it uses Over-Relaxed Jacobi */

        /* velocity induced by all panels with their present sources (this was
           obtained from the ""Tree Code"!) + freestream velocity.
           ... in absolute coordinate system:  */

        VVV(Veltot, = Vel(bp), +Uext(bp));

        /* ... in local coordinate system: */

        uzl = Dot(Ez(bp), Veltot);

        /* Change in source

           Also, Recall that dzlsphi=2*pi...
        */
        /* GSW Apr16 93: So far, it uses Over-Relaxed Jacobi ! */

        /* sign?  There's only two possibilities... */
        bp->dsigma = relax * uzl / (6.283185308);
        bp->sigma += bp->dsigma;
        sumsq += uzl * uzl * Ip(bp);
        sumarea += Ip(bp);
        /* This is as good a place as any to re-zero the 'nterms' counter */
        /* and errsum accumulators in each body */
        bp->nterms = 1;
        bp->errsum = 0.;
        bp->errsum2 = 0.;
    }
    MPMY_ICombine_Init(&req);
    MPMY_ICombine(&sumsq, &sumsq, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine(&sumarea, &sumarea, 1, MPMY_DOUBLE, MPMY_SUM, req);
    MPMY_ICombine_Wait(req);
    *residual = sqrt(sumsq / sumarea);
}
