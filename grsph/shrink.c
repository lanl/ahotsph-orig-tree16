/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include "physics_sph.h"
#include "stk.h"
#include "vop.h"

void ShrinkBtab(body **btabp, int *nobj, float r_limit) {
    body *btab = *btabp;
    body *p;
    Stk s;
    body *q;
    float r2;

    StkInitEz(&s);
    r2 = r_limit * r_limit;

    for (p = btab; p < btab + *nobj; p++) {
        if (Dot(p->pos, p->pos) >= r2) { /* acceptable */
            q = StkPush(&s, sizeof(body));
            *q = *p;
        }
    }
    Free(btab);
    StkCrunch(&s);
    *nobj = StkSz(&s) / sizeof(body);
    btab = StkBase(&s);
    *btabp = Realloc(btab, *nobj * sizeof(body));
}
