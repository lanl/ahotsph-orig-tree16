/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#define NO_MSGS
/* Msgs here really slow things down, even if they aren't activated */
#include "Msgs.h"
#include "fastflpt.h"
#include "physics_n.h"
#include "tensop.h"
#include "timers.h"
#include "vop.h"

Counter_t GravCnt;

static float Eps2;

void set_eps(float eps) { Eps2 = eps * eps; }

void set_body(void *o, void *p) { memcpy(o, p, TBODYSZ); }

void do_grav2(void *p0, void *list, int bsize, int n) {
    body *p = p0;
    int ncut;
    float total_mass;
    char *clist = (char *)list;
#ifdef __PARAGON__
    int nn;
#endif

    AddCounter(&GravCnt, n);

#ifdef __PARAGON__
    nn = n - (n % 3);
    /* Use the interface to the fast assembly code */
    if ((int)list & 07 || bsize & 07)
        Error("btab not aligned for asm code\n");
    do_grav_fast(
        list, (float *)(clist + nn * bsize), p->pos, &total_mass, p->acc, &p->phi, &Eps2, &ncut);
    if (n % 3)
        do_grav((char *)list + nn * bsize,
                (float *)(clist + n * bsize),
                p->pos,
                &total_mass,
                p->acc,
                &p->phi,
                &Eps2,
                &ncut);
#else
    do_grav(list, (float *)(clist + n * bsize), p->pos, &total_mass, p->acc, &p->phi, &Eps2, &ncut);
#endif
}
