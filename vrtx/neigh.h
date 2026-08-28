/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#ifndef NEIGH_dot_H
#define NEIGH_dot_H

/* structures and prototypes for the neighbor-finding code */

#include "physics_vrtx.h"
#include "timers.h"

#ifdef NDIM
#if NDIM != 3
#error this code assumes NDIM is 3
#endif
#else
#define NDIM 3
#endif

typedef struct {
    float pos[NDIM];
    float r2;
} ncell;

typedef struct {
    float center[NDIM];
    float size;
} ncofm;

typedef struct {
    body *bp;
} nsink;


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void NeighCofmFromDaugh(hcell *hptr, hcell **dlist);
void NeighCellFromCofm(ncell *cellp, ncofm *cofmp);
void NeighInherit(const nsink *from, nsink *to, hcell *pp);
void NeighMACv(nsink *sink, const hcell **srcs, int *results, int nsrc);

extern Counter_t NfindTestsCnt;
extern Counter_t NfindAcceptsCnt;

#ifdef __cplusplus
}
#endif

#endif
