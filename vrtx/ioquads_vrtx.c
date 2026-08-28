/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

/* Ensure that physics_vrtx.h won't be included again by quads_vrtx.h */
#include "physics_vrtx.h"

#define body iobody
#define bodyptr iobodyptr

typedef struct {
    float pos[3];
    int nterms;
    float psi[3];
    float vel[3];
    float gradvel[3][3];
    float dstr[3];
    float errsum, errsum2;
} body, *bodyptr;

/* And tweak all the external symbols in quads_vrtx.c */
#define InteractCell InteractIOCell
#define CellCnt IOCellCnt

#include "quads_vrtx.c"
