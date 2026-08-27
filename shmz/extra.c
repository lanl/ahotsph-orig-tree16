/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include "extra.h"

#include "bigmalloc.h"
#include "chn.h"
#include "error.h"
#include "physics.h"

#define N_PER_BLOCK 16
#define NMAX 16
static Chn ids[NMAX];
static int sizes[NMAX];

static void *chn_extra_alloc(int n) {
    int i;
    if (n <= 0)
        Error("Bad size\n");
    if (n % 8)
        Error("Alignment problems possible\n");

    for (i = 0; i < NMAX; i++) {
        if (sizes[i] == n)
            return ChnAlloc(ids + i);
        else if (sizes[i] == 0) {
            sizes[i] = n;
            ChnInit(ids + i, n, N_PER_BLOCK, Realloc_f);
            return ChnAlloc(ids + i);
        }
    }
    Error("Too many different sizes, increase NMAX\n");
}

void CellExtraFree(void) {
    int i;

    for (i = 0; i < NMAX; i++) {
        if (sizes[i] != 0)
            ChnTerminate(ids + i);
    }
}

int CellExtraSz(void *pp) {
    int sz;
    cell *p = pp;
    sz = p->nu * p->nv * sizeof(complex);
    if (sz < 0 || sz > 1024 * 1024)
        Error("Bad sz (%d)\n", sz);
    return (sz);
}

void *CellExtraPtr(void *pp) {
    cell *p = pp;
    return (p->ffsf);
}

void *CellExtraAlloc(void *pp) {
    cell *p = pp;
    void *ptr;
    if (p->nu > 256 || p->nu < 0)
        Error("Bad nu (%d)\n", p->nu);
    if (p->nv > 256 || p->nv < 0)
        Error("Bad nv (%d)\n", p->nv);
    ptr = chn_extra_alloc(p->nu * p->nv * sizeof(complex));
    p->ffsf = ptr;
    return (ptr);
}
