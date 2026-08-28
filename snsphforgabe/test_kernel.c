/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <stdio.h>

/* External Fortran linkage */
#define Fortran(x) x##_
/* GNU Fortran adds two underscores if there is an underscore in the name */
#ifdef __GNUC__
#define Fortran2(x) x##__
#else
#define Fortran2(x) x##_
#endif


#define ITABLE 80000
extern void Fortran(ktable)(void);
/* common /table/ wij(0:itable), grwij(0:itable), dvtable */
extern float Fortran(table)[ITABLE * 2];

main(int argc, char *argv) {
    int i;
    float *wij, *grwij;

    Fortran(ktable)();
    wij = Fortran(table);
    grwij = Fortran(table) + ITABLE + 1;
    for (i = 0; i < ITABLE; i++) { printf("%5d %8g %8g\n", i, wij[i], grwij[i]); }
}
