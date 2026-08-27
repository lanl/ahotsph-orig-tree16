/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a) {
    void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
    int j;

    for (j = 1; j <= m; j++) spline(x2a, ya[j], n, 1.0e30, 1.0e30, y2a[j]);
}
