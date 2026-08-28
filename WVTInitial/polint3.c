/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <math.h>


void polint(double xa[], double ya[], int n, double x, double *y, double *dy) {
    int i, m, ns = 1;
    double den, dif, dift, ho, hp, w;
    double *c, *d;

    dif = fabs(x - xa[1]);
    c = vector(1, n);
    d = vector(1, n);
    for (i = 1; i <= n; i++) {
        if ((dift = fabs(x - xa[i])) < dif) {
            ns = i;
            dif = dift;
        }
        c[i] = ya[i];
        d[i] = ya[i];
    }
    *y = ya[ns--];
    for (m = 1; m < n; m++) {
        for (i = 1; i <= n - m; i++) {
            ho = xa[i] - x;
            hp = xa[i + m] - x;
            w = c[i + 1] - d[i];
            if ((den = ho - hp) == 0.0)
                error("Error in routine polint");
            den = w / den;
            d[i] = hp * den;
            c[i] = ho * den;
        }
        *y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
    }
    free_vector(d, 1, n);
    free_vector(c, 1, n);
}
