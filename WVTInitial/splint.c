/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

void splint(float xa[], float ya[], float y2a[], int n, float x, float *y) {
    void nrerror(char error_text[]);
    int klo, khi, k;
    float h, b, a;

    klo = 1;
    khi = n;
    while (khi - klo > 1) {
        k = (khi + klo) >> 1;
        if (xa[k] > x)
            khi = k;
        else
            klo = k;
    }
    h = xa[khi] - xa[klo];
    if (h == 0.0)
        nrerror("Bad xa input to routine splint");
    a = (xa[khi] - x) / h;
    b = (x - xa[klo]) / h;
    *y = a * ya[klo] + b * ya[khi]
         + ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
}
