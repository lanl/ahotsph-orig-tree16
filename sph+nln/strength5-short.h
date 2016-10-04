/* n.b., Fortran wants real*4 passed as double's, and all arguments passed by reference */
void plastic_(float *sxxi, float *syyi, float *sxyi, float *sxzi, float *syzi, float *ui, float *dmi, float *umelti, float *yiei, float *vonmises);

void straintensor_(double *grpmrj, double *dvx, double *dvy, double *dvz, double *dx, double *dy, double *dz, double *depsxxi, double *depsyyi, double *depszzi, double *depsxyi, double *depsxzi, double *depsyzi, double *drxyi, double *drxzi, double *dryzi);

void deviator_(float *xmui, float *sxxi, float *syyi, float *szzi, float *sxyi, float *sxzi, float *syzi, float *epsxxi, float *epsyyi, float *epszzi, float *epsxyi, float *epsxzi, float *epsyzi, float *rxyi, float *rxzi, float *ryzi, float *dsxxi, float *dsyyi, float *dsxyi, float *dsxzi, float *dsyzi);

void strengthforce_(float *grpmj, float *rhoij, float *sxxi, float *syyi, float *sxyi, float *sxzi, float *syzi, float *sxxj, float *syyj, float *sxyj, float *sxzj, float *syzj, float *dmi, float *dmj, float *dx, float *dy, float *dz, float *dfxi, float *dfyi, float *dfzi);

void fracture_(float *sxxi, float *syyi, float *sxyi, float *sxzi, float *syzi, float *pri, float *dmi, int *nflawi, int *ifrac, float *youngi, float *epsmini, float *xmi, float *acoefi, float *ddmi);

void strengthdu_(float *rhoi, float *sxxi, float *syyi, float *sxyi, float *sxzi, float *syzi, float *epsxxi, float *epsyyi, float *epsxyi, float *epsxzi, float *epsyzi, float *epszzi, float *dmi, float *dudt);
