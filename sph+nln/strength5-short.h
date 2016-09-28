/* n.b., Fortran wants real*4 passed as double's, and all arguments passed by reference */
void plastic_(double *sxxi, double *syyi, double *sxyi, double *sxzi, double *syzi, double *ui, double *dmi, double *umelti, double *yiei, double *vonmises);

void straintensor_(double *grpmrj, double *dvx, double *dvy, double *dvz, double *dx, double *dy, double *dz, double *depsxxi, double *depsyyi, double *depszzi, double *depsxyi, double *depsxzi, double *depsyzi, double *drxyi, double *drxzi, double *dryzi);

void deviator_(double *xmui, double *sxxi, double *syyi, double *szzi, double *sxyi, double *sxzi, double *syzi, double *epsxxi, double *epsyyi, double *epszzi, double *epsxyi, double *epsxzi, double *epsyzi, double *rxyi, double *rxzi, double *ryzi, double *dsxxi, double *dsyyi, double *dsxyi, double *dsxzi, double *dsyzi);

void strengthforce_(double *grpmj, double *rhoij, double *sxxi, double *syyi, double *sxyi, double *sxzi, double *syzi, double *sxxj, double *syyj, double *sxyj, double *sxzj, double *syzj, double *dmi, double *dmj, double *dx, double *dy, double *dz, double *dfxi, double *dfyi, double *dfzi);

void fracture_(double *sxxi, double *syyi, double *sxyi, double *sxzi, double *syzi, double *pri, double *dmi, int *nflawi, int *ifrac, double *youngi, double *epsmini, double *xmi, double *acoefi, double *ddmi);

void strengthdu_(double *rhoi, double *sxxi, double *syyi, double *sxyi, double *sxzi, double *syzi, double *epsxxi, double *epsyyi, double *epsxyi, double *epsxzi, double *epsyzi, double *epszzi, double *dmi, double *dudt);
