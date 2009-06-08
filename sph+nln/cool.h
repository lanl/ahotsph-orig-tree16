//void calc_lcool(double temp, double lcool)
//double calc_lcool(double temp,double uint,double dens,int numu,int numD);
double calc_lcool1(double temp,int extrapolate);
//double calc_lcool2(double temp,double uint,double dens,int numu,int numD);

double analytic_cool(double temp);

void locate(float xx[], long Nel, float x, long *j);

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

void polin2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y, double *dy);

void init_CoolTable(void);

//void init_ionfrac(void);

