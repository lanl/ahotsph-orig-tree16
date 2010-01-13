/*
static const int NISO = 5;

typedef struct{
	float abund[NISO];
	int np[NISO];
	int nn[NISO];
} SPHbody;
*/

double calc_lcool1(float abundarr[],int nparr[],int nnarr[],double rho,double temp,int Gridpts,int Nel,int extrapolate);
double calc_lcool2(double temp,double uint,double dens,int numu,int numD);

double analytic_cool(double temp);

void locate(float xx[], long Nel, float x, long *j);

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

void polin2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, double x2, double *y, double *dy);

void init_CoolTable(int *Gridpts, int *Nel);
