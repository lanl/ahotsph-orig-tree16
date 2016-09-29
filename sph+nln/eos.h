
extern double eos_n;
extern double eos_u;

typedef struct material_s {
	float rho0;
	float A, B, a, b, alpha, beta;
	float Eiv, Ecv;
	float u0, umelt;
	float mu, yield, pweib, cweib;
} Material_t;

/* In eos.c */
double uvst(double t);
double duvst(double t);
double liquid_eos (double k_bulk, double eta);
double murnaghan_eos(double k_bulk, double n_M, double eta);
void setconst1(Material_t *m);
void tillotson_eos (float rho, float u, Material_t *m, float *pressure, float *cs);
double anton_schmidt_eos(double k_bulk, double power_n, double eta);

