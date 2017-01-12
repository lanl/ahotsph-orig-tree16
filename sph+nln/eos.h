
extern double eos_n;
extern double eos_u;

#ifndef EOS_H
#define EOS_H
typedef struct material_s {
	double rho0;
	double A, B, a, b, alpha, beta;
	double Eiv, Ecv;
	double u0, umelt;
	double mu, yield, pweib, cweib;
} Material_t;
#endif
