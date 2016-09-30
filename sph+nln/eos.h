
extern double eos_n;
extern double eos_u;

#ifndef EOS_H
#define EOS_H
typedef struct material_s {
	float rho0;
	float A, B, a, b, alpha, beta;
	float Eiv, Ecv;
	float u0, umelt;
	float mu, yield, pweib, cweib;
} Material_t;
#endif
