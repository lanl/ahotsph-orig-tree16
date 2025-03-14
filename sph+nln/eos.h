
extern double eos_n;
extern double eos_u;

#ifndef EOS_H
#define EOS_H
typedef struct material_s {
    double Vol0;      /* initial target volume */
    double rho0;      /* original density */
    double A;         /* bulk modulus, units of pressure */
    double B;         /* Tillotson parameter, units of pressure */
    double a;         /* Tillotson parameter */
    double b;         /* Tillotson parameter */
    double alpha;     /* Tillotson parameter */
    double beta;      /* Tillotson parameter */
    double Eiv;       /* energy of incipient vaporization */
    double Ecv;       /* energy of complete vaporization */
    double u0;        /* initial material specific energy density */
    double umelt;     /* melt energy */
    double mu;        /* strain, rho/rho0 - 1 */
    double yield;     /* yield strength */
    double pweib;     /* some Weibull parameter */
    double cweib;     /* some Weibull parameter */
    double G_shear;   /* shear modulus */
    double E_Young;   /* Young's modulus */
    float material_k; /* material coefficient for flaws */
    float material_m; /* powerlaw index for flaws */
} Material_t;
#endif
