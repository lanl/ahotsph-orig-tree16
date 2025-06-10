#ifndef STRENGTH_DOT_H
#define STRENGTH_DOT_H

#include "SDF.h"

/*#include "physics_sph.h"*/

extern double *flaw_actv_tbl;
extern double vol_scaling;
extern int *flaw_actv_tbl_lookup;

void init_defects_table(
    int gnobj, int Nflaws, double **eps, int **flaws_tbl_lookup, float kVol, float m);
void read_defects_table(SDF *sdfp, int *nflaws, double **eps, int **flaws_tbl_lookup);
void write_defects_table(char *name, int gnobj, int nflaws, double *eps, int *flaws_tbl_lookup);
// int has_strength(const SPHbody *p);
// void strength_force(double *grpmj,
//                     double *rhoij,
//                     double *sxxi,
//                     double *syyi,
//                     double *sxyi,
//                     double *sxzi,
//                     double *syzi,
//                     double *sxxj,
//                     double *syyj,
//                     double *sxyj,
//                     double *sxzj,
//                     double *syzj,
//                     double *dmi,
//                     double *dmj,
//                     double *dx,
//                     double *dy,
//                     double *dz,
//                     double *dfxi,
//                     double *dfyi,
//                     double *dfzi);

double equiv_strain(float const s[]);
double equiv_stress(float const s[]);

typedef struct material_s {
    double Vol0;      /* initial target volume */
    double rho0;      /* original density */
    double mAtomic;   /* atomic mass */
    double A;         /* bulk modulus, units of pressure */
    double B;         /* Tillotson parameter, units of pressure */
    double a;         /* Tillotson parameter */
    double b;         /* Tillotson parameter */
    double alpha;     /* Tillotson parameter */
    double beta;      /* Tillotson parameter */
    double Eiv;       /* energy of incipient vaporization */
    double Ecv;       /* energy of complete vaporization */
    double C_v;       /* specific heat capacity */
    double u0;        /* initial material specific energy density */
    double umelt;     /* melt energy */
    double tmelt;     /* melting temperature */
    double mu;        /* strain, rho/rho0 - 1 */
    double yield;     /* yield strength */
    double pweib;     /* some Weibull parameter */
    double cweib;     /* some Weibull parameter */
    double G_shear;   /* shear modulus */
    double E_Young;   /* Young's modulus */
    double chi;       /* parameter for PTW */
    double sgB;       /* parameter for PTW */
    double ptw_alpha; /* parameter for PTW */
    double ptw_beta;  /* parameter for PTW */
    float material_k; /* material coefficient for flaws */
    float material_m; /* powerlaw index for flaws */
} Material_t;
#endif /*STRENGTH_DOT_H*/