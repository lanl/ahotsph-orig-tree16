#pragma once

typedef struct params_s {
    double theta;
    double p;
    double s0;
    double sInf;
    double kappa;
    double lgamma;
    double y0;
    double yInf;
    double y1;
    double y2;
} params_t;

typedef struct consts_s {
    const double alpha = 0.2;
    const double beta = 0.33;
    const double mAtomic = 45.9;
    const double TMelt0 = 2110.0;
    const double rho0 = 4.419;
    const double Cv0 = 0.525e-5;
    const double G0 = 0.4;
    const double chi = 1.0;
    const double sgB = 6.44e-4;
} consts_t;

consts_t consts;
params_t params;

void* ptw(double* dt);
double* calc_specific_heat();
void *update_T(
    double *temp,
    const double *stress,
    const double *edot,
    const double *dt,
    const double *C_v,
    const double *rho);
void* calc_tmelt();
double* calc_shear_modulus(const double *temp, const double *tmelt);
void* calc_flow_stress();