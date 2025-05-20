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
    const double alpha;
    const double beta;
    const double mAtomic;
    const double TMelt0;
    const double rho0;
    const double Cv0;
    const double G0;
    const double chi;
    const double sgB;
} consts_t;

typedef struct state_s {
    double* stress;
    double* strain;
    double* strain_rate;
    double* temp;
    double* G;
    double* rho;
    double* time;
    double* Tmelt;
} state_t;

extern const consts_t mat_consts;
extern params_t ptw_params;

double ptw(const double* edot,
           const double* temp,
           const double* tmelt,
           const double* shear,
           const double* eps);
double calc_specific_heat();
void update_T(double* temp,
              const double* stress,
              const double* edot,
              const double* dt,
              const double* C_v,
              const double* rho);
double calc_tmelt();
double calc_shear_modulus(const double* temp, const double* tmelt);
void* calc_flow_stress();