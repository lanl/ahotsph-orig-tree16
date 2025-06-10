#pragma once

#include "strength.h"

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
    double alpha;
    double beta;
} plasticity_params_t;

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

double ptw(const double* edot,
           const double* temp,
           const double* shear,
           const double* eps,
           const plasticity_params_t* ptw_params,
           const Material_t* mat_consts);
double calc_specific_heat(const Material_t* mat_consts);
void update_T(double* temp,
              const double* stress,
              const double* edot,
              const double* dt,
              const double* C_v,
              const double* rho,
              const double* chi);
double calc_tmelt(const Material_t* mat_consts);
double calc_shear_modulus(const double* temp, const double* tmelt, const Material_t* mat_consts);
void* calc_flow_stress();