/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

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
    double alpha;
    double beta;
} plasticity_params_t;

typedef struct consts_s {
    const double mAtomic;
    const double TMelt0;
    const double rho0;
    const double Cv0;
    const double G0;
    const double chi;
    const double sgB;
} mat_consts_t;

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
           const mat_consts_t* mat_consts);
double calc_specific_heat(const mat_consts_t* mat_consts);
void update_T(double* temp,
              const double* stress,
              const double* edot,
              const double* dt,
              const double* C_v,
              const double* rho,
              const mat_consts_t* mat_consts);
double calc_tmelt(const mat_consts_t* mat_consts);
double calc_shear_modulus(const double* temp, const mat_consts_t* consts);
void* calc_flow_stress();