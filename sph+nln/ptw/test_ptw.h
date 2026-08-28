/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#pragma once

#include <stdbool.h>

#include "ptw.h"

bool check_all(const state_t* state, const int n);
bool check_ptw(const double* ptw_value,
               const double* edot,
               const double* temp,
               const double* shear,
               const double* eps,
               const plasticity_params_t* params,
               const mat_consts_t* consts);
bool check_shear_modulus(const double* temp, const mat_consts_t* consts);
bool check_update_T(const double* stress,
                    const double* edot,
                    const double* dt,
                    const double* C_v,
                    const double* rho,
                    const mat_consts_t* consts);
