#pragma once

#include <stdbool.h>
#include "ptw.h"

bool check_all(const state_t* state, const int n);
bool check_ptw(const double* ptw_value, const double* edot, const double* temp, const double* tmelt, const double* shear, const double* eps, const plasticity_params_t* params);
bool check_shear_modulus(const double* temp, const double* tmelt);
bool check_update_T(const double* stress, const double* edot, const double* dt, const double* C_v, const double* rho);
