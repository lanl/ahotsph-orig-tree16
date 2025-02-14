#include "ptw.h"

extern consts_t consts;
extern params_t params;

void* ptw(double* dt);
void* calc_specific_heat();
void* update_T();
void* calc_tmelt();
void* calc_shear_modulus();
void* calc_flow_stress();
