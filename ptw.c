#include "ptw.h"

extern consts_t consts;
extern params_t params;

void* ptw(double* dt);
void* calc_specific_heat();
void* update_T();
void calc_tmelt()
{
    return consts.TMelt0;
};
double* calc_shear_modulus(const double* temp, const double* tmelt)
{
    double shear = consts.G0 * (1.0 - consts.alpha * (*temp / *tmelt));
    return &shear;
};
void* calc_flow_stress();
