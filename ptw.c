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
    // Stein Shear Modulus
    double aterm = 0.0;
    double bterm = consts.sgB * (*temp - 300.0);
    double gnow = consts.G0 * (1.0 + aterm - bterm);
    if (*temp > *tmelt) gnow = 0.0;
    if (gnow < 0.0) gnow = 0.0;
    return &gnow;
};
void* calc_flow_stress();
