#include "ptw.h"

extern consts_t consts;
extern params_t params;

void* ptw(double* dt);
double* calc_specific_heat()
{
    return &(consts.Cv0);
};
void update_T(
    double* temp, 
    const double *stress, 
    const double* edot, 
    const double* dt,
    const double* C_v,
    const double* rho
)
{
    const double edotcrit = 1.e-6;
    int cond = (int)(*edot > edotcrit);
    *temp += (double)cond * consts.chi * (*stress) * (*edot) * (*dt) / ((*C_v) * (*rho));
};
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
