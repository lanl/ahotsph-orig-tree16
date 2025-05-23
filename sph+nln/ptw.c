#include "ptw.h"
#include "strength.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

double ptw(const double* edot,
           const double* temp,
           const double* shear,
           const double* eps,
           const plasticity_params_t* ptw_params,
           const Material_t* mat_consts) {
    double bad_scaled_stress = -999.0;
    double scaled_stress = bad_scaled_stress;
    bool good = (ptw_params->sInf < ptw_params->s0) * (ptw_params->yInf < ptw_params->y0) * (ptw_params->y0 < ptw_params->s0)
                * (ptw_params->yInf < ptw_params->sInf) * (ptw_params->y1 > ptw_params->s0) * (ptw_params->y2 > mat_consts->ptw_beta);
    if (!good) {
        printf("PTW bad val.");
        return bad_scaled_stress;
    };

    // convert to 1/s strain rate since PTW rate is in that unit
    double edot_scaled = *edot * 1.0e6;
    double t_hom = (*temp) / (mat_consts->tmelt);
    double afact = (4.0 / 3.0) * M_PI * mat_consts->rho0 / mat_consts->mAtomic;
    double ainv = pow(afact, (1.0 / 3.0));
    double xfact = sqrt(*shear / mat_consts->rho0);
    double xiDot = 0.5 * ainv * xfact * pow(6.022e29, (1.0 / 3.0)) * 1.0e4;
    double argErf = ptw_params->kappa * t_hom * (ptw_params->lgamma + log(xiDot / edot_scaled));
    double saturation1 = ptw_params->s0 - (ptw_params->s0 - ptw_params->sInf) * erf(argErf);
    double saturation2 = ptw_params->s0 * exp(mat_consts->ptw_beta * (-ptw_params->lgamma + log(edot_scaled / xiDot)));
    double tau_s;
    if (saturation1 > saturation2) {
        tau_s = saturation1;
    } else {
        tau_s = saturation2;
    };
    double ayield = ptw_params->y0 - (ptw_params->y0 - ptw_params->yInf) * erf(argErf);
    double byield = ptw_params->y1 * exp(-ptw_params->y2 * (ptw_params->lgamma + log(xiDot / edot_scaled)));
    double cyield = ptw_params->s0 * exp(-mat_consts->ptw_beta * (ptw_params->lgamma + log(xiDot / edot_scaled)));
    double dyield;
    if (byield < cyield) {
        dyield = byield;
    } else {
        dyield = cyield;
    };
    double tau_y;
    if (ayield > dyield) {
        tau_y = ayield;
    } else {
        tau_y = dyield;
    };
    const double small = 1.0e-10;
    scaled_stress = tau_s;
    //?? - use the commented if block for this
    int ind = (int)((ptw_params->p > small) * (fabs(tau_s - tau_y) > small));
    double eArg1 = (ptw_params->p * (tau_s - tau_y) / (ptw_params->s0 - tau_y));
    double eArg2 = ((*eps) * ptw_params->p * ptw_params->theta) / (ptw_params->s0 - tau_y) / (exp(eArg1) - 1.0);
    double check_val = 1.0 - (1.0 - exp(-eArg1)) * exp(-eArg2);
    if ((check_val <= 0.0) || (check_val != check_val)) {
        printf("bad\n");
        return bad_scaled_stress;
    }
    double theLog = log(check_val);
    scaled_stress = (tau_s + (ptw_params->s0 - tau_y) * theLog / ptw_params->p);
    int ind2 = (int)((ptw_params->p <= small) * tau_s > tau_y);
    scaled_stress = (tau_s - (tau_s - tau_y) * exp(-(*eps) * ptw_params->theta / (tau_s - tau_y)));
    return (scaled_stress * (*shear) * 2.0);
};
double calc_specific_heat(const Material_t* mat_consts) { return (mat_consts->C_v); };
void update_T(double* temp,
              const double* stress,
              const double* edot,
              const double* dt,
              const double* C_v,
              const double* rho,
              const double *chi) {
    const double edotcrit = 1.e-6;
    int cond = (int)(*edot > edotcrit);
    *temp += (double)cond * (*chi) * (*stress) * (*edot) * (*dt) / ((*C_v) * (*rho));
};
double calc_tmelt(const Material_t* mat_consts) { return mat_consts->tmelt; };
double calc_shear_modulus(const double* temp, const double* tmelt, const Material_t* mat_consts) {
    // Stein Shear Modulus
    double aterm = 0.0;
    double bterm = mat_consts->sgB * (*temp - 300.0);
    double gnow = mat_consts->G_shear * (1.0 + aterm - bterm);
    if (*temp > *tmelt)
        gnow = 0.0;
    if (gnow < 0.0)
        gnow = 0.0;
    return gnow;
};
void* calc_flow_stress();
