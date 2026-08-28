#include "ptw.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

double ptw(const double* edot,
           const double* temp,
           const double* shear,
           const double* eps,
           const plasticity_params_t* ptw_params,
           const mat_consts_t* mat_consts) {
    double scaled_stress = -999.0;
    bool good = (ptw_params->sInf < ptw_params->s0) * (ptw_params->yInf < ptw_params->y0)
                * (ptw_params->y0 < ptw_params->s0) * (ptw_params->yInf < ptw_params->sInf)
                * (ptw_params->y1 > ptw_params->s0) * (ptw_params->y2 > ptw_params->beta);
    if (!good) {
        printf("PTW bad val.");
        return scaled_stress;
    };

    // convert to 1/s strain rate since PTW rate is in that unit - N.b.: for Flag data
    double edot_scaled = *edot;  // * 1.0e6;
    double t_hom = (*temp) / (mat_consts->TMelt0);
    double afact = (4.0 / 3.0) * M_PI * mat_consts->rho0 / mat_consts->mAtomic;
    double ainv = pow(afact, (1.0 / 3.0));
    double xfact = sqrt(*shear / mat_consts->rho0);
    /* LA-UR-04-0305 eqn. 3 */
    double xiDot = 0.5 * ainv * xfact * pow(6.022e29, (1.0 / 3.0)) * 1.0e4;
    double argErf = ptw_params->kappa * t_hom * (ptw_params->lgamma + log(xiDot / edot_scaled));
    /* LA-UR-04-0305 eqn. 7 */
    double saturation1 = ptw_params->s0 - (ptw_params->s0 - ptw_params->sInf) * erf(argErf);
    /* LA-UR-04-0305 eqn. 8 */
    double saturation2
        = ptw_params->s0 * exp(ptw_params->beta * (-ptw_params->lgamma + log(edot_scaled / xiDot)));
    double tau_s;
    /* LA-UR-04-0305 eqn. 9 */
    if (saturation1 > saturation2) {
        tau_s = saturation1;
    } else {
        tau_s = saturation2;
    };
    /* LA-UR-04-0305 eqn. 6 */
    double ayield = ptw_params->y0 - (ptw_params->y0 - ptw_params->yInf) * erf(argErf);
    /* LA-UR-04-0305 eqn. 10, but xiDot and edot_scaled are flipped? */
    double byield
        = ptw_params->y1 * exp(-ptw_params->y2 * (ptw_params->lgamma + log(xiDot / edot_scaled)));
    /* LA-UR-04-0305 eqn. 8, since at very high strain rates tauhat_y = tauhat_s */
    double cyield
        = ptw_params->s0 * exp(-ptw_params->beta * (ptw_params->lgamma + log(xiDot / edot_scaled)));
    double dyield;
    /* LA-UR-04-0305 eqn. 11 */
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
    // int ind = (int)((ptw_params->p > small) * (fabs(tau_s - tau_y) > small));
    if (ptw_params->p > 0.0) {
        if (fabs(tau_s - tau_y) < small) {
            scaled_stress = tau_s;
        } else {
            double eArg1 = (ptw_params->p * (tau_s - tau_y) / (ptw_params->s0 - tau_y));
            double eArg2 = ((*eps) * ptw_params->p * ptw_params->theta) / (ptw_params->s0 - tau_y)
                           / (exp(eArg1) - 1.0);
            double theLog = log(1.0 - (1.0 - exp(-eArg1)) * exp(-eArg2));
            scaled_stress = (tau_s + (ptw_params->s0 - tau_y) * theLog / ptw_params->p);
        }
    } else {
        if (tau_s > tau_y) {
            scaled_stress
                = (tau_s - (tau_s - tau_y) * exp(-(*eps) * ptw_params->theta / (tau_s - tau_y)));
        } else {
            scaled_stress = tau_s;
        }
    }
    return (scaled_stress * (*shear) * 2.0);
};

double calc_specific_heat(const mat_consts_t* consts) { return (consts->Cv0); };

void update_T(double* temp,
              const double* stress,
              const double* edot,
              const double* dt,
              const double* C_v,
              const double* rho,
              const mat_consts_t* consts) {
    const double edotcrit = 1.0;
    int cond = (int)(*edot > edotcrit);
    *temp += (double)cond * consts->chi * (*stress) * (*edot) * (*dt) / ((*C_v) * (*rho));
};

double calc_tmelt(const mat_consts_t* consts) { return consts->TMelt0; };

double calc_shear_modulus(const double* temp, const mat_consts_t* consts) {
    // Stein Shear Modulus
    double aterm = 0.0;
    double bterm = consts->sgB * (*temp - 300.0);
    double gnow = consts->G0 * (1.0 + aterm - bterm);
    if (*temp > consts->TMelt0)
        gnow = 0.0;
    if (gnow < 0.0)
        gnow = 0.0;
    return gnow;
};

void* calc_flow_stress();
