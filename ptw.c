#include "ptw.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

const consts_t consts = {.alpha = 0.2,
                         .beta = 0.33,
                         .mAtomic = 45.9,
                         .TMelt0 = 2110.0,
                         .rho0 = 4.419,
                         .Cv0 = 0.525e-5,
                         .G0 = 0.4,
                         .chi = 1.0,
                         .sgB = 6.44e-4};
params_t params;

double ptw(const double* edot,
           const double* temp,
           const double* tmelt,
           const double* shear,
           const double* eps) {
    double scaled_stress = -999.0;
    bool good = (params.sInf < params.s0) * (params.yInf < params.y0) * (params.y0 < params.s0)
                * (params.yInf < params.sInf) * (params.y1 > params.s0) * (params.y2 > consts.beta);
    if (!good) {
        printf("PTW bad val.");
        return scaled_stress;
    };

    // convert to 1/s strain rate since PTW rate is in that unit
    double edot_scaled = *edot * 1.0e6;
    double t_hom = (*temp) / (*tmelt);
    double afact = (4.0 / 3.0) * M_PI * consts.rho0 / consts.mAtomic;
    double ainv = pow(afact, (1.0 / 3.0));
    double xfact = sqrt(*shear / consts.rho0);
    /* LA-UR-04-0305 eqn. 3 */
    double xiDot = 0.5 * ainv * xfact * pow(6.022e29, (1.0 / 3.0)) * 1.0e4;
    double argErf = params.kappa * t_hom * (params.lgamma + log(xiDot / edot_scaled));
    /* LA-UR-04-0305 eqn. 7 */
    double saturation1 = params.s0 - (params.s0 - params.sInf) * erf(argErf);
    /* LA-UR-04-0305 eqn. 8 */
    double saturation2 = params.s0 * exp(consts.beta * (-params.lgamma + log(edot_scaled / xiDot)));
    double tau_s;
    /* LA-UR-04-0305 eqn. 9 */
    if (saturation1 > saturation2) {
        tau_s = saturation1;
    } else {
        tau_s = saturation2;
    };
    /* LA-UR-04-0305 eqn. 6 */
    double ayield = params.y0 - (params.y0 - params.yInf) * erf(argErf);
    /* LA-UR-04-0305 eqn. 10, but xiDot and edot_scaled are flipped? */
    double byield = params.y1 * exp(-params.y2 * (params.lgamma + log(xiDot / edot_scaled)));
    /* LA-UR-04-0305 eqn. 8, since at very high strain rates tauhat_y = tauhat_s */
    double cyield = params.s0 * exp(-consts.beta * (params.lgamma + log(xiDot / edot_scaled)));
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
    // int ind = (int)((params.p > small) * (fabs(tau_s - tau_y) > small));
    double eArg1 = (params.p * (tau_s - tau_y) / (params.s0 - tau_y));
    double eArg2 = ((*eps) * params.p * params.theta) / (params.s0 - tau_y) / (exp(eArg1) - 1.0);
    double check_val = 1.0 - (1.0 - exp(-eArg1)) * exp(-eArg2);
    if ((check_val <= 0.0) || (check_val != check_val))
        printf("bad\n");
    double theLog = log(check_val);
    scaled_stress = (tau_s + (params.s0 - tau_y) * theLog / params.p);
    // int ind2 = (int)((params.p <= small) * tau_s > tau_y);
    scaled_stress = (tau_s - (tau_s - tau_y) * exp(-(*eps) * params.theta / (tau_s - tau_y)));
    return (scaled_stress * (*shear) * 2.0);
};

double calc_specific_heat() { return (consts.Cv0); };

void update_T(double* temp,
              const double* stress,
              const double* edot,
              const double* dt,
              const double* C_v,
              const double* rho) {
    const double edotcrit = 1.e-6;
    int cond = (int)(*edot > edotcrit);
    *temp += (double)cond * consts.chi * (*stress) * (*edot) * (*dt) / ((*C_v) * (*rho));
};

double calc_tmelt() { return consts.TMelt0; };

double calc_shear_modulus(const double* temp, const double* tmelt) {
    // Stein Shear Modulus
    double aterm = 0.0;
    double bterm = consts.sgB * (*temp - 300.0);
    double gnow = consts.G0 * (1.0 + aterm - bterm);
    if (*temp > *tmelt)
        gnow = 0.0;
    if (gnow < 0.0)
        gnow = 0.0;
    return gnow;
};

void* calc_flow_stress();
