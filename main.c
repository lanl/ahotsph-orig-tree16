#include <stdio.h>
#include <stdlib.h>

#include "ptw.h"

int main(int argc, char* argv[]) {
    plasticity_params_t params;
    params.theta = 0.1;
    params.p = 2.0;
    params.s0 = 0.02;
    params.sInf = 0.01;
    params.kappa = 0.3;
    params.lgamma = -12.0;
    params.y0 = 0.01;
    params.yInf = 0.003;
    params.y1 = 0.09;
    params.y2 = 0.7;

    state_t state;

    const int nhist = 100;
    double edot = 2500.0 * 1.e-6;
    double emax = 0.6;
    double tmax = emax / edot;
    double dt = tmax / (nhist - 1);
    const double edotcrit = 1.e-6;

    state.time = malloc(nhist * sizeof(double));
    state.stress = malloc(nhist * sizeof(double));
    state.strain_rate = malloc(nhist * sizeof(double));
    state.strain = malloc(nhist * sizeof(double));
    state.temp = malloc(nhist * sizeof(double));
    state.G = malloc(nhist * sizeof(double));
    state.rho = malloc(nhist * sizeof(double));
    state.Tmelt = malloc(sizeof(double));

    /* initialize state */
    state.temp[0] = 1000.0;
    state.strain[0] = 0.0;
    state.strain_rate[0] = edot;
    state.Tmelt[0] = consts.TMelt0;
    state.G[0] = calc_shear_modulus(&(state.temp[0]), &(state.Tmelt[0]));
    state.rho[0] = consts.rho0;
    state.stress[0] = ptw(
        &(state.strain_rate[0]), &(state.temp[0]), state.Tmelt, &(state.G[0]), &(state.strain[0]), &params);
    printf("%f %.10f %f %.10f %.10f %.10f %f\n",
           state.time[0],
           state.strain[0],
           state.strain_rate[0],
           state.stress[0],
           state.temp[0],
           state.G[0],
           state.rho[0]);

    for (int i = 1; i < nhist; i++) {
        state.time[i] = state.time[i - 1] + dt;
        double C_v = calc_specific_heat();
        state.rho[i] = consts.rho0;
        state.temp[i] = state.temp[i - 1];
        state.strain[i] = state.strain[i - 1] + edot * dt;
        state.strain_rate[i] = (state.strain[i] - state.strain[i - 1]) / dt;
        update_T(&(state.temp[i]),
                 &(state.stress[i - 1]),
                 &(state.strain_rate[i - 1]),
                 &dt,
                 &C_v,
                 &(state.rho[i-1]));
        state.Tmelt[i] = calc_tmelt();
        state.G[i] = calc_shear_modulus(&(state.temp[i]), &(state.Tmelt[i]));
        state.stress[i] = ptw(&(state.strain_rate[i]),
                              &(state.temp[i]),
                              &(state.Tmelt[i]),
                              &(state.G[i]),
                              &(state.strain[i]),
                              &params);
        printf("%f %.10f %f %.10f %.10f %.10f %f\n",
               state.time[i],
               state.strain[i],
               state.strain_rate[i],
               state.stress[i],
               state.temp[i],
               state.G[i],
               state.rho[i]);
    }
}