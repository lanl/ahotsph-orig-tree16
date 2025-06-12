#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "test_ptw.h"

bool check_ptw(const double* ptw_value, const double* edot, const double* temp, const double* shear, const double* eps, const plasticity_params_t* params, const mat_consts_t* consts) {
    double ptw_calc = ptw(edot, temp, shear, eps, params, consts);
    if (fabs(*ptw_value - ptw_calc) / (*ptw_value) > 1.0e-10) {
        printf("Error: ptw_test = %f, ptw_calc = %f\n", *ptw_value, ptw_calc);
        return false;
    }
    return true;
}

bool approx_equal(const double* expected, const double* actual, const double* rel_err) {
    return fabs(*expected - *actual) / *expected < *rel_err;
}

bool check_all(const state_t* state, const int n) {
    // open and read in test_physics.dat
    FILE* fp = fopen("../test_physics.dat", "r");

    if (fp ==  NULL) {
        printf("Error opening file\n");
        return false;
    }
    // determine number of lines in file
    int nlines = 0;
    char c;
    while ((c = fgetc(fp)) != EOF) {
        if (c == '\n') {
            nlines++;
        }
    }
    fseek(fp, 0, SEEK_SET);

    // subtract header line
    nlines--;

    if (n != nlines) {
        printf("Warning: Number of iterations not equal to verification data.\n");
        printf("%d vs %d\n", nlines, n);
        // return false;
    }

    double times[nlines], strain[nlines], strain_rate[nlines], stress[nlines], temp[nlines], G[nlines], rho[nlines];

    char header[255];
    fgets(header, sizeof(header), fp);
    printf("%s\n", header);

    for (int i = 0; i < nlines; ++i) {
        fscanf(fp, "%lf %lf %lf %lf %lf %lf %lf", &(times[i]), &(strain[i]), &(strain_rate[i]), &(stress[i]), &(temp[i]), &(G[i]), &(rho[i]));
        /* data files from Impala are in usec. */
        strain_rate[i] *= 1.0e6;
        times[i] *= 1.0e-6;
    }            

    double rel_err = 1.0e-4;
    for (int idx = 0; idx < nlines; ++idx) {
        if (! approx_equal(&(state->time[idx+1]), &(times[idx]), &rel_err)) 
            printf("Error: ptw test failed at iteration %i: time not equal. Expected: %f, actual: %f\n",
                idx, times[idx], state->time[idx+1]);
        if (! approx_equal(&(state->stress[idx+1]), &(stress[idx]), &rel_err)) 
            printf("Error: ptw test failed at iteration %i: stress not equal. Expected: %f, actual: %f\n",
                idx, stress[idx], state->stress[idx+1]);
        if (! approx_equal(&(state->temp[idx+1]), &(temp[idx]), &rel_err)) 
            printf("Error: ptw test failed at iteration %i: temperature not equal. Expected: %f, actual: %f\n",
                idx, temp[idx], state->temp[idx+1]);
        if (! approx_equal(&(state->G[idx+1]), &(G[idx]), &rel_err)) 
            printf("Error: ptw test failed at iteration %i: G not equal. Expected: %f, actual: %f\n",
                idx, G[idx], state->G[idx+1]);        
    }
}

bool check_shear_modulus(const double* temp, const mat_consts_t* consts) {
    double G_calc = calc_shear_modulus(temp, consts);
    double G_test = 0.0;
    if (*temp > consts->TMelt0) {
        G_test = 0.0;
    } else {
        G_test = consts->G0 * (1.0 + 0.0 - consts->sgB * (*temp - 300.0));
    }
    if (abs(G_test - G_calc) > 1.0e-10) {
        printf("Error: G_test = %f, G_calc = %f\n", G_test, G_calc);
        return false;
    }
    return true;
}

bool check_update_T(const double* stress, const double* edot, const double* dt, const double* C_v, const double* rho, const mat_consts_t* consts) {
    double temp_calc = 0.0;
    double temp_test = 0.0;
    double edotcrit = 1.0e-6;
    int cond = (int)(*edot > edotcrit);
    temp_calc = temp_test + (double)cond * consts->chi * (*stress) * (*edot) * (*dt) / ((*C_v) * (*rho));
    if (abs(temp_calc - temp_test) > 1.0e-10) {
        printf("Error: temp_calc = %f, temp_test = %f\n", temp_calc, temp_test);
        return false;
    }
    return true;
}