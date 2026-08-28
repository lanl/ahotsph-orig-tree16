/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

 extern double *flaw_actv_tbl;
extern double vol_scaling;
extern int *flaw_actv_tbl_lookup;

void init_defects_table(
    int gnobj, int Nflaws, double **eps, int **flaws_tbl_lookup, float kVol, float m);
void read_defects_table(SDF *sdfp, int *nflaws, double **eps, int **flaws_tbl_lookup);
void write_defects_table(char *name, int gnobj, int nflaws, double *eps, int *flaws_tbl_lookup);
int has_strength(SPHbody p);
void strength_force(double *grpmj,
                    double *rhoij,
                    double *sxxi,
                    double *syyi,
                    double *sxyi,
                    double *sxzi,
                    double *syzi,
                    double *sxxj,
                    double *syyj,
                    double *sxyj,
                    double *sxzj,
                    double *syzj,
                    double *dmi,
                    double *dmj,
                    double *dx,
                    double *dy,
                    double *dz,
                    double *dfxi,
                    double *dfyi,
                    double *dfzi);
