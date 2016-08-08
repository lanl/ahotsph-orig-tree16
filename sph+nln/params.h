#include "SDF.h"
#include "ndim.h"

typedef struct {
    char name[256]; /* "datafile" */
    char SPHdatafile[256]; /* "SPHdatafile" */
	char template_name[256]; /* wind template */
	char winddata_name[256]; /* "winddata_name" */
    int timeout;
    int fail_if_slow;
    int do_restart;
    int do_periodic;
    int cosmology;
    int set_id;
    int setpvel;
    int do_sph;
    int do_diffusion;
    int do_cooling;
    int do_burning;
    int do_grav;
    int do_winds;
    int do_point_mass;
    int do_point_mass2;
    int do_boundary;
    int do_absorbing_bndry;
    int do_drag;
    int has_grav_data;
    int windpartpershell;
    int old_winds;
    int const_winds;
    int nonconst_winds;
    int accreting_winds;
    int do_DL;
    int do_BH;
    int do_Bmax;
    int do_Arel;
    int do_output;
    int nsteps;
    int log_time;
    int comov_eps;
    int save_first;
    int ntimer_detail;
    int exact_rho;
    int nbrcut_max;
    int nbrcut_min;
    int adaptive_dt;
    int independent_dt;
    int dark_independent_dt;
    int default_nterms;
    float new_h;
    float new_u;
    float r_inner;
    float r_outer;
    float centmass;
    float openangle_wind;
    float omega_wind;
    float r_wind;
    float t_wind;
    float v_wind;
    float mdot_wind;
    float u_wind;
    float eps;
    float tol;
    float frac_tol;
    float CWfac;
    float SPHCWfac;
	float dt;
    float dark_dt;
    float comov_eps_epoch;
    float visc_alpha;
    float visc_beta;
    float visc_epsilon;
    float heat_f1;
    float min_h;
    float max_h;
    float nbrcut_fac;
    float Gamma;
    float courant_number;
    float fmassCF;
    float flenCF;
    float ftimeCF;
} setup_params_t;

typedef struct{
    float pos[NDIM];
    float vel[NDIM];
    float p[NDIM];
    float l[NDIM];
    float mass;
    float r;
} bndry_t;

extern setup_params_t params;

void read_initial_ctl (SDF *sdfp, setup_params_t *params);
void print_initial_ctl (setup_params_t params);
void read_absorb_bndry (SDF *sdfp, bndry_t *bndry);
void print_absorb_bndry (bndry_t bndry);
