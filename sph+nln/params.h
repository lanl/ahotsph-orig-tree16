#include "SDF.h"

typedef struct {
    char name[256]; /* "datafile" */
    char SPHdatafile[256]; /* "SPHdatafile" */
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
    float new_h;
    float new_u;
    float r_inner;
    float r_outer;
    float centmass;
} setup_params_t;

extern setup_params_t params;

void read_initial_ctl(SDF *sdfp, setup_params_t *params);
void print_initial_ctl(setup_params_t params);
