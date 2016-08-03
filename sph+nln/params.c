#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SDF.h"
#include "singlio.h"
#include "error.h"
#include "params.h"
#include "ndim.h"

void read_initial_ctl(SDF *sdfp, setup_params_t *params) {

	/* default values as found in main.c */
	params->setpvel = 0;

	/* read in settings */
    SDFgetintOrDefault(sdfp, "timeout", &(params->timeout), 600);
//    if (params->timeout > 0) MPMY_TimeoutSet(params->timeout);
#ifdef __PARAGON__
    {
	int fail_if_slow;
	SDFgetintOrDefault(sdfp, "fail_if_slow", &fail_if_slow, 0);
	chk_slow(fail_if_slow);
    }
#endif

    SDFgetstring(sdfp, "datafile", params->name, sizeof(params->name));
    SDFgetintOrDefault(sdfp, "do_restart", &(params->do_restart), 0);
    SDFgetintOrDefault(sdfp, "do_periodic", &(params->do_periodic), 0);
    SDFgetintOrDefault(sdfp, "cosmology", &(params->cosmology), 0);
    SDFgetintOrDefault(sdfp, "set_id", &(params->set_id), 0);
    SDFgetintOrDefault(sdfp, "setpvel", &(params->setpvel), 0);
    SDFgetintOrDefault(sdfp, "do_sph", &(params->do_sph), 1);
    SDFgetintOrDefault(sdfp, "do_diffusion", &(params->do_diffusion), 0);
    SDFgetintOrDefault(sdfp, "do_cooling", &(params->do_cooling), 0);
    SDFgetintOrDefault(sdfp, "do_burning", &(params->do_burning), 0);
    SDFgetintOrDefault(sdfp, "do_grav", &(params->do_grav), 0);
    SDFgetintOrDefault(sdfp, "do_winds", &(params->do_winds), 0);
    SDFgetintOrDefault(sdfp, "do_point_mass", &(params->do_point_mass), 0);
    SDFgetintOrDefault(sdfp, "do_point_mass2", &(params->do_point_mass2), 0);
    SDFgetintOrDefault(sdfp, "do_boundary", &(params->do_boundary), 0);
    SDFgetintOrDefault(sdfp, "do_absorbing_bndry", &(params->do_absorbing_bndry), 0);
    SDFgetintOrDefault(sdfp, "do_drag", &(params->do_drag), 0);
    SDFgetintOrDefault(sdfp, "has_grav_data", &(params->has_grav_data), params->do_grav);

    if (SDFhasname("SPHdatafile", sdfp))
        SDFgetstring(sdfp, "SPHdatafile", params->SPHdatafile, sizeof(params->SPHdatafile));

    SDFgetfloatOrDefault(sdfp, "new_h", &(params->new_h), 0.0);
    SDFgetfloatOrDefault(sdfp, "new_u", &(params->new_u), 0.0);

    if (params->do_point_mass) 
        SDFgetfloatOrDie(sdfp, "r_inner", &(params->r_inner));

    if (params->do_point_mass2) {
        SDFgetfloatOrDie(sdfp, "r_inner", &(params->r_inner));
        SDFgetfloatOrDie(sdfp, "centmass", &(params->centmass));
	}

	if (params->do_boundary) {
	    SDFgetfloatOrDie(sdfp, "r_inner", &(params->r_inner));
	    SDFgetfloatOrDie(sdfp, "r_outer", &(params->r_outer));
	    SDFgetfloatOrDie(sdfp, "centmass", &(params->centmass));
	}
}

void print_initial_ctl(setup_params_t params) {
    singlPrintf("printing params structure\n");

	singlPrintf("int timeout = %d;\n", params.timeout);
	singlPrintf("string datafile[] = %s;\n", params.name);
	singlPrintf("int do_restart = %d;\n", params.do_restart);
    singlPrintf("int do_periodic = %d;\n", params.do_periodic);
    singlPrintf("int cosmology = %d;\n", params.cosmology);
	singlPrintf("int set_id = %d;\n", params.set_id);
    singlPrintf("int setpvel = %d;\n", params.setpvel);
	singlPrintf("int do_sph = %d;\n", params.do_sph);
    singlPrintf("int do_diffusion = %d;\n", params.do_diffusion);
    singlPrintf("int do_cooling = %d;\n", params.do_cooling);
    singlPrintf("int do_burning = %d;\n", params.do_burning);
	singlPrintf("int do_grav = %d;\n", params.do_grav);
    singlPrintf("int do_winds = %d;\n", params.do_winds);
	singlPrintf("int do_point_mass = %d;\n", params.do_point_mass);
	singlPrintf("int do_point_mass2 = %d;\n", params.do_point_mass2);
	singlPrintf("int do_boundary = %d;\n", params.do_boundary);
	singlPrintf("int do_absorbing_bndry = %d;\n", params.do_absorbing_bndry);
    singlPrintf("int do_drag = %d;\n", params.do_drag);
	singlPrintf("int has_grav_data = %d;\n", params.has_grav_data);

	singlPrintf("string SPHdatafile[] = %s;\n", params.SPHdatafile);
    if (params.do_point_mass || params.do_point_mass2) {
        singlPrintf("float r_inner = %f;\n", params.r_inner);
        singlPrintf("float centmass = %e;\n", params.centmass);
    }

    if (params.do_boundary) {
        singlPrintf("float r_inner = %f;\n", params.r_inner);
        singlPrintf("float r_outer = %f;\n", params.r_outer);
        singlPrintf("float centmass = %e;\n", params.centmass);
    }

    singlPrintf("end printing params structure\n");
}

void print_absorb_bndry(bndry_t bndry) {
    if (params.do_absorbing_bndry) {
		/* print position */
		singlPrintf("float bndry[] = [ %g", bndry.pos[0]);
#if NDIM>=2
		singlPrintf(",%g ", bndry.pos[1]);
#if NDIM>=3
		singlPrintf(",%g ", bndry.pos[2]);
#endif
#endif
		singlPrintf("];\n");

		/* print velocity */
		singlPrintf("float bndry_vel[] = [ %g", bndry.vel[0]);
#if NDIM>=2
		singlPrintf(",%g ", bndry.vel[1]);
#if NDIM>=3
		singlPrintf(",%g ", bndry.vel[2]);
#endif
#endif
		singlPrintf("];\n");

		/* print linear momentum */
		singlPrintf("float bndry_p[] = [ %g", bndry.p[0]);
#if NDIM>=2
		singlPrintf(",%g ", bndry.p[1]);
#if NDIM>=3
		singlPrintf(",%g ", bndry.p[2]);
#endif
#endif
		singlPrintf("];\n");

		/* print angular momentum */
		singlPrintf("float bndry_l[] = [ %g", bndry.l[0]);
#if NDIM>=2
		singlPrintf(",%g ", bndry.l[1]);
#if NDIM>=3
		singlPrintf(",%g ", bndry.l[2]);
#endif
#endif
		singlPrintf("];\n");

        singlPrintf("float bndry_mass = %g;\n", bndry.mass);
        singlPrintf("float bndry_r = %g;\n", bndry.r);
    }
}

void read_absorb_bndry(SDF *sdfp, bndry_t *bndry) {
    SDFgetfloatOrDie(sdfp, "bndry_x", &(bndry->pos[0]));
    SDFgetfloatOrDie(sdfp, "bndry_vx", &(bndry->vel[0]));
    SDFgetfloatOrDefault(sdfp, "bndry_lx", &(bndry->l[0]), 0.0);
    SDFgetfloatOrDefault(sdfp, "bndry_px", &(bndry->p[0]), 0.0);
#if NDIM>=2        
    SDFgetfloatOrDie(sdfp, "bndry_y", &(bndry->pos[1]));
    SDFgetfloatOrDie(sdfp, "bndry_vy", &(bndry->vel[1]));
    SDFgetfloatOrDefault(sdfp, "bndry_ly", &(bndry->l[1]), 0.0);
    SDFgetfloatOrDefault(sdfp, "bndry_py", &(bndry->p[1]), 0.0);
#if NDIM>=3
    SDFgetfloatOrDie(sdfp, "bndry_z", &(bndry->pos[2]));
    SDFgetfloatOrDie(sdfp, "bndry_vz", &(bndry->vel[2]));
    SDFgetfloatOrDefault(sdfp, "bndry_lz", &(bndry->l[2]), 0.0);
    SDFgetfloatOrDefault(sdfp, "bndry_pz", &(bndry->p[2]), 0.0);
#endif
#endif
    SDFgetfloatOrDie(sdfp, "bndry_mass", &(bndry->mass));
    SDFgetfloatOrDie(sdfp, "bndry_r", &(bndry->r));

}
