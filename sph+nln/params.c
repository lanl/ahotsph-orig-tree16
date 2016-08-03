#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "SDF.h"
#include "singlio.h"
#include "error.h"
#include "params.h"

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

    singlPrintf("end printing params structure\n");
}
