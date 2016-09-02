#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

void init_defects_table(int gnobj, int Nflaws, float **eps, int ***flaws_tbl_lookup, float kVol, float m) {
	int i, j, nflaws_i;
	int part_id, cur_nflaws, index;
	int ran;
	float **eps_act;
	float iv_kVol, iv_m;

	/* initialize lookup table */
	nflaws_i = (int)(Nflaws/gnobj);
	if (nflaws_i * gnobj < Nflaws)
		nflaws_i += 1;
	for (i = 0; i < gnobj; i++) {
		(*flaws_tbl_lookup)[i][1] = i * nflaws_i; /* starting array index in eps[] */
		(*flaws_tbl_lookup)[i][0] = 0; /* number of flaws in particle */
	}

	/* create a ton of space for all the eps_activation values */
	eps_act = (float **)malloc ( gnobj * sizeof (float *));
	for (i = 0; i < gnobj; i++) {
		eps_act[i] = (float *) malloc (3 * nflaws_i * sizeof (float));
		for (j = 0; j < 3 * nflaws_i; j++)
			eps_act[i][j] = -1.;
		/*
		*/
	}

	iv_kVol = 1./kVol;
	iv_m = 1./m;
	srand(time(NULL));
	fprintf(stdout, "RAND_MAX = %d\n", RAND_MAX);

	/* create temporary table of activation thresholds */
	for (j = 1; j <= Nflaws; j++) {
		/* select random particle, add flaw */
		ran = rand();
		part_id = (int) ((float) ran / (float) (RAND_MAX - 1) * gnobj);
		cur_nflaws = (*flaws_tbl_lookup)[part_id][0];
		(*flaws_tbl_lookup)[part_id][0] += 1;
		eps_act[part_id][cur_nflaws] = pow ((float) j * iv_kVol, iv_m);
	}

	/* create actual table with flaw activation thresholds */
	index = 0;
	for (i = 0; i < gnobj; i++) {
		cur_nflaws = (*flaws_tbl_lookup)[i][0];
		(*flaws_tbl_lookup)[i][1] = index;
	//	memcpy (&((*eps)[index]), &(eps_act[i][0]), (size_t) cur_nflaws);
		for (j = 0; j < cur_nflaws; j++)
			(*eps)[index + j] = eps_act[i][j];
		index += cur_nflaws;
	}

	for (i = 0; i < gnobj; i++)
		free(eps_act[i]);
	free(eps_act);
}
