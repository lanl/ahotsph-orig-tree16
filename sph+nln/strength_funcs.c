#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>

void init_defects_table(int gnobj, int Nflaws, double **eps, int ***flaws_tbl_lookup, float kVol, float m) {
	int i, j, nflaws_i;
	int part_id, cur_nflaws, index;
	int ran;
	double **eps_act;
	float iv_kVol, iv_m;

	/* initialize lookup table */
	nflaws_i = (int)(Nflaws/gnobj);
	if (nflaws_i * gnobj < Nflaws)
		nflaws_i += 1;
	nflaws_i *= 3; /* generally, the max-nflaws_i seems to stay at < 3* ave-nflaws_i */
	for (i = 0; i < gnobj; i++) {
		(*flaws_tbl_lookup)[i][1] = -1; /* starting array index in eps[] */
		(*flaws_tbl_lookup)[i][0] = 0; /* number of flaws in particle */
	}

	/* create a ton of space for all the eps_activation values */
	eps_act = (double **)malloc ( gnobj * sizeof (double *));
	for (i = 0; i < gnobj; i++) {
		eps_act[i] = (double *) malloc (nflaws_i * sizeof (double));
	}

	iv_kVol = 1./kVol;
	iv_m = 1./m;
	srand(time(NULL));

	/* create temporary table of activation thresholds */
	for (j = 1; j <= Nflaws; j++) {
		/* select random particle, add flaw */
		ran = rand();
		/* this is technically part_id-1 */
		part_id = (int) ((double) ran / (double) (RAND_MAX - 1) * gnobj);
		cur_nflaws = (*flaws_tbl_lookup)[part_id][0];
		if (cur_nflaws >= nflaws_i) {
			fprintf(stdout, "too many flaws in particle\n");
			exit(1);
		}
		(*flaws_tbl_lookup)[part_id][0] += 1;
		eps_act[part_id][cur_nflaws] = pow ((double) j * iv_kVol, iv_m);
	}

	/* create actual table with flaw activation thresholds */
	index = 0;
	for (i = 0; i < gnobj; i++) {
		cur_nflaws = (*flaws_tbl_lookup)[i][0];
		if (cur_nflaws > 0)
			(*flaws_tbl_lookup)[i][1] = index;
	//	memcpy (&((*eps)[index]), &(eps_act[i][0]), (size_t) cur_nflaws);
		if (index + cur_nflaws > Nflaws) {
			fprintf(stdout, "too many flaws\n");
			exit(1);
		}
		for (j = 0; j < cur_nflaws; j++)
			(*eps)[index + j] = eps_act[i][j];
		index += cur_nflaws;
	}

	for (i = 0; i < gnobj; i++)
		free(eps_act[i]);
	free(eps_act);
}
