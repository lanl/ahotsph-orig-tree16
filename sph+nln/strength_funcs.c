#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "SDF.h"
#include "malloc.h"
#include "bigmalloc.h"
#include "Msgs.h"

void init_defects_table(int gnobj, int Nflaws, double **eps, int **flaws_tbl_lookup, float kVol, float m) {
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
		(*flaws_tbl_lookup)[i * 2] = -1; /* starting array index in eps[] */
		(*flaws_tbl_lookup)[i * 2 + 1] = 0; /* number of flaws in particle */
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
		cur_nflaws = (*flaws_tbl_lookup)[part_id * 2 + 1];
		if (cur_nflaws >= nflaws_i) {
			fprintf(stdout, "too many flaws in particle\n");
			exit(1);
		}
		(*flaws_tbl_lookup)[part_id * 2 + 1] += 1;
		eps_act[part_id][cur_nflaws] = pow ((double) j * iv_kVol, iv_m);
	}

	/* create actual table with flaw activation thresholds */
	index = 0;
	for (i = 0; i < gnobj; i++) {
		cur_nflaws = (*flaws_tbl_lookup)[i * 2 + 1];
		if (cur_nflaws > 0)
			(*flaws_tbl_lookup)[2 * i] = index;
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

void read_defects_table(SDF *sdfp, int *nflaws, double **eps, int **flaws_tbl_lookup) {
	if (!sdfp)
		Error("Unable to access file with defects table\n");

	int npart, index, count;
	int i, j;
	int *part_ids;

	SDFgetintOrDie(sdfp, "npart", &npart);
	SDFgetintOrDie(sdfp, "nflaws", nflaws);

	(*eps) = (double *) malloc (*nflaws * sizeof (double));
	(*flaws_tbl_lookup) = (int *) malloc (2 * npart * sizeof (int));
	part_ids = (int *) malloc (*nflaws * sizeof (int));

	SDFrdvecs(sdfp, "part_id", nflaws, &part_ids, SDFtype_sizes[SDF_INT],
			"eps_actv", nflaws, eps, SDFtype_sizes[SDF_DOUBLE]);

	/* create the lookup table */
	index = 0;
	for (i = 0, j = 0; i < npart; i++) {
		(*flaws_tbl_lookup)[part_ids[j] * 2] = index;
		while (part_ids[j] == part_ids[j + 1]) {
			j++;
		}
		j++;
		(*flaws_tbl_lookup)[i * 2 + 1] = j - index;
		index = j;
	}

	free (part_ids);
}

void write_defects_table (char *name, int gnobj, int nflaws, double *eps, int *flaws_tbl_lookup) {
	FILE *fp = NULL;
	int i, j, nflaws_i;
	int index;

	fp = fopen (name, "w");
	
	/* write header */
    fprintf (fp, "# SDF\n");
    fprintf (fp, "parameter byteorder = %#x;\n", SDFcpubyteorder());
	fprintf (fp, "int npart = %d;\n", gnobj);
	fprintf (fp, "int nflaws = %d;\n", nflaws);
	fprintf (fp, "struct {\n");
	fprintf (fp, "\tint part_id;\n");
	fprintf (fp, "\tdouble eps_actv;\n");
    fprintf (fp, "}[%d];\n", nflaws);
    fprintf (fp, "#\n");
    fprintf (fp, "# SDF-EOH\n");

	/* write data */
	for (i = 0; i < gnobj; i++) {
		index = flaws_tbl_lookup[i * 2];
		nflaws_i = flaws_tbl_lookup[i * 2 + 1];
		for (j = 0; j < nflaws_i; j++){
			fwrite (&i, sizeof (int), 1, fp);
			fwrite (&eps[index + j], sizeof (double), 1, fp);
		}
	}

	fclose (fp);
}
