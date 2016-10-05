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
#include "physics_sph.h"

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
		eps_act[part_id][cur_nflaws] = pow (((double) j * iv_kVol), iv_m);
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
	int part_id;
	void **addrs;
	void *data; /* need continuous space in memory for all the data */
	int stride = sizeof (int) + sizeof (double);
	int offset = sizeof (int);

	SDFgetintOrDie(sdfp, "npart", &npart);
	SDFgetintOrDie(sdfp, "nflaws", nflaws);

	(*eps) = (double *) malloc (*nflaws * sizeof (double));
	(*flaws_tbl_lookup) = (int *) malloc (2 * npart * sizeof (int));
	addrs = (void **) malloc (2 * sizeof (void *));
	data = (void *) malloc (*nflaws * stride);

	addrs[0] = (char *)data;
	addrs[1] = (char *)data + sizeof (int);

	count = SDFseekrdvecs(sdfp,
			"part_id", 0, *nflaws, addrs[0], stride,
			"eps_actv", 0, *nflaws, addrs[1], stride, 
			NULL);

	/* create the lookup table */
	index = 0;
	for (i = 0, j = 0; i < npart, j < *nflaws; i++) {
		part_id = * (int *) (data + stride * j); /* potentially, part_id != i */
		/* while this part_id same as previous part_id */
		do {
			(*eps)[j] = *(double *) (data + offset + stride * j);
			j++;
		} while (j < *nflaws && 
				*(int *) (data + stride * j) == 
				*(int *) (data + stride * (j - 1))); 
		
		(*flaws_tbl_lookup)[part_id * 2] = index;
		(*flaws_tbl_lookup)[part_id * 2 + 1] = j - index;
		index = j;
	}

	free (data);
	free (addrs);
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

int has_strength(SPHbody p) {
	if (p.data.strengthbody.is_strength == 0) /* redundant??? */
		return 0;
	if (p.u >= params.umelt)
		return 0;
	if (p.data.strengthbody.crack_len >= p.h * 4.)
		return 0;
	return 1;
}

void strength_force (double *grpmj, double *rhoij,
     double *sxxi,double *syyi,double *sxyi,double *sxzi,double *syzi,double *sxxj,double *syyj,double *
     sxyj,double *sxzj,double *syzj,double *dmi,double *dmj,double *dx,double *dy,double *dz,double *dfxi,double *dfyi,double *dfzi) {

	double redi, redj;
	double sigxxij, sigyyij, sigzzij, sigxyij, sigxzij, sigyzij;
	double tx, ty, tz;
//  material strength if pairs belong to same object

      redi=1.-*dmi * *dmi * *dmi;
      redj=1.-* dmj * *dmj * *dmj;
      sigxxij=(redi* *sxxi + redj* *sxxj)/ *rhoij;
      sigyyij=(redi* *syyi + redj* *syyj)/ *rhoij;
      sigzzij=(redi* *sxxi - redi* *syyi + redj* *sxxj - redj* *syyj)/ *rhoij;
      sigxyij=( *sxyi + redj* *sxyj)/ *rhoij;
      sigxzij=( *sxzi + redj* *sxzj)/ *rhoij;
      sigyzij=( *syzi + redj* *syzj)/ *rhoij;
      tx=sigxxij* *dx + sigxyij* *dy + sigxzij* *dz;
      ty=sigxyij* *dx + sigyyij* *dy + sigyzij* *dz;
      tz=sigxzij* *dx + sigyzij* *dy + sigzzij* *dz;

      *dfxi= *grpmj*tx;
      *dfyi= *grpmj*ty;
      *dfzi= *grpmj*tz;
}
