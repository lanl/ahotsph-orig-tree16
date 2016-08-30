#include <math.h>
#include "randoms.h"
#include "physics.h"
#include "physics_sph.h"
#include "bigmalloc.h"
#include "vop.h"
#include "singlio.h"
#include "fastflpt.h"
#include "mpmy.h"
#include "gc.h"
#include "Msgs.h"
#include "SDF.h"
#include "SDFread.h"
#include "SDFreadf.h"
#include "nrutil.h"
#include "cool.h"

	void *
InitRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
		SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, 
		int set_id, int setpvel, float new_h, float new_u)
{
	int i;
	SDF *sdfp;
	int massconf, xconf, yconf, zconf;
	int vxconf, vyconf, vzconf;
	int identconf;
	body *btab, *p; 
	int nobj, gnobj;
	SPHbody *SPHbtab, *q; 
	float hubble;

	singlPrintf("Reading \"%s\"\n", name);
	sdfp = SDFread(csdfp, (void **)btabp, gnobjp, nobjp, sizeof(body),
			"mass", offsetof(body, mass), &massconf,
			"x", offsetof(body, pos[0]), &xconf,
			"y", offsetof(body, pos[1]), &yconf,
			"z", offsetof(body, pos[2]), &zconf,
			"vx", offsetof(body, vel[0]), &vxconf,
			"vy", offsetof(body, vel[1]), &vyconf,
			"vz", offsetof(body, vel[2]), &vzconf,
			"ident", offsetof(body, ident), &identconf,
			NULL);
	nobj = *nobjp;
	gnobj = *gnobjp;
	btab = *(body **)btabp;
	Msgf(("Data read, nobj=%d, gnobj=%d\n", *nobjp, *gnobjp));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
				MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
		SinglError("Could not find %s %s %s %s in data file!\n",
				(massconf==0)? "mass" : "",
				(xconf==0)? "x" : "",
				(yconf==0)? "y" : "",
				(zconf==0)? "z" : "");
	}
	if (vxconf != vyconf || vxconf != vzconf){
		if (setpvel) SinglError("Missing velocity components!\n");
	}
	if (identconf == 0 || set_id){
		SinglWarning("No \"ident\" in file, numbering sequentially\n");
		FixId(btab, nobj, gnobj);
	}

	SPHbtab = Malloc(nobj * sizeof(SPHbody));
	singlPrintf("Setting h to %f\n", new_h);
	singlPrintf("Setting u to %f\n", new_u);
	SDFgetfloatOrDie(sdfp, "hubble",  &hubble);
	for (i = 0; i < nobj; i++) {
		p = btab+i;
		q = SPHbtab+i;
		q->mass = p->mass * 0.1; /* 10 percent baryons */
		p->mass *= 0.9;
		VV(q->pos, = p->pos);
		/* Offset a little so tree build doesn't fail due to identical pos */
		q->pos[0] += 3.9;
		VV(q->vel, = hubble * p->pos);
		q->ident = p->ident + gnobj;
		q->h = new_h;
		q->u = new_u;
	}
	*SPHgnobjp = gnobj;
	*SPHnobjp = nobj;
	*SPHbtabp = SPHbtab;
	return sdfp;
}

	void
DarkSPHTestData(void *csdfp, void **btabp, int *gnobjp, int *nobjp, 
		SPHbody **SPHbtabp, int *SPHgnobjp, int *SPHnobjp, int periodic)
{
	int i;
	ran_state ranstate;
	int seed, cencon;
	int start;
	int gnobj, nobj;
	body *btab, *p;
	int SPHgnobj, SPHnobj;
	SPHbody *SPHbtab, *q;
	float new_u;
	float h, rsq;

	singlPrintf("Generating random dataset\n");
	if (SDFgetint(csdfp, "nobj", &gnobj))
		SinglError("Sorry, you've got to have an \"nobj\"\n");
	SDFgetintOrDefault(csdfp, "seed", &seed, 123);
	SDFgetintOrDefault(csdfp, "cencon", &cencon, 0);
	SDFgetfloatOrDefault(csdfp, "new_u", &new_u, 0.0);
	singlPrintf("int seed = %d;\n", seed);
	singlPrintf("int cencon = %d;\n", cencon);

	NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start);
	btab = (body *) Malloc(nobj*sizeof(SPHbody));
	SPHgnobj = gnobj;
	SPHnobj = nobj;
	SPHbtab = (SPHbody *) Malloc(SPHnobj*sizeof(SPHbody));
	ran_init(seed*(MPMY_Procnum()+1), &ranstate);
	for (p = &btab[0]; p < &btab[nobj]; p++) {
#ifdef __PARAGON__
		clear_tregs();	/* avoid system bug */
#endif
		p->mass = 1.0 / gnobj;		 /*   set masses equal */
		if (periodic)
			rsq = cube_rand(&ranstate, NDIM, p->pos);
		else
			rsq = sphere_rand(&ranstate, NDIM, p->pos);
		VS(p->vel, = 0.0);
	}
	h = pow((float)8.5/SPHgnobj, .333333);
	for (i = 0; i < nobj; i++) {
		p = btab+i;
		q = SPHbtab+i;
		q->mass = p->mass * 0.1; /* 10 percent baryons */
		p->mass = p->mass * 0.9;
		/* Offset a little so tree build doesn't fail due to identical pos */
		VVS(q->pos, = p->pos, + .001);
		VS(q->vel, = (float)0.0);
		q->ident = p->ident + gnobj;
		q->h = h;
		q->u = new_u;
	}
	singlPrintf("Extracted 10%% baryons from dark matter input\n");
	FixId(btab, nobj, gnobj);
	FixNterms(btab, nobj);
	SPHFixId(SPHbtab, SPHnobj, SPHgnobj);
	SPHFixNterms(SPHbtab, SPHnobj);
	*gnobjp = gnobj;
	*nobjp = nobj;
	*btabp = btab;
	*SPHgnobjp = SPHgnobj;
	*SPHnobjp = SPHnobj;
	*SPHbtabp = SPHbtab;
}


	void *
SPHRead(char *name, void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp,
		int set_id, int setpvel, float new_h, float new_u)
{
	SDF *sdfp;
	int massconf, xconf, yconf, zconf;
	int vxconf, vyconf, vzconf;
	int hconf, uconf;
	int identconf, windidconf;
	SPHbody *btab, *p; 
	int nobj, gnobj;

	singlPrintf("SPHReading \"%s\"\n", name);
	sdfp = SDFreadf(name, (void **)btabp, gnobjp, nobjp, sizeof(SPHbody),
			"mass", offsetof(SPHbody, mass), &massconf,
			"x", offsetof(SPHbody, pos[0]), &xconf,
			"y", offsetof(SPHbody, pos[1]), &yconf,
			"z", offsetof(SPHbody, pos[2]), &zconf,
			"vx", offsetof(SPHbody, vel[0]), &vxconf,
			"vy", offsetof(SPHbody, vel[1]), &vyconf,
			"vz", offsetof(SPHbody, vel[2]), &vzconf,
			"u", offsetof(SPHbody, u), &uconf,
			"h", offsetof(SPHbody, h), &hconf,
			"ident", offsetof(SPHbody, ident), &identconf,
			"windid", offsetof(SPHbody, windid), &windidconf,
			NULL);
	nobj = *nobjp;
	gnobj = *gnobjp;
	btab = *btabp;
	Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", *nobjp, *gnobjp));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
				MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
		SinglError("Could not find %s %s %s %s in data file!\n",
				(massconf==0)? "mass" : "",
				(xconf==0)? "x" : "",
				(yconf==0)? "y" : "",
				(zconf==0)? "z" : "");
	}
	if (vxconf != vyconf || vxconf != vzconf){
		if (setpvel) SinglError("Missing velocity components!\n");
	}
	if (identconf == 0 || set_id){
		SinglWarning("No \"ident\" in file, numbering sequentially\n");
		SPHFixId(btab, nobj, gnobj);
	}
	if (windidconf == 0) {
		SinglWarning("No \"windid\" in file; are you using wind source?\n");
	}
	if (new_h != (float)0.0) {
		singlPrintf("Setting h to %f\n", new_h);
		for (p = btab; p < btab+nobj; p++) p->h = new_h;
	} else if (hconf == 0) {
		SinglError("No h in data file\n");
	}
	if (new_u != (float)0.0) {
		singlPrintf("Setting u to %f\n", new_u);
		for (p = btab; p < btab+nobj; p++)  p->u = new_u;
	} else if (uconf == 0) {
		SinglError("No u in data file\n");
	}
	return sdfp;
}

/*read in file with abundance data */
	void *
SPHRead_nw(char *name, void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp,
		int set_id, int setpvel, float new_h, float new_u)
{
	SDF *sdfp;
	int massconf, xconf, yconf, zconf;
	int vxconf, vyconf, vzconf;
	int hconf, uconf;
	int identconf, tempconf;
	int f1conf, p1conf, m1conf, Y_elconf;
	SPHbody *btab, *p; 
	int nobj, gnobj;

	singlPrintf("SPHReading \"%s\"\n", name);
	singlPrintf("if you get Y_el related error, check in file sphinit.c\n");
	sdfp = SDFreadf(name, (void **)btabp, gnobjp, nobjp, sizeof(SPHbody),
			"mass", offsetof(SPHbody, mass), &massconf,
			"x", offsetof(SPHbody, pos[0]), &xconf,
			"y", offsetof(SPHbody, pos[1]), &yconf,
			"z", offsetof(SPHbody, pos[2]), &zconf,
			"vx", offsetof(SPHbody, vel[0]), &vxconf,
			"vy", offsetof(SPHbody, vel[1]), &vyconf,
			"vz", offsetof(SPHbody, vel[2]), &vzconf,
			"u", offsetof(SPHbody, u), &uconf,
			"h", offsetof(SPHbody, h), &hconf,
			"ident", offsetof(SPHbody, ident), &identconf,
			"temp", offsetof(SPHbody, temp), &tempconf,
			/* comment this line out if no Y_el in sdf file */
			/*"Y_el", offsetof(SPHbody, Y_el), &Y_elconf;*/
			/*"f1", offsetof(SPHbody, composition[0].abund), &f1conf,*/
			/*"p1", offsetof(SPHbody, composition[0].np), &p1conf,*/
			/*"m1", offsetof(SPHbody, composition[0].nn), &m1conf,*/
			/* or .... */
			/* let's do this the dummest, ugliest, most painful way possible. he. he. he. -CIE */
			"f1", offsetof(SPHbody, data.nucnetw.abund[0]), &f1conf,
			"f2", offsetof(SPHbody, data.nucnetw.abund[1]), &f1conf,
			"f3", offsetof(SPHbody, data.nucnetw.abund[2]), &f1conf,
			"f4", offsetof(SPHbody, data.nucnetw.abund[3]), &f1conf,
			"f5", offsetof(SPHbody, data.nucnetw.abund[4]), &f1conf,
			"f6", offsetof(SPHbody, data.nucnetw.abund[5]), &f1conf,
			"f7", offsetof(SPHbody, data.nucnetw.abund[6]), &f1conf,
			"f8", offsetof(SPHbody, data.nucnetw.abund[7]), &f1conf,
			"f9", offsetof(SPHbody, data.nucnetw.abund[8]), &f1conf,
			"f10", offsetof(SPHbody, data.nucnetw.abund[9]), &f1conf,
			"f11", offsetof(SPHbody, data.nucnetw.abund[10]), &f1conf,
			"f12", offsetof(SPHbody, data.nucnetw.abund[11]), &f1conf,
			"f13", offsetof(SPHbody, data.nucnetw.abund[12]), &f1conf,
			"f14", offsetof(SPHbody, data.nucnetw.abund[13]), &f1conf,
			"f15", offsetof(SPHbody, data.nucnetw.abund[14]), &f1conf,
			"f16", offsetof(SPHbody, data.nucnetw.abund[15]), &f1conf,
			"f17", offsetof(SPHbody, data.nucnetw.abund[16]), &f1conf,
			"f18", offsetof(SPHbody, data.nucnetw.abund[17]), &f1conf,
			"f19", offsetof(SPHbody, data.nucnetw.abund[18]), &f1conf,
			"f20", offsetof(SPHbody, data.nucnetw.abund[19]), &f1conf,
			NULL);
	nobj = *nobjp;
	gnobj = *gnobjp;
	btab = *btabp;
	Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", *nobjp, *gnobjp));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
				MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
		SinglError("Could not find %s %s %s %s in data file!\n",
				(massconf==0)? "mass" : "",
				(xconf==0)? "x" : "",
				(yconf==0)? "y" : "",
				(zconf==0)? "z" : "");
	}
	if (vxconf != vyconf || vxconf != vzconf){
		if (setpvel) SinglError("Missing velocity components!\n");
	}
	if (tempconf == 0) {
		SinglWarning("Temperature not found in input data file.\n");
	}
	if (identconf == 0 || set_id){
		SinglWarning("No \"ident\" in file, numbering sequentially\n");
		SPHFixId(btab, nobj, gnobj);
	}
	if (new_h != (float)0.0) {
		singlPrintf("Setting h to %f\n", new_h);
		for (p = btab; p < btab+nobj; p++) p->h = new_h;
	} else if (hconf == 0) {
		SinglError("No h in data file\n");
	}
	if (new_u != (float)0.0) {
		singlPrintf("Setting u to %f\n", new_u);
		for (p = btab; p < btab+nobj; p++)  p->u = new_u;
	} else if (uconf == 0) {
		SinglError("No u in data file\n");
	}
	return sdfp;
}

/*read in file with strength data */
	void *
SPHReadStrength(char *name, void *csdfp, Strengthbody **btabp, int *gnobjp, int *nobjp,
		int set_id, int setpvel, float new_h, float new_u)
{
	SDF *sdfp;
	int massconf, xconf, yconf, zconf;
	int vxconf, vyconf, vzconf;
	int hconf, uconf;
	int identconf;
	int dmgconf, ddmgdtconf, stressconf[9], dstressdtconf[9];
	Strengthbody *btab, *p; 
	int nobj, gnobj;
	int i;

	singlPrintf("SPHReading \"%s\"\n", name);
	sdfp = SDFreadf(name, (void **)btabp, gnobjp, nobjp, sizeof(Strengthbody),
			"mass", offsetof(Strengthbody, mass), &massconf,
			"x", offsetof(Strengthbody, pos[0]), &xconf,
			"y", offsetof(Strengthbody, pos[1]), &yconf,
			"z", offsetof(Strengthbody, pos[2]), &zconf,
			"vx", offsetof(Strengthbody, vel[0]), &vxconf,
			"vy", offsetof(Strengthbody, vel[1]), &vyconf,
			"vz", offsetof(Strengthbody, vel[2]), &vzconf,
			"u", offsetof(Strengthbody, u), &uconf,
			"h", offsetof(Strengthbody, h), &hconf,
			"ident", offsetof(Strengthbody, ident), &identconf,
			"dmg", offsetof(Strengthbody, dmg), &dmgconf,
			"ddmgdt", offsetof(Strengthbody, ddmgdt), &ddmgdtconf,
			"stressxx", offsetof(Strengthbody, stress[0]), &stressconf[0],
			"stressxy", offsetof(Strengthbody, stress[1]), &stressconf[1],
			"stressxz", offsetof(Strengthbody, stress[2]), &stressconf[2],
			"stressyx", offsetof(Strengthbody, stress[3]), &stressconf[3],
			"stressyy", offsetof(Strengthbody, stress[4]), &stressconf[4],
			"stressyz", offsetof(Strengthbody, stress[5]), &stressconf[5],
			"stresszx", offsetof(Strengthbody, stress[6]), &stressconf[6],
			"stresszy", offsetof(Strengthbody, stress[7]), &stressconf[7],
			"stresszz", offsetof(Strengthbody, stress[8]), &stressconf[8],
			"dstressxxdt", offsetof(Strengthbody, dstressdt[0]), &dstressdtconf[0],
			"dstressxydt", offsetof(Strengthbody, dstressdt[1]), &dstressdtconf[1],
			"dstressxzdt", offsetof(Strengthbody, dstressdt[2]), &dstressdtconf[2],
			"dstressyxdt", offsetof(Strengthbody, dstressdt[3]), &dstressdtconf[3],
			"dstressyydt", offsetof(Strengthbody, dstressdt[4]), &dstressdtconf[4],
			"dstressyzdt", offsetof(Strengthbody, dstressdt[5]), &dstressdtconf[5],
			"dstresszxdt", offsetof(Strengthbody, dstressdt[6]), &dstressdtconf[6],
			"dstresszydt", offsetof(Strengthbody, dstressdt[7]), &dstressdtconf[7],
			"dstresszzdt", offsetof(Strengthbody, dstressdt[8]), &dstressdtconf[8],
			NULL);
	nobj = *nobjp;
	gnobj = *gnobjp;
	btab = *btabp;
	Msgf(("Data read, SPHnobj=%d, SPHgnobj=%d\n", *nobjp, *gnobjp));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
				MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
		SinglError("Could not find %s %s %s %s in data file!\n",
				(massconf==0)? "mass" : "",
				(xconf==0)? "x" : "",
				(yconf==0)? "y" : "",
				(zconf==0)? "z" : "");
	}
	if (vxconf != vyconf || vxconf != vzconf){
		if (setpvel) SinglError("Missing velocity components!\n");
	}
	if (identconf == 0 || set_id){
		SinglWarning("No \"ident\" in file, numbering sequentially\n");
		SPHFixId(btab, nobj, gnobj);
	}
	if (new_h != (float)0.0) {
		singlPrintf("Setting h to %f\n", new_h);
		for (p = btab; p < btab+nobj; p++) p->h = new_h;
	} else if (hconf == 0) {
		SinglError("No h in data file\n");
	}
	if (new_u != (float)0.0) {
		singlPrintf("Setting u to %f\n", new_u);
		for (p = btab; p < btab+nobj; p++)  p->u = new_u;
	} else if (uconf == 0) {
		SinglError("No u in data file\n");
	}
	for (i = 0; i < 9; i++) {
		if (stressconf[i] == 0)
			SinglWarning("Missing stress tensor component %d.\n", i);
		if (dstressdtconf[i] == 0)
			SinglWarning("Missing stress tensor rate of change component %d.\n", i);
	}
	return sdfp;
}

	void *
DarkRead(char *name, void *csdfp, void **btabp, int *gnobjp, int *nobjp,
		int set_id, int setpvel)
{
	SDF *sdfp;
	int massconf, xconf, yconf, zconf;
	int vxconf, vyconf, vzconf;
	int lxconf, lyconf, lzconf, accmassconf;
	int identconf;
	body *btab;
	int nobj, gnobj;

	singlPrintf("Reading \"%s\"\n", name);
	sdfp = SDFreadf(name, (void **)btabp, gnobjp, nobjp, sizeof(body),
			"mass", offsetof(body, mass), &massconf,
			"x", offsetof(body, pos[0]), &xconf,
			"y", offsetof(body, pos[1]), &yconf,
			"z", offsetof(body, pos[2]), &zconf,
			"vx", offsetof(body, vel[0]), &vxconf,
			"vy", offsetof(body, vel[1]), &vyconf,
			"vz", offsetof(body, vel[2]), &vzconf,
			"lx", offsetof(body, l[0]), &lxconf,
			"ly", offsetof(body, l[1]), &lyconf,
			"lz", offsetof(body, l[2]), &lzconf,
			"accmass", offsetof(body, accmass), &accmassconf,
			"ident", offsetof(body, ident), &identconf,
			NULL);
	nobj = *nobjp;
	gnobj = *gnobjp;
	btab = *btabp;
	Msgf(("Data read, nobj=%d, gnobj=%d\n", *nobjp, *gnobjp));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
				MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if (massconf==0 || xconf==0 || yconf==0 || zconf==0) {
		SinglError("Could not find %s %s %s %s in data file!\n",
				(massconf==0)? "mass" : "",
				(xconf==0)? "x" : "",
				(yconf==0)? "y" : "",
				(zconf==0)? "z" : "");
	}
	if (vxconf != vyconf || vxconf != vzconf){
		if (setpvel) SinglError("Missing velocity components!\n");
	}
	if (accmassconf==0 || lxconf==0 || lyconf==0 || lzconf==0) {
		SinglError("Could not find %s %s %s %s in data file!\n", 
				(lxconf==0)? "lx" : "",
				(lyconf==0)? "ly" : "",
				(lzconf==0)? "lz" : "",
				(accmassconf==0)? "accmass" : "");
	}
	if (identconf == 0 || set_id){
		SinglWarning("No \"ident\" in file, numbering sequentially\n");
		FixId(btab, nobj, gnobj);
	}
	return sdfp;
}

	void
SPHTestData(void *csdfp, SPHbody **btabp, int *gnobjp, int *nobjp, int periodic)
{
	ran_state ranstate;
	int seed, cencon;
	int start;
	int gnobj, nobj;
	SPHbody *btab, *p;
	float new_u;
	float h;
	/*	 float rsq; */

	singlPrintf("Generating random dataset\n");
	if (SDFgetint(csdfp, "nobj", gnobjp))
		SinglError("Sorry, you've got to have an \"nobj\"\n");
	gnobj = *gnobjp;
	SDFgetintOrDefault(csdfp, "seed", &seed, 123);
	SDFgetintOrDefault(csdfp, "cencon", &cencon, 0);
	SDFgetfloatOrDefault(csdfp, "new_u", &new_u, 0.0);
	singlPrintf("int seed = %d;\n", seed);
	singlPrintf("int cencon = %d;\n", cencon);

	NobjInitial(gnobj, MPMY_Nproc(), MPMY_Procnum(), &nobj, &start);
	btab = (SPHbody *) Malloc(nobj*sizeof(SPHbody));
	ran_init(seed*(MPMY_Procnum()+1), &ranstate);
	h = pow((float)8.5/gnobj, .333333);
	for (p = &btab[0]; p < &btab[nobj]; p++) {
#ifdef __PARAGON__
		clear_tregs();	/* avoid system bug */
#endif
		p->mass = 1.0 / gnobj;		 /*   set masses equal */
		/* Removed because cube_rand and sphere_rand both expect positions
		   as floats anyway... */
		/* 	if (periodic) { */
		/* 		rsq = cube_rand(&ranstate, NDIM, p->pos); */
		/* 	} else { */
		/* 		rsq = sphere_rand(&ranstate, NDIM, p->pos); */
		/* 	} */
		/* 	if (cencon == 1) { */
		/* 		rsq = -1.0/sqrt(rsq); */
		/* 		VV(p->vel, = rsq * p->pos); */
		/* 	} else { */
		/* 		VS(p->vel, = 0.0); */
		/* 	} */
		p->h = h;
		p->u = new_u;
	}
	SPHFixId(btab, nobj, gnobj);
	SPHFixNterms(btab, nobj);
	*nobjp = nobj;
	*btabp = btab;
}

	void *
WindRead(char *name, void *csdfp, windbody **btabp, int *gnobjp, int *nobjp)
{
	SDF *sdfp;
	int xconf, yconf, zconf;
	int vxconf, vyconf, vzconf;
	int rhoconf, vwindconf, uwindconf, identconf;
	windbody *btab;
	int nobj, gnobj;

	singlPrintf("Reading \"%s\"\n", name);
	sdfp = SDFreadwind(name, (void **)btabp, gnobjp, nobjp, sizeof(windbody),
			"xwind", offsetof(windbody, pos[0]), &xconf,
			"ywind", offsetof(windbody, pos[1]), &yconf,
			"zwind", offsetof(windbody, pos[2]), &zconf,
			"vxwind", offsetof(windbody, vel[0]), &vxconf,
			"vywind", offsetof(windbody, vel[1]), &vyconf,
			"vzwind", offsetof(windbody, vel[2]), &vzconf,
			"rhowind", offsetof(windbody, rhowind), &rhoconf,
			"vwind", offsetof(windbody, vwind), &vwindconf,
			"uwind", offsetof(windbody, uwind), &uwindconf,
			"identwind", offsetof(windbody, ident), &identconf,
			NULL);
	nobj = *nobjp;
	gnobj = *gnobjp;
	btab = *btabp;
	Msgf(("Wind data read, windnobj=%d, windgnobj=%d\n", *nobjp, *gnobjp));
	Msgf(("Nproc:%d, Procnum: %d, Doc: %d\n",
				MPMY_Nproc(), MPMY_Procnum(), ilog2(MPMY_Nproc())));
	if (xconf==0 || yconf==0 || zconf==0) {
		SinglError("Could not find %s %s %s in wind file!\n",
				(xconf==0)? "x" : "",
				(yconf==0)? "y" : "",
				(zconf==0)? "z" : "");
	}
	if (vxconf != vyconf || vxconf != vzconf){
		SinglError("Missing velocity components!\n");
	}
	if (identconf == 0){
		SinglError("No \"ident\" in wind file, aborting\n");
	}
	return sdfp;
}
