extern double *flaw_actv_tbl;
extern double vol_scaling;
extern int *flaw_actv_tbl_lookup;
extern double G_shear;
extern double E_Young;
extern double K_bulk;
extern double YieldStr;
extern double u_melt;

void init_defects_table(int gnobj, int Nflaws, double **eps, int **flaws_tbl_lookup, float kVol, float m);
void read_defects_table(SDF *sdfp, int *nflaws, double **eps, int **flaws_tbl_lookup);
void write_defects_table (char *name, int gnobj, int nflaws, double *eps, int *flaws_tbl_lookup);
int has_strength(SPHbody p);
