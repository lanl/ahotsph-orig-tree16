extern double *flaw_actv_tbl;
extern int *flaw_actv_tbl_lookup;

void init_defects_table(int gnobj, int Nflaws, double **eps, int **flaws_tbl_lookup, float kVol, float m);
void read_defects_table(SDF *sdfp, int *nflaws, double **eps, int **flaws_tbl_lookup);
void write_defects_table (char *name, int gnobj, int nflaws, double *eps, int *flaws_tbl_lookup);
