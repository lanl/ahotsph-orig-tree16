# Do I want float or double output? 0=only double 1=only float 2=both 
int do_floatoutput=0;
int do_debugoutput=0;

# Define interpolation range for the startup grid
# if omitted the size of the simulation will be used
struct {double outrmin;}[3] = {-3.4e3 ,-3.4e3 ,-3.4e3 };
struct {double outrmax;}[3] = { 3.4e3 , 3.4e3 , 3.4e3 };
int n_x = 32;
int n_y = 32;
int n_z = 32;
int dogrid = 0;
#char Msg_turn_on[] = "decomp.c,tree.c,sph.c";

int targetnobj=200000; 
double outerbound= 3.4e3;
double innerbound= 1.0e-2;
int keepcenterfixed=0;

#Define the type of input data 
# 1: Radial inputmodel in ascii table, 
# 2: Cylindrical grid in binary format 
# 3: Radial Power law
# 4: Cartesian grid in binary format 
int inputoption=1;

# Cartesian grid data
int cart_dimx=557;                    # Number of radial grid points
double cart_minx=-0.76734252;        # First x grid point  
double cart_maxx=0.066657471;        # Last x grid point

int cart_dimy=478;                    # Number of z grid points
double cart_miny=-0.35769973;	      # First z grid point
double cart_maxy=0.35780027;          # Last z grid point

int cart_dimz=437;                    # Number of theta grid points
double cart_minz=-0.32639764;         # First theta grid point
double cart_maxz=0.32760235;          # Last theta grid point

double cart_xcenter=-0.35039370;      # Offset for grid
double cart_ycenter=0.;               # (SPH sphere will be centered
double cart_zcenter=0.;               # around these coordinates)

char cartfile_rho[]="donor_rho.cdat"; # Binary file holding grid rho values
char cartfile_h[]="donor_h.cdat";     # Binary file holding grid h values


# Cylindrical grid data
int dimr=127;                    # Number of radial grid points
double cyl_minr=0.0039525692;    # First radial grid point  
double cyl_maxr=1.;              # Last radial grid point

int dimz=49;                     # Number of z grid points
double cyl_minz=0.0039525692;    # First z grid point
double cyl_maxz=0.38339921;      # Last z grid point

int dimtheta=256;                # Number of theta grid points
double cyl_mintheta=0.;          # First theta grid point
double cyl_maxtheta=6.2586416;   # Last theta grid point

double cyl_xcenter=-0.4008;      # Offset for grid
double cyl_ycenter=0.;           # (SPH sphere will be centered
double cyl_zcenter=0.;           # around these coordinates)

# EOS
int do_eospolytrope=0;
double kpolytrope=0.0372;


# Number of loops WVT should go through
int nloop=40;

# Maximal number of loops the smoothing length determination should go through
# to keep the number of neighbors approximately at targetneighbors
double targetneighbors=40.;
int nhloop=10;

# Maximal number of loop the mass determination goes through. 
# WARNING: Make sure you know what you are doing here. In general, nmassloop=0
# should be sufficient. If you are resolution limited in the center, anything 
# else may give you the wrong answer!
int nmassloop=1;
double rhomin=1e-30;

# If set, u is calculated from pressure profile, in case we can't exactly 
# match rho(r). Default is 0.
int do_hydrostatic=0;

# Central Particle properties
int do_center=0;
int center_dual=0;
int center_sphfixed=0;
int center_posfixed=0;
double center_h=0.1;
double center_grav_mass=0.391973;

# New Spline Kernel Coefficients from Chris
# if omitted, the default is the Monaghan kernel 
#int kernel_ncoef1 = 5;
#struct {double kernel_coef1;}[5] = {1.2798, 0.0, -3.25703558174757,
#3.0425, -0.83492};
#int kernel_ncoef2 = 5;
#struct {double kernel_coef2;}[5] = {2.1164, -3.3596555, 1.797,
#-0.32368, 0.000272};

# Output directory
char outdir[] = "./";

int do_externalstart=0;
char startfile[]="oldwvtfinal.sdf";


# SDF-EOF
