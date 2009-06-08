/*purpose: take in any necessary values from calling routine, 
and calculate the cooling, and return the cooling value to 
the calling routine-note
-CE*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
//#include "physics_sph.h"
//#include "vop.h"
//#include "fastflpt.h"
//#include "timers.h"
//#include "error.h"
#include "cool.h"
#include "nrutil.h" //ok to have this in

#ifndef M_1_PI
#define	M_1_PI 0.31830988618379067154
#endif
#define NKERNEL_TABLE 80000
#define MAX_INDEX (NKERNEL_TABLE+2)

#define NO_UPDATE 2

/******************************** NOTE: *******************************
 * take all command line outputs like error or status messages out for 
 * now. eventually should have a flag in the .ctl files to turn 
 * debugging on (i.e. with output messages) or off
 **********************************************************************/
//arrays in C: array[nrows][ncolumns]

double calc_lcool1(double temp,int extrapolate)
{
    extern float *tablep;	/*global array with cooling values*/
    extern float *ionfracp;	/*global array with ion fractions*/
    static double lcool;	/*holds final cooling term;returned*/
    double *lcoolp;	/*pointer to lcool*/
    double dy,df;	/*measure of accuracy returned from interp.*/
    double *dyp;	/*pointer to dy*/
    double *dfp;	/*pointer to df*/
    double ioncool=0.0;		/*hold intermediate cooling term*/
    double fracn;	/*holds ion fraction returned by interp. */
    double *fracnp;	/*pointer to frac*/
    double rowarr[2],interp[2],fracns[2],temps[2];	/*temporary arrays for interp.*/
    double logtemp; 	/*log(temp) for ionfracp interpolation*/
    int Nel=ionfracp[1];	/*number of elements in table*/
    int NMax=ionfracp[0];	/*number of grid points per ion in table*/
    long j;	/*holds index returned by locate routine*/
    long *jp;	/*pointer to j*/
    long j_prev;
    int n,m,offset,index; 	/*some indices for looping and arrays*/
	
    //initialize things so I don't get stupid warnings all the time:
    j=-999;
    dy=-999;
    jp=&j;
    dyp=&dy;
    lcoolp=&lcool;
    dfp=&df;
    fracnp=&fracn;
    offset=NMax+Nel;//offset=51 temperatures + 30 elements
	

    //locate the indices of the table
    //same for both tables as they go over the same range/grid points
    locate(&tablep[30], NMax, temp, jp); 
    j_prev=j;
	
    //if we're outside the table, do analytic cooling if extrapolate=0
    //or extrapolate if extrapolate=1:
    //extrapolation is still a little funky - CE
    //above table= j=-2, below table= j=-99
    if ((j==-2) || (j==-99))
    {
	if (extrapolate)
	{
	    //printf("extrapolating......\n");
	    //reset j for extrapolation
	    if (j==-2) j=0;
	    //now do interpolation, should automatically extrapolate as set up below
	    //reset j for extrapolation
	    if (j==-99) j=NMax-1;
	    //now do interpolation, should automatically extrapolate as set up below
	}
	else
	    return analytic_cool(temp); //we're done here
    } 
    //technically, need to interpolate over every ion in every element
    //then weigh by equil fraction at current T, and mix according to 
    //composition
    logtemp=log10(temp);
    //for (n=0;n<Nel;n++)	//loop over elements
    for (n=7;n<8;n++)	//loop over oxygen (no abundance tracking yet)
    {
        for (m=0;m<=(n+1);m++)	//loop over ions for each element,dont skip any
	{
	    index=(((n+1)*(n+2)>>1)-1+m)*51;	/*find row of element n, ion m*/
	    //re-assign table values so interpolation can be done in 
	    //double precision (table is floats).
	    rowarr[0]=tablep[j+Nel];
	    rowarr[1]=tablep[j+Nel+1];
	    interp[0]=tablep[index+offset+j];
	    interp[1]=tablep[index+offset+j+1];
	    //interpolate
	    polint(rowarr,interp,2,temp,lcoolp,dyp);
//THIS IS JUST A QUICK DIRTY FIX, MAKE MORE ROBUST!!!!
	    //reset value if extrapolated to negative value
	    if (lcool<0.0)
	    {
		if ((interp[0]-interp[1]) <0) lcool=0.0;
		//else fracn=1.0;//this case shouldn't happen
	    }
		
	    //find equil fraction of this ion:
	    temps[0]=ionfracp[j+2];
	    temps[1]=ionfracp[j+2+1];
	    //break-down of ionfracp index:
	    //-first two elements contain # grid points (NMax) and # of elements
	    //-then NMax temperature grid points, then the data.
	    //The data includes neutrals, skip these! (see read_ioncool-b.c) done
	    //with NMax*(n+1). change to *n to skip bare ions. 
	    fracns[0]=ionfracp[2+NMax+index+j];
	    fracns[1]=ionfracp[2+NMax+index+j+1];
	    //interpolate
	    polint(temps,fracns,2,logtemp,fracnp,dfp);
//THIS IS JUST A QUICK DIRTY FIX, MAKE MORE ROBUST!!!!
	    //reset value if extrapolated to unphysical value
	    if (fracn<0.0)
	    {
	        if ((fracns[0]-fracns[1]) <0) fracn=1.0e-30;
		else fracn=1.0;
	    }
	    //sum over all ions of given species. 
	    //also need to find equil fraction of ions (interpolate also)
	    //ioncool+=lcool*fracn*X_el[n];
	    ioncool+=lcool*fracn; //assume pure O composition
	}
    }
    return ioncool;
} /*end calc_lcool*/



double analytic_cool(double temp)
{
    /* From Chris's email; fit to Dalgarno and McCray (ARA&A
	 1972, 10, 375) and Sutherland and Dopita (ApJS, 88,
	 253) 
	 */
	double lcool;
	
	if (temp < 1.0e4)
		lcool = 1.0e-26 * exp(-1.0e5/temp) * sqrt(temp);//guessed; for O - CE
		//lcool = 1.0e-27 * exp(-1.0e2/temp) * sqrt(temp);
	else if (temp < 3.0e5)
		lcool=1.0e-21;
	else
		lcool=1.0e-21/(3.0*(log10(temp)-5.5)+1.0);
		//lcool=1.0e-21/(3.0*(log10(temp)-5.5)+1.0);
	
	return lcool;
}//end analytic_cool



/* this bisection routine is from NR for C: */
/*"Given an array xx[1..n], and given an value x, returns a value
j such that x is between xx[j] and xx[j+1]. xx must be monotonic,
either increasing or decreasing. j=-2 or j=-99 is returned to indicate
that x is out of range.*/
/*includes: none*/
/*call syntax: locate(&xx, N, x, j) */
/*~~~~~~~~~~~~~~~~WORKS!!! 04/23/2009~~~~~~~~~~~~~*/
void 
locate(float xx[], long Nel, float x, long *j)
{
   //floats should be enough for finding correct indices
   long ju, jm, jl;
   int ascnd,sign;

   jl=-1;
   ju=Nel;
   //check whether the table is in increasing (ascnd=1) or 
   //decreasing (ascnd=0) order
   ascnd=(xx[Nel-1] >= xx[0]);   //what does this line do - check whether
                             //xx[N]>xx[1] or not -CE
   if (ascnd) sign=1;
   if (!ascnd) sign=-1;

   while (ju-jl > 1)
   {
      jm=(ju+jl) >> 1;
      if ((x >= xx[jm]) == ascnd)   
         jl=jm;
      else
         ju=jm;
   }
   //if (true && true) && true (1 is true) then below (j<0) table
   if ( (((x - xx[0])*sign <0) && ((x-xx[Nel-1])*sign <0)) && 1)
	   *j=-2; 
   //if (not false && not false) && true (1 is true) then above (j>Nel) table
   else if ( (!((x - xx[0])*sign <0) && !((x-xx[Nel-1])*sign <0)) && 1)
	   *j=-99; 
   else *j=jl;
} //end locate



/*this is a polynomial interpolation routine from NR for C (S3.1):
"Given arrays xa[1..n] and and ya[1..n], and given a value x, this 
routine returns a value y, and an error estimate dy. If P(xP) is the 
polynomial of degree N-1 such that P(xa_i)=ya_i, i=1,...,n, then 
the returned value y=P(x).*/
/*includes: math.h, "nrutil.h"*/
/*call syntax: polint(&xx[14],&yy[14],4,x,yp,dyp) for 4-point 
interpolation between tabulated points [14..17]*/
/*~~~~~~~~~~~~~~~~~ WORKS!!!! 04/23/2009 ~~~~~~~~~~~~~~~~~*/
void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   //ZERO!!! offset is assumed in all indices
   int i,m,ns=0,size;
   double den,dif,dift,ho,hp,w;
   double *c, *d;

   //printf("interpolating between %E and %E .... get ", ya[0], ya[1]);
   //printf("xa= {%.1E, %.1E} \n", xa[0], xa[1]);
   //printf("ya= {%.1E, %.1E} \n", ya[0], ya[1]);

   //if we're exactly at one grid point, just return that value*/
   if ((x-xa[0]) == 0.0) 
	   *y=ya[0];
   else if ((x-xa[1]) == 0.0) 
	   *y=ya[1];
   else	//else interpolate*/
   {
   dif=fabs(x-xa[0]);
   c=dvector(0,n-1);
   d=dvector(0,n-1);
   //find the index ns of the closest table entry
   for (i=0;i<n;i++)
   {
      if ((dift=fabs(x-xa[i])) < dif)
      {
         ns=i;
         dif=dift;
      }
      //initialize the tableau of c's and d's
      c[i]=ya[i];
      d[i]=ya[i];
   }
   //initial approximation to y
   *y=ya[ns--];
   //for each column of the tableau ...
   for (m=0;m<n-1;m++)
   {
      //... loop over current c's and d's and update them
      for (i=0;i<n-m-1;i++)
      {
         ho=xa[i]-x;
         hp=xa[i+m+1]-x;
         w=c[i+1]-d[i];
	 den=ho-hp;
         //error message: two input xa's are identical to within roundoff
         //if (den==0.0) nerror("Error in routine polint");
         den=w/den;
         //update c's and d's
         d[i]=hp*den;
         c[i]=ho*den;
      }
      /*After each column in the tableau is completed, we decide wich 
        correction, c or d, we want to add to our accumulating value of
        y, i.e., which path to take through the tableau - forking up or 
        down. We do this in such a way as to take the most "straight 
        line" route through the tableau to its apex, updating ns 
        accordingly to keep track of where we are. This route keeps the 
        partial approximations centered (insofar as possible) on the 
        target x. the last dy added is thus the error indication. */
      *y += (*dy=(2*(ns+1)<(n-(m+1)) ? c[ns+1] : d[ns--]));
      //printf("%E\n",*y);
   }
   //printf("%E\n", *y);
   free_dvector(d,0,n-1);
   free_dvector(c,0,n-1);
   }
} //end polint



/*bilinear polynomial interpolation in 2 dimensions from NR for C: */
/*Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a 
submatrix of function values ya[1..m][1..n], tabulated at the grid points
defined by x1a and x2a; and given values x1 and x2 of the independent 
variables; this routine teturns as interpolated function value y, and an 
accuracy idication dy (based only on the interpolation in the x1 direction)*/
/*includes: "nrutil.h" */
/* call syntax: polin2d(&x1a[jj],&x2a[kk],&yap[0][0],2,2,x1,x2,yp,dyp) 
 * where x1a,x1, refers to rows and jj is the row number from which to 
 * start the interpolation, x2a,x2 refers to columns and kk is the column
 * number from which to start the interpolation, 
 * and yap is an array of pointers (which point to the values
 * to interpolate) like so:
 * yap[0][0]=&table[jj][kk]
 * yap[1][0]=&table[jj][kk+1]
 * yap[0][1]=&table[jj+1][kk]
 * yap[1][1]=&table[jj+1][kk+1]
 * (really, whatever the order, polin2d treats the first index of yap as 
 * rows and interpolates over those first)
 */
/*~~~~~~~~~~~~~~~~ WORKS!!!! 04/23/2009 ~~~~~~~~~~~~~~*/
void
polin2d(double x1a[], double x2a[], double **ya, int m, int n, double x1, 
        double x2, double *y, double *dy)
{
   //ZERO!!! offset is assumed in all indices
   void polint(double xa[], double ya[], int n, double x, double *y, double *dy);

   int j;
   double *ymtmp;
 
   ymtmp=dvector(0,m-1);
   //loop over rows
   for (j=0;j<m;j++)
   {
      //interpolate over the 'rows' of the grid square containing the point
      //(x1,x2), i.e. over (x1a[j],x2). Put answer into temporary storage
      polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
   }
   //do the final interpolation, i.e. sort of interpolate over the 
   //interpolated 'rows'
   polint(x1a,ymtmp,m,x1,y,dy);
   free_dvector(ymtmp,0,m-1);
} //end polin2d

