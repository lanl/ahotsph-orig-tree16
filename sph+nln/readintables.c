#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<strings.h>
#include<assert.h>
#include<time.h>
#include<errno.h>
#include"error.h"
#include"fastflpt.h"
#include"Msgs.h"
#include"physics.h"
#include"physics_sph.h"
#include"stk.h"
#include"vop.h"
#include"singlio.h"
#include"mpmy.h"
#include"timers.h"
#include"cool.h"
#include"nrutil.h"

extern float *tablep;
extern float *ionfracp;
//float X_el[30];

/*
void init_X(void)
{
    extern float X_el[30];
    int n;
    double sum=1.0;

    for (n=0;n<30;n++) X_el[n]=0.0;
    X_el[0]=0.0;//0.715; //H mass fraction from Lodders 2003
    //X_el[1]=0.271.; //He fraction
    X_el[1]=0.0;//pow(10,-1.01); //He, abundance fraction
    //Z=0.014= sum(mass of all metals)/M(H)
    X_el[5]=0.0;//pow(10,-3.44); //C, abundance fraction
    X_el[6]=0.0;//pow(10,-3.95); //N, abundance fraction
    X_el[7]=1.0;//pow(10,-3.07); //O, abundance fraction
    X_el[9]=0.0;//pow(10,-3.91); //Ne, abundance fraction
    X_el[10]=0.0;//pow(10,-5.67); //Na, abundance fraction
    X_el[11]=0.0;//pow(10,-4.42); //Mg, abundance fraction
    X_el[12]=0.0;//pow(10,-5.53); //Al, abundance fraction
    X_el[13]=0.0;//pow(10,-4.45); //Si, abundance fraction
    X_el[15]=0.0;//pow(10,-4.79); //S, abundance fraction
    X_el[16]=0.0;//pow(10,-6.50); //S, abundance fraction
    X_el[17]=0.0;//pow(10,-5.44); //Ar, abundance fraction
    X_el[19]=0.0;//pow(10,-5.64); //Ca, abundance fraction
    X_el[25]=0.0;//pow(10,-4.33); //Fe, abundance fraction
    X_el[27]=0.0;//pow(10,-5.75); //Ni, abundance fraction
    //
    //for (n=1;n<30;n++) sum+=X_el[n];
    //sum=1.0/sum;
    //printf("nH/N=%E\n",sum);
    //X_el[0]=sum;
    //for (n=1;n<30;n++) X_el[n]*sum;
    //for (n=0;n<30;n++) printf("X_el[%d]=%E\n",n,X_el[n]);
    //
}
*/



void init_CoolTable(void)
{
    FILE *File1p;	//pointer to file with cooling curves
    FILE *file2p;	//pointer to file with ion fractions
    long lSize;		//holds file size (number of characters in file)
    extern float *tablep;	//global array for cooling curves
    extern float *ionfracp;	//global array for ion fractions
    int i,j,k;		//indices for looping through arrays
    int iMax;		//max number of grid points
    int jMax;		//max number of
    int ZMax;		//max number of elements in tables
    int index;		//index to access correct array element
    int myint;		//holds un-needed integers read in from files
    char mychar,myline[50];	//holds new-lines and text read in from files

    //open table with cooling curves
    File1p = fopen("/home/cellinge/SNSPH.dir/tree16/sph+nln/CHIANTI-COOLING.dat", "r");
    if (File1p == NULL) 
	Error("error opening cooling curves: %s\n",strerror(errno));

    //open table with ion fractions
    file2p=fopen("/home/cellinge/SNSPH.dir/tree16/sph+nln/mazzotta_etal_9.ioneq","r");
    if (File1p == NULL) 
	Error("error opening ion fractions: %s\n",strerror(errno));

    //obtain file size for cooling curves:
    fseek (File1p, 0, SEEK_END);
    lSize = ftell (File1p);
    rewind (File1p);

    //allocate memory to contain the whole file:
    //lSize contains the total number of characters in the file. the numbers
    //are in single precision, scientific notation, and are ~11 characters 
    //long on average. divide lSize by 12 to get more appropriate size for 
    //total array.
    tablep = (float*) malloc(sizeof(float) *lSize/12);
    //tablep is a pointer to the beginning of the allocated memory block.
    //individual elements can be accessed with tablep[nelem].
    //if (tablep == NULL) {fputs ("memory allocation for cooling curves failed",stderr); exit(2);}

    //obtain file size for ion fractions:
    fseek (file2p,0,SEEK_END);
    lSize=ftell(file2p);
    rewind(file2p);

    //allocate memory to contain the whole file:
    ionfracp= (float*) malloc(sizeof(float)*lSize/10);
    //if (ionfracp==NULL) {fputs("memory allocation for ionfrac failed",stderr); exit(2);}

    fgets(myline,50,File1p);//read in first line of text in cooling curves
    fgets(myline,50,File1p);//read in second line of text in cooling curves

    //get max number of rows/columns. for-loop won't work with tablep[0]
    //directly.
    fscanf(file2p, "%f %f", &ionfracp[0],&ionfracp[1]);
    iMax=(int)(ionfracp[0]);
    //max number(Z) of elements in table
    ZMax=(int)(ionfracp[1]);

    //read in number (Z) of element for cooling curves
    for (i=0;i<ZMax;i++) fscanf(File1p,"%E",&tablep[i]);

    mychar=fgetc(File1p);//read in extra new-line in cooling curves
    fgets(myline,50,File1p);//read in line "temperatures...." in cooling curves

    //read in temperatures for cooling curves; index=ZMax+i
    for (i=0;i<iMax;i++) fscanf(File1p,"%E",&tablep[i+ZMax]);

    //read in log(temperatures) for ion fractions
    for (i=0;i<iMax;i++) fscanf(file2p,"%E",&ionfracp[i+2]);
	
    fgets(myline,50,File1p);//get trailing new-line in cooling curves
    fgets(myline,50,File1p);//get trailing new-line in cooling curves

    //read in data
    for (j=0;j<ZMax;j++)//loop over elements
    { 
        fgets(myline,50,File1p);//get line with element name
	//jMax=tablep[j];
	for (k=0;k<=(j+1);k++)//loop over ions in current element
	{
	    fscanf(file2p,"%i",&myint);
	    fscanf(file2p,"%d",&myint);
            for (i=0;i<iMax;i++)//loop over data point in current ion in current element 
	    {
	        //to get correct index:sum from 0 to number (Z) of element, 
		//add number of ionization state (e.g. 0 for neutral, 1 for 
		//singly ionized, etc), add offset (=ZMax+iMax=30+51)
	        index=((((j+1)*(j+2))>>1)-1+k)*iMax+i;
    	        fscanf(File1p,"%E",&tablep[2*iMax+ZMax+index]);
	        fscanf(file2p,"%E",&ionfracp[iMax+2+index]);
	/* in CHIANTI-file, neutrals (bare ions?) seem to be missing, but are
	 * present in Mazzotta et al. ion fractions file. In the above set-up
	 * the neutrals from the ion-fractions file are read in, and the ions
	 * in tablep and ionfracp are lined up (i.e. each row of ionfrac
	 * and tablep contain the same ion of the same element), and the 
	 * neutral atoms are read in as zeros in tablep).
	 */
	    }
	}
	mychar=fgetc(File1p);//get first trailing new-line
	mychar=fgetc(File1p);//get second trailing new-line
    }

    //close file
    fclose(File1p);
    fclose(file2p);
    //DO NOT free(tablep); UNTIL THE END OF THE WHOLE PROGRAM!!!!
}

