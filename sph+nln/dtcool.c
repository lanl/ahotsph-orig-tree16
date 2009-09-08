#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "cool.h"
#include "nrutil.h"
#include ""
#include ""
#include ""

void subcycle_lcool(double frac,float dt, float *lcool, float temp, float u,float udot,float rho)
{

/*checked in calling routine whether subcycling is necessary*/
/*divide timestep in half and calc udot from each. see if new 
udot at each new timestep meets some criterion. calc temp at 
each new timestep. if not, repeat subdividing until met. */
/* pass back or update p->udot */

    float dt0,temp0,u0,lcool0,u0; /*save initial values*/
    float dt1,temp1,u1,lcool1,udot1,dt1_tot;
    float kboltz,mh,l,m,t;
    int smallenough=0,done=0;

    kboltz=1.38065e-16;
    mh=1.67372267e-24;
    l=3.085678e18;
    m=1.988900e33;
    t=3.156300e7;

    u0=u;
    dt1=dt;
    dt1_tot=0.0;

    dt1=dt1/2.0;
    while(dt1_tot < dt) {
         if (abs(dt1*udot) < u*frac) { /*yes; update udot, calc new temperature */
              /*for this sub-timestep*/
              dt1_tot += dt1;	/*keep track of time*/
              u -= udot*dt1;	/*new u*/
              temp = 2.0*mh *u /(2.5 *kboltz);/*temp after dt1 timestep*/
              n = rho * m /(l*l*l) /(2.*mh);
              /*for next sub-timestep*/
              lcool = calc_lcool1(temp,1);	/*lcool at this temp (after dt1)*/
              udot -= lcool*n /(2.0*mh) *t*t*t/(l*l);	/*udot after this timestep*/
              udot_tot += udot;
	  } 
	  else { /*no; subdivide further*/
              dt1=dt1/2.0;
	  }
    }
    /*output udot_tot and (u0-u)/dt for comparison*/
}
