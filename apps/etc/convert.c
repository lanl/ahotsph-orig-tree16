/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDF.h>
#include "bigmalloc.h"
#include "macr.h"
#include "protos.h"

int
main(int argc, char *argv[])
{
   SDF *sdfp;
   FILE *outfp;
   char *infile;
   short image;
   int nx, ny;
   int ret;

   infile = argv[1];

   Fopen(outfp, "topo_bs.img", "w");
   sdfp = SDFopen("img.hdr", infile);
   
   if( SDFgetint(sdfp, "nx", &nx) )
     Error("Sorry, you've got to have an \"nx\"\n");
   
   if( SDFgetint(sdfp, "ny", &ny) )
     Error("Sorry, you've got to have an \"ny\"\n");
   
   for (i = 0; i < nx*ny; i++) {
     ret = SDFseekrdvecs(sdfp, "image", i, 1, &image, NULL);
     if (ret) 
       Error("SDFseekrdvecs failed, %s\n", SDFerrstring);
     Fwrite(image, sizeof(short), 1, outfp);
   }

   Fclose(outfp);

   exit(0);
}
