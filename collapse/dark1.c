/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */

#include <stdio.h>
#include <stdlib.h>

#include "SDF.h"
#include "SDFwrite.h"
#include "physics.h"

void main(int argc, char *argv[]) {
    outbody b;

    b.mass = 1.0;

    b.pos[0] = 1.0;
    b.pos[1] = 1.0;
    b.pos[2] = 0.0;

    b.vel[0] = 0.0;
    b.vel[1] = -5.0;
    b.vel[2] = 0.0;

    b.ident = 11569;

    SDFwrite("dark1", 1, 1, &b, sizeof(outbody), OUTBODYDESC, "npart", SDF_INT, 1, NULL);

    exit(0);
}
