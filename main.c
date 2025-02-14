#include<stdlib.h>
#include<stdio.h>

#include "ptw.h"

int main(int argc, char* argv) {

    double dt = 0.1;
    params.theta = 0.1;
    params.p = 2.0;
    params.s0 = 0.02;
    params.sInf = 0.01;
    params.kappa = 0.3;
    params.lgamma = -12.0;
    params.y0 = 0.01;
    params.yInf = 0.003;
    params.y1 = 0.09;
    params.y2 = 0.7;

    ptw(&dt);
}