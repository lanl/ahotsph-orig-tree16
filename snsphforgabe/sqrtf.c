/*
 * © 2026. Triad National Security, LLC. All rights reserved.
 * This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.
 */


float sqrtf_fast(float x) {
    float __value, __arg = (x);
    asm("fsqrts %0,%1" : "=f"(__value) : "f"(__arg));
    return __value;
}


float recip8bit(float x) {
    float __value, __arg = (x);
    asm("fres %0,%1" : "=f"(__value) : "f"(__arg));
    return __value;
}


/* Really 5bit */
float recipsqrt8bit(float x) {
    float __value, __arg = (x);
    asm("frsqrte %0,%1" : "=f"(__value) : "f"(__arg));
    return __value;
}
