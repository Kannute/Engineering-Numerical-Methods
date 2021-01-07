#include "functions.c"
#include <stdlib.h>
#include <math.h>

int main(){

    double* steps = malloc(3*sizeof(double));
    steps[0] = 0.01; 
    steps[1] = 0.1; 
    steps[2] = 1.0;

    MJ_Euler(1.0 ,-1.0, 0.0, 5.0, steps);
    
    MJ_RK2(1.0 ,-1.0, 0.0, 5.0, steps);

    MJ_RK4(1.0 ,-1.0, 0.0, 5.0, steps);
    

    double delta_t = pow(10, -4);

    RRZ_Tier2(100.0, 0.1, 0.001, delta_t);

    free(steps);

    return 0;
}