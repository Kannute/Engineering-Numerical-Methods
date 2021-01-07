#include "functions.c"
#include <stdlib.h>
#include <stdio.h>

int main(){

    double omega_g;

    omega_g = 1.0;
    FILE * integral = fopen("integral_1.0_relglob.dat", "w");
    FILE * map = fopen("map_1.0_relglob.dat", "w");
    FILE * error = fopen("error_1.0_relglob.dat", "w");
    
    relaksacja_glob(omega_g, integral, map, error);

    fclose(integral);
    fclose(map);
    fclose(error);


    omega_g=0.6;
    integral = fopen("integral_0.6_relglob.dat", "w");
    map = fopen("map_0.6_relglob.dat", "w");
    error = fopen("error_0.6_relglob.dat", "w");

    relaksacja_glob(omega_g, integral, map, error);

    fclose(integral);
    fclose(map);
    fclose(error);


    double omega_l;

    omega_l = 1.0;
    integral = fopen("integral_1.0_relloc.dat", "w");
    relaksacja_lok(omega_l, integral);
    fclose(integral);

    omega_l = 1.4;
    integral = fopen("integral_1.4_relloc.dat", "w");
    relaksacja_lok(omega_l, integral);
    fclose(integral);

    omega_l = 1.8;
    integral = fopen("integral_1.8_relloc.dat", "w");
    relaksacja_lok(omega_l, integral);
    fclose(integral);

    omega_l = 1.9;
    integral = fopen("integral_1.9_relloc.dat", "w");
    relaksacja_lok(omega_l, integral);
    fclose(integral);


    return 0;
}
