#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "functions.c"

int main(){

    double delta = 0.1;
    int nx = 4;
    int ny = 4;
    int eps1 = 1;
    int eps2 = 1;
    int V1 = 10;
    int V3 = 10;
    int V2 = -10;
    int V4 = -10;
    double xmax = 0;
    double ymax = 0;
    

    FILE * vector = NULL;
    FILE * matrix = NULL;
    FILE * map = NULL;

    /** Sprawdzenie poprawnosci wypelnienia wektora i macierzy **/
    puts("------------Sprawdzenie poprawnosci wypelnienia wektora i macierzy------------");
    vector = fopen("vector_check.dat", "w");
    matrix = fopen("matrix_check.dat", "w");

    check_poisson(matrix, vector, delta, xmax, ymax, nx, ny, eps1, eps2, V1, V2, V3, V4);

    fclose(vector);
    fclose(matrix);

    /** Mapy potencjalu **/
    puts("------------Mapy Potencjalu------------");
    //podpunkt 1.a
    puts("podpunkt 1.a");
    nx = 50;
    ny = 50;
    
    
    map = fopen("map1a.dat", "w");
    algebrae_poisson(map, delta, xmax, ymax, nx, ny, eps1, eps2, V1, V2, V3, V4);
    fclose(map);

    //podpunkt 1.b
    puts("podpunkt 1.b");
    nx = 100;
    ny = 100;
   
    map = fopen("map1b.dat", "w");
    algebrae_poisson(map, delta, xmax, ymax, nx, ny, eps1, eps2, V1, V2, V3, V4);
    fclose(map);

    //podpunkt 1.c
    puts("podpunkt 1.c");
    nx = 200;
    ny = 200;
    
    map = fopen("map1c.dat", "w");
    algebrae_poisson(map, delta, xmax, ymax, nx, ny, eps1, eps2, V1, V2, V3, V4);
    fclose(map);


    /** Rozklady potencjalow **/
    puts("------------Rozklady potencjalow------------");
    nx = ny = 100;
	V1=V2=V3=V4=0;
	xmax = delta*nx;
	ymax = delta*ny;

    //podpunkt 2.a
    puts("podpunkt 2.a");

    eps1 = 1;
    eps2 = 1;

    map = fopen("map2a.dat", "w");
    algebrae_poisson(map, delta, xmax, ymax, nx, ny, eps1, eps2, V1, V2, V3, V4);
    fclose(map);

    //podpunkt 2.b
    puts("podpunkt 2.b");
    eps1 = 1;
    eps2 = 2;

    map = fopen("map2b.dat", "w");
    algebrae_poisson(map, delta, xmax, ymax, nx, ny, eps1, eps2, V1, V2, V3, V4);
    fclose(map);

    //podpunkt 2.c
    puts("podpunkt 2.c");
    eps1 = 1;
    eps2 = 10;

    map = fopen("map2c.dat", "w");
    algebrae_poisson(map, delta, xmax, ymax, nx, ny, eps1, eps2, V1, V2, V3, V4);
    fclose(map);

    
}