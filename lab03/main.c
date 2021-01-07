#include "functions.c"
#include <stdlib.h>

int main(){

    /*Utworzenie wskaznikow do plikow, w ktorych zapisane beda dane*/
    double TOL;
    FILE *zad1 = fopen("zad1.dat", "w");
    FILE *zad2 = fopen("zad2.dat", "w");

    /*RK2 dla dwoch roznych wartosci TOL*/
    TOL = pow(10,-2);
    kontrola_kroku_czasowego(TOL, 1, zad1);
    fprintf(zad1, "\n\n");
    TOL = pow(10,-5);
    kontrola_kroku_czasowego(TOL, 1, zad1);

    /*Metoda Trapezow dla dwoch roznych wartosci TOL*/
    TOL = pow(10,-2);
    kontrola_kroku_czasowego(TOL, 2, zad2);
    fprintf(zad2, "\n\n");
    TOL = pow(10,-5);
    kontrola_kroku_czasowego(TOL, 2, zad2);

    fclose(zad1);
    fclose(zad2);

    return 0;
}