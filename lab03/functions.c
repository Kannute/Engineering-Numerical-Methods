#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Definicja makra zwracajacego maksymalna wartosc
#define max(x, y) (((x) > (y)) ? (x) : (y))

//Struktura reprezentujaca zwracana przez schematy funkcje pare liczb
struct vals{
    double x;
    double v;
};

//Definicja stalych
double x0 = 0.01;
double v0 = 0.0;
double delta_t0 = 1.0;
double s = 0.75;
double p = 2.0;
double t_max = 40.0;
double alfa = 5.0;

//Deklaracja funkcji przeprowadzjacej kontrole kroku czasowego
void kontrola_kroku_czasowego(double, int, FILE *);

//Deklaracja schematow numerycznych
struct vals schemat_numerycznyMT(double, double, double, double);
struct vals schemat_numerycznyRK2(double, double, double, double);

//Deklaracja funkcji pomocniczych
double f(double);
double g(double, double, double);
double F(double, double, double, double, double);
double G(double, double, double, double, double, double);
double a11();
double a12(double);
double a21(double, double, double, double);
double a22(double, double, double);



/******************************* implementacja metod *****************************/

double f(double v){
    return v;
}

double g(double x, double v, double alfa){
    return (alfa*(1 - pow(x,2))*v - x);
}

double F(double x_n1, double x_n, double v_n1, double v_n, double delta_t){
    return x_n1 - x_n - (delta_t/2.0) * ( f(v_n) + f(v_n1) );
}

double G(double x_n1, double x_n, double v_n1, double v_n, double delta_t, double alfa){
    return v_n1 - v_n  - (delta_t/2.0) * ( g(x_n,v_n, alfa) + g(x_n1, v_n1, alfa) );
}

double a11(){
    return 1.0;
}

double a12(double delta_t){
    return -1.0 * (delta_t/2.0);
}

double a21(double delta_t, double x_n1_k, double v_n1_k, double alfa){
    return (-1.0)*(delta_t/2.0)*( (-2.0*alfa * x_n1_k * v_n1_k) - 1 );
}


double a22(double delta_t, double x_n1_k, double alfa){
    return 1.0 - (delta_t/2.0)*alfa*(1-pow(x_n1_k,2));
}


struct vals schemat_numerycznyMT(double x_n, double v_n, double delta_t, double alfa){
    
    double delta = pow(10,-10);
    double x_n1 = x_n;
    double v_n1 = v_n;
    double x_n1_k = x_n;
    double v_n1_k = v_n;
    double delta_x;
    double delta_v;

    do{
        delta_x = 0.0;
        delta_v = 0.0;
        x_n1_k = x_n1;
        v_n1_k = v_n1;

        delta_x =  ((-1.0)*F(x_n1,x_n, v_n1, v_n, delta_t)*a22(delta_t,x_n1_k,alfa) - (-1.0)*G(x_n1,x_n,v_n1,v_n,delta_t, alfa)*a12(delta_t))
                                        /(a11()*a22(delta_t,x_n1_k,alfa) - a12(delta_t)*a21(delta_t,x_n1_k,v_n1_k,alfa));

        delta_v = (a11()*(-1.0*G(x_n1, x_n, v_n1, v_n, delta_t, alfa)) - a21(delta_t,x_n1_k, v_n1_k, alfa)*(-1.0*F(x_n1, x_n, v_n1, v_n, delta_t)))
                                        /(a11()*a22(delta_t,x_n1_k,alfa) - a12(delta_t)*a21(delta_t,x_n1_k,v_n1_k,alfa));

        x_n1 = x_n1_k + delta_x;
        v_n1 = v_n1_k + delta_v;

    }while(fabs(delta_x) < delta && fabs(delta_v) < delta);

    struct vals tmp;
    tmp.x = x_n + (delta_t/2.0) * (f(v_n) + f(v_n1));
    tmp.v = v_n + (delta_t/2.0 )* (g(x_n, v_n, alfa) + g(x_n1, v_n1, alfa));

    return tmp;

}


struct vals schemat_numerycznyRK2(double x_n, double v_n, double delta_t, double alfa){

    double k_1x = f(v_n);
    double k_1v = g(x_n, v_n, alfa);
    double k_2x = f(v_n+delta_t*k_1v);
    double k_2v = g(x_n + delta_t*k_1x, v_n + delta_t*k_1v, alfa);

    struct vals tmp;
    tmp.x = x_n + (delta_t/2.0) * (k_1x + k_2x);
    tmp.v = v_n + (delta_t/2.0) * (k_1v + k_2v);
    return tmp;
}
void kontrola_kroku_czasowego(double TOL, int flag, FILE *fp){


    double E_x = 0.0;
    double E_v = 0.0;
    struct vals temp1;
    struct vals temp2;

    double t = 0;
    double delta_t = delta_t0;
    double x_n = x0;
    double v_n = v0;

        
    if(flag == 1)
    {
        fprintf(fp, "%lf %lf %lf %lf \n", t, delta_t, x_n, v_n);

        do{
            temp2 = schemat_numerycznyRK2(x_n, v_n, delta_t, alfa);
            temp2 = schemat_numerycznyRK2(temp2.x, temp2.v, delta_t, alfa);

        
            temp1 = schemat_numerycznyRK2(x_n, v_n, 2.0*delta_t, alfa);

        
            E_x = (temp2.x - temp1.x)/(pow(2,p) - 1);
            E_v = (temp2.v - temp1.v)/(pow(2,p) - 1);

        if(max(fabs(E_x),fabs(E_v)) < TOL ){
            t += 2*delta_t;
            x_n = temp2.x;
            v_n = temp2.v;
            fprintf(fp, "%f %f %f %f \n", t, delta_t, x_n, v_n);

        }
        delta_t = (pow(s*TOL/(max(fabs(E_x),fabs(E_v))),(1.0/(p+1.0)))*delta_t);

        }while(t < (t_max - delta_t));

    }

    else if( flag == 2)
    {

    fprintf(fp, "%lf %lf %lf %lf \n", t, delta_t, x_n, v_n);

    do{
        temp2 = schemat_numerycznyMT(x_n, v_n, delta_t, alfa);
        temp2 = schemat_numerycznyMT(temp2.x, temp2.v, delta_t, alfa);

        
        temp1 = schemat_numerycznyMT(x_n, v_n, 2.0*delta_t, alfa);

        
        E_x = (temp2.x - temp1.x)/(pow(2,p) - 1);
        E_v = (temp2.v - temp1.v)/(pow(2,p) - 1);

        if(max(fabs(E_x),fabs(E_v)) < TOL ){
            t += 2*delta_t;
            x_n = temp2.x;
            v_n = temp2.v;
            fprintf(fp, "%f %f %f %f \n", t, delta_t, x_n, v_n);

        }
        delta_t = pow((s*TOL)/(max(fabs(E_x),fabs(E_v))),(1.0/(p+1.0))) * delta_t;

    }while(t < (t_max - delta_t));

 
    }

  
}