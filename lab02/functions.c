#include <stdlib.h>
#include <math.h>
#include <stdio.h>


/* deklaracja stalych */
const double beta = 0.001;
const double N = 500;
const double gmm = 0.1;
const int tmax = 100;
const double delta_t = 0.1;
const int mu_max = 20;
const int u0 = 1;
const double TOL = pow(10,-6);



/** deklaracja funkcji reprezentujacych metody numeryczne oraz funkcji pomocniczych **/

void metodaTrapezow_itPicarda();
void metodaTrapezow_itNewtona();
double f(double);


void metodaNiejawna_RK2();
double F1(double, double, double, double, double);
double F2(double, double, double, double, double);
double m11(double, double);
double m12(double, double);
double m21(double, double);
double m22(double, double);



/*************************************************************************************/
/*************************************************************************************/



/** implementacja funkcjii **/

double f(double u){
    return ((beta*N - gmm)*u - beta*pow(u,2));
}


void metodaTrapezow_itPicarda(){

    FILE *fp;
    fp = fopen("zad1.dat", "w");

    double t = 0.0;
    double u_n = u0;
    double u_n1 = u0;
    double u_mu = 0.0;
    int min_it = 0;
    int max_it = 0;

    fprintf(fp, "%lf %lf %lf \n", t, u_n1, N-u_n1);

    for(t=0.1; t<tmax; t+= delta_t){

        min_it = 0;
        u_n = u_n1;
        u_mu = 0.0;

        while(fabs(u_n1-u_mu) > TOL && min_it++ < mu_max){
            u_mu = u_n1;
            u_n1 = u_n + (delta_t/2.0) * (f(u_n)+f(u_mu));
            if(min_it>max_it)
                max_it = min_it;
        }

        fprintf(fp, "%lf %lf %lf \n", t, u_n1, N-u_n1);

    }

    fclose(fp);

}

void metodaTrapezow_itNewtona(){

    FILE *fp;
    fp = fopen("zad2.dat", "w");

    double t = 0.0;
    double u_n = u0;
    double u_n1 = u0;
    double u_mu = 0.0;
    int min_it = 0;
    int max_it = 0;

    fprintf(fp, "%lf %lf %lf \n", t, u_n1, N-u_n1);

    for(t=0.1; t<tmax; t+= delta_t){

        min_it = 0;
        u_n = u_n1;
        u_mu = 0.0;

        while(fabs(u_n1-u_mu) > TOL && min_it++ < mu_max){
            u_mu = u_n1;
            u_n1 = u_mu - (u_mu - u_n - (delta_t/2.0)*(f(u_n)+f(u_mu)))/
                        (1 - (delta_t/2.0)*(beta*N - gmm  - 2*beta*u_mu));
            if(min_it>max_it)
                max_it = min_it;
        }

        fprintf(fp, "%lf %lf %lf \n", t, u_n1, N-u_n1);

    }

    fclose(fp);
}



double F1(double u1, double u2, double u_n, double a11, double a12){
    return (u1 - u_n - delta_t*(a11*((beta*N - gmm)*u1-beta*pow(u1,2)) + a12*((beta*N - gmm)*u2 - beta*pow(u2,2))));
}

double F2(double u1, double u2, double u_n, double a21, double a22){
    return (u2 - u_n - delta_t*(a21*((beta*N - gmm)*u1-beta*pow(u1,2)) + a22*((beta*N - gmm)*u2 - beta*pow(u2,2))));
}

double m11(double u1, double a11){
    return (1-delta_t*a11*(beta*N - gmm - 2*beta*u1));
}
double m12(double u2, double a12){
    return (-1*delta_t*a12*(beta*N - gmm - 2*beta*u2));
}
double m21(double u1, double a21){
    return (-1*delta_t*a21*(beta*N - gmm - 2*beta*u1));
}
double m22(double u2, double a22){
    return (1-delta_t*a22*(beta*N - gmm - 2*beta*u2));
}


void metodaNiejawna_RK2(){
    FILE *fp = fopen("zad3.dat", "w");

    double a11 = 0.25;
    double a12 = 0.25 + sqrt(3)/6.0;
    double a21 = 0.25 + sqrt(3)/6.0;
    double a22 = 0.25;

    double b1 = 0.5;
    double b2 = 0.5;

    double u1 = 0.0;
    double u2 = 0.0;
    double u1_mu = 0.0;
    double u2_mu = 0.0;

    double t = 0.0;
    double u_n = u0;
    double u_n1 = u0;
    int min_it = 0;
    int max_it = 0;

    double delta_u1;
    double delta_u2;

    fprintf(fp, "%lf %lf %lf \n", t, u_n1, N-u_n1);

    for(t=0.1;t<tmax;t+=delta_t){

        min_it = 0;
        u_n = u_n1;
        u1_mu = 0.0;
        u2_mu = 0.0;
        u1 = u_n;
        u2 = u_n;

        while((fabs(u1 - u1_mu) > TOL || fabs(u2 - u2_mu) > TOL) && min_it++ < mu_max){
            u1_mu = u1;
            u2_mu = u2;

            delta_u1 = (F2(u1,u2,u_n,a21,a22)*m12(u2,a12) - F1(u1,u2,u_n,a11,a12)*m22(u2,a22))/
                        (m11(u1,a11)*m22(u2,a22) - m12(u2,a12)*m21(u1,a21));

            delta_u2 = (F1(u1,u2,u_n,a11,a12)*m21(u1,a21) - F2(u1,u2,u_n,a21,a22)*m11(u1,a11))/
                        (m11(u1,a11)*m22(u2,a22) - m12(u2,a12)*m21(u1,a21));

            u1 = u1_mu + delta_u1;
            u2 = u2_mu + delta_u2;

            if(min_it > max_it)
                max_it = min_it;
        }

        u_n1  = u_n + delta_t*(b1*f(u1) + b2*f(u2));
        fprintf(fp, "%lf %lf %lf \n", t, u_n1, N-u_n1);

    }

    fclose(fp);
}