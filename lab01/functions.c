#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*Metoda jawna Eulera*/
void MJ_Euler(double, double, double, double, double*);

/*Metoda jawna RK2 (trapezów)*/
void MJ_RK2(double, double, double, double, double*);

/*Metoda jawna RK4*/
void MJ_RK4(double, double, double, double, double*);

/*Rozwiazanie rownania RRZ rzedu 2 z pomoca metody jawnej RK4*/
void RRZ_Tier2(double, double, double, double);



/******** FUNKCJE POMOCNICZE ********/

/*   Obliczenie y(t) = e^(lamda*t)    */
double calc_y(double lambda, double t){
    return exp(lambda*t);
}

/*  Obliczenie V(t) */
double calc_V(double omega, double t){
    return 10*sin(omega*t);
}

/*  Obliczenie g(t,Q,I) */
double calc_g(double R, double L, double C, double Q, double I, double t, double omega){
    return (calc_V(omega, t)/L - (R/L)*I - Q/(L*C));
}



/**************     IMPLEMENTACJA METOD     ***************/

void MJ_Euler(double _y, double lambda, double t_min, double t_max, double *steps){
    
    FILE *fp1, *fp2;
    fp1 = fopen("zad1_1.dat", "w");
    fp2 = fopen("zad1_2.dat", "w");

    for(int i=0;i<3;i++)
    {
        double y[2] = {_y, _y};
        for(double t = t_min; t <= t_max; t += steps[i])
        {
            //printf("%lf %lf \n", t, y[0]);
            //printf("%lf %lf \n", t, y[0]-calc_y(lambda, t));

            fprintf(fp1, "%lf %lf\n", t, y[0]);
            fprintf(fp2, "%lf %lf\n", t, y[0] -calc_y(lambda, t));

            y[0] = y[1] + steps[i] * lambda * y[1];
            y[1] = y[0];
        }
        fprintf(fp1, "\n\n");
        fprintf(fp2, "\n\n");

    }

    for(double t = t_min; t <= t_max; t += steps[0])
        fprintf(fp1, "%lf %lf\n",t, calc_y(lambda,t));
    
 
    fclose(fp1);
    fclose(fp2);
}

void MJ_RK2(double _y, double lambda, double t_min, double t_max, double *steps){

    FILE *fp1, *fp2;
    fp1 = fopen("zad2_1.dat", "w");
    fp2 = fopen("zad2_2.dat", "w");

    for(int i=0;i<3;i++)
    {
        double y[2] = {_y, _y};
        for(double t = t_min; t <= t_max; t += steps[i])
        {
            //printf("%f %f \n", t, y[0]);
            //printf("%f %f \n", t, y[0]-calc_y(lambda, t));

            fprintf(fp1, "%f %f\n", t, y[0]);
            fprintf(fp2, "%f %f\n", t, y[0]-calc_y(lambda, t));

            double k1 = lambda * y[0];
            double k2 = lambda * (y[0] + steps[i] * k1);
            y[0] = y[1] + (steps[i]/2.0)*(k1+k2);
            y[1] = y[0];
        }
        fprintf(fp1, "\n\n");
        fprintf(fp2, "\n\n");
    }
    for(double t = t_min; t <= t_max; t += steps[0])
        fprintf(fp1, "%lf %lf\n",t, calc_y(lambda,t));
    

    fclose(fp1);
    fclose(fp2);
}

void MJ_RK4(double _y, double lambda, double t_min, double t_max, double *steps){

    FILE *fp1, *fp2;
    fp1 = fopen("zad3_1.dat", "w");
    fp2 = fopen("zad3_2.dat", "w");

    for(int i=0;i<3;i++)
    {
        double y[2] = {_y, _y};
        for(double t = t_min; t <= t_max; t += steps[i])
        {
            //printf("%lf %lf \n", t, y[0]);
            //printf("%lf %lf \n", t, y[0]-calc_y(lambda, t));

            fprintf(fp1,"%lf %lf\n", t, y[0]);
            fprintf(fp2,"%lf %lf\n", t, y[0]-calc_y(lambda, t));

            double k1 = lambda * y[0];
            double k2 = lambda * (y[0] + (steps[i] / 2.0) * k1);
            double k3 = lambda * (y[0] + (steps[i] / 2.0) * k2);
            double k4 = lambda * (y[0] + steps[i] * k3);

            y[0] = y[1] + (steps[i]/6.0)*(k1+ (2*k2) + (2*k3) + k4);
            y[1] = y[0];
        }
        fprintf(fp1, "\n\n");
        fprintf(fp2, "\n\n");
    }
    
    for(double t = t_min; t <= t_max; t += steps[0])
        fprintf(fp1, "%lf %lf\n",t, calc_y(lambda,t));
    
    fclose(fp1);
    fclose(fp2);
}


void RRZ_Tier2(double R, double L, double C, double delta_t){

    FILE *fp1, *fp2;
    fp1 = fopen("zad4_1.dat", "w");
    fp2 = fopen("zad4_2.dat", "w");

    double omega_0 = 1.0/sqrt(L*C);
    double T_0 = 2*3.14159265/omega_0;
    double t_min = 0.0;
    double t_max = 4.0*T_0;
    double omega_V[4] = {0.5*omega_0, 0.8*omega_0, 1.0*omega_0, 1.2*omega_0};

    for(int i=0; i<4; i++)
    {
        double Q[2] = {0.0,0.0};
        double I[2] = {0.0,0.0};
        for(double t=t_min; t<t_max; t+=delta_t)
        {
            //printf("%lf %lf\n", t, Q[1]);
            //printf("%lf %lf\n", t, I[1]);

            fprintf(fp1,"%lf %lf\n", t, Q[1]);
            fprintf(fp2,"%lf %lf\n", t, I[1]);


            double kQ1 = I[0];
            double kI1 = calc_g(R,L,C,Q[0],I[0], t, omega_V[i]);

            double kQ2 = I[0] + (delta_t/2.0)*kI1;
            double kI2 = calc_g(R,L,C,Q[0]+(delta_t/2.0)*kQ1,I[0]+(delta_t/2.0)*kI1,t+delta_t/2.0,omega_V[i]);

            double kQ3 = I[0] + (delta_t/2.0)*kI2;
            double kI3 = calc_g(R,L,C,Q[0]+(delta_t/2.0)*kQ2,I[0]+(delta_t/2.0)*kI2,t+delta_t/2.0,omega_V[i]);

            double kQ4 = I[0] + (delta_t/2.0)*kI3;
            double kI4 = calc_g(R,L,C,Q[0]+(delta_t/2.0)*kQ3,I[0]+(delta_t/2.0)*kI3,t+delta_t/2.0,omega_V[i]);

            Q[1] = Q[0] + (delta_t/2.0)*(kQ1 + 2*kQ2 + 2*kQ3 + kQ4);
            I[1] = I[0] + (delta_t/2.0)*(kI1 + 2*kI2 + 2*kI3 + kI4);

            Q[0] = Q[1];
            I[0] = I[1];
        }
        fprintf(fp1, "\n\n");
        fprintf(fp2, "\n\n");
    }

    
    fclose(fp1);
    fclose(fp2);
}