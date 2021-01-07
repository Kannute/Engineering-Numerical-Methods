#include <stdlib.h>
#include <math.h>
#include <stdio.h>

/*** Stałe ***/
#define n_x 150
#define n_y 100
#define delta 0.1

double eps = 1;
double V1 = 10;
double V2 = 0;

const double TOL = pow(10,-8);
const double x_max  = delta * n_x;
const double y_max = delta * n_y;

double sigma_x = 0.1 * 15;
double sigma_y = 0.1 * 10;

/* Deklaracja tablicy gestosci */
double ro[n_x+1][n_y+1];

/*** Deklaracja funkcji ***/
double warunek_stopu(double V[][n_y + 1]);
void fill_ro();
void relaksacja_glob(double omega_g, FILE * integral, FILE * map, FILE * error);
void relaksacja_lok(double omega_l, FILE * integral);
void zero_tab(double V[][n_y + 1]);
void warunki_brzeg(double V[][n_y + 1]);
void copy_tabs(double V_s[][n_y + 1], double V_n[][n_y + 1]);
void rel_glob_loop(double omega_g, double V_s[][n_y + 1], double V_n[][n_y + 1]);
void rel_lok_loop(double omega_l, double V[][n_y + 1]);
double solution_error(double V[][n_y + 1], int i, int j);
/*****************************************IMPLEMENTACJA*******************************************/

//Wyzerowanie tablicy
void zero_tab(double V[][n_y + 1]){

    for(int i = 0 ; i <= n_x; i++){
        for(int j = 0; j <= n_y; j++){
            V[i][j] = 0.0;
        }
    }
}

//WB tablicy
void warunki_brzeg(double V[][n_y + 1]){

    for(int i = 0; i <= n_x; i++){
        V[i][0] = V1;
        V[i][n_y] = V2;
    }

}

//Kopiowanie tablicy V_n do tablicy V_s
void copy_tabs(double V_s[][n_y + 1], double V_n[][n_y + 1]){

    for(int i = 0 ; i <= n_x; i++)
            for(int j = 0; j <= n_y; j++)
                V_s[i][j] = V_n[i][j];
    
}

//Wypelnienie tablicy gestosci
void fill_ro(){

    for(int i = 0; i <= n_x; i++){
        for(int j = 0; j <= n_y; j++){
			ro[i][j] =  exp( -pow(delta*i - 0.35 * x_max, 2) / (pow(sigma_x,2)) - pow(delta*j - 0.5 * y_max, 2) / (pow(sigma_y,2)) );
            ro[i][j] += (-1.0)*exp( -pow(delta*i - 0.65 * x_max, 2) / (pow(sigma_x,2)) - pow(delta*j - 0.5 * y_max, 2) / (pow(sigma_y,2)) ) ;
            }
        }
}

//Warunek stopu
double warunek_stopu(double V[][n_y + 1]){

    double S = 0.0;

    for(int i = 0; i < n_x; i++) 
		for(int j = 0; j < n_y; j++) 
				S += pow(delta,2) * ( 0.5*pow((V[i+1][j] - V[i][j])/delta, 2) +  0.5*pow((V[i][j+1] - V[i][j])/delta, 2) - ro[i][j]*V[i][j]);
                	
    return S;
}
//Wykonywanie krokow dla relaksacji globalnej
void rel_glob_loop(double omega_g, double V_s[][n_y + 1], double V_n[][n_y + 1]){

        
    for(int i = 1; i < n_x; i++)
        for(int j = 1; j < n_y; j++)
            V_n[i][j] = 0.25 * ( V_s[i+1][j] + V_s[i-1][j] + V_s[i][j+1] + V_s[i][j-1] +  pow(delta,2)/eps * ro[i][j] );
        
        

    for(int j = 1; j < n_y; j++){
        V_n[0][j] = V_n[1][j];
        V_n[n_x][j] = V_n[n_x - 1][j];
    }

    for(int i = 0 ; i <= n_x ; i++)
        for(int j = 0; j <= n_y; j++)
            V_s[i][j] = (1-omega_g)*V_s[i][j] + omega_g*V_n[i][j];      

}
//Wykonywanie krokow dla relaksacji lokalnej
void rel_lok_loop(double omega_l, double V[][n_y + 1]){

    for(int i = 1; i < n_x; i++)
            for(int j = 1; j < n_y; j++)
                V[i][j] = (1 - omega_l)* V[i][j] +(omega_l/4.0) * ( V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] +  pow(delta,2)/eps * ro[i][j] );

    for(int j = 1; j < n_y; j++){
        V[0][j] = V[1][j];
        V[n_x][j] = V[n_x - 1][j];
    }

}

//Obliczenie bledu rozwiazania
double solution_error(double V[][n_y + 1], int i, int j){
    return ( (V[i+1][j] - 2*V[i][j] + V[i-1][j])/(pow(delta,2)) + (V[i][j+1] - 2*V[i][j] + V[i][j-1])/(pow(delta,2)) + ro[i][j]/eps );
}

/***  Relaksacja globalna dla danej omegi wpisujaca dane do plikow  ***/
void relaksacja_glob(double omega_g, FILE * integral, FILE * map, FILE * error){

    double s_it[2] = {0.0, 0.0};

    double V_n[n_x+1][n_y+1];
    double V_s[n_x+1][n_y+1];

    fill_ro();

    zero_tab(V_n);
    
    warunki_brzeg(V_n);

    copy_tabs(V_s, V_n);


    s_it[1] = warunek_stopu(V_n);
   
    int iteration = 0;

    do{
    
        iteration+=1;

        rel_glob_loop(omega_g, V_s, V_n);

        s_it[1] = s_it[0];
        s_it[0] = warunek_stopu(V_s);

        
        fprintf(integral, "%d %f \n", iteration, s_it[0]);

    }while(fabs( (s_it[0] - s_it[1])/s_it[1] ) > TOL);


    for(int i = 0; i <= n_x; i++){
        for(int j = 0; j <= n_y; j++){
            
            fprintf(map, "%f %f %f \n", delta*i, delta*j, V_s[i][j]);

            
            if(i > 0 && j > 0 && i < n_x && j < n_y)
                fprintf(error, "%f %f %f \n", delta*i, delta*j, solution_error(V_s, i, j) );
             else
                fprintf(error, "%f %f %f \n", delta*i, delta*j, 0.0);
        }
        fprintf(map, "\n");
        fprintf(error, "\n");

    }

    

}

/***  Relaksacja lokalna dla danej omegi wpisujaca dane do pliku  ***/
void relaksacja_lok(double omega_l, FILE * integral){

    double s_it[2] = {0.0, 0.0};

    double V[n_x+1][n_y+1];

    zero_tab(V);
    
    warunki_brzeg(V);

    s_it[1] = warunek_stopu(V);

    int iteration = 0;
    do{
        
        iteration+=1;

        rel_lok_loop(omega_l, V);        

        s_it[1] = s_it[0];
        s_it[0] = warunek_stopu(V);

        
        fprintf(integral, "%d %f \n", iteration, s_it[0]);
        
    }while(fabs( (s_it[0] - s_it[1])/s_it[1] ) > TOL);

}




