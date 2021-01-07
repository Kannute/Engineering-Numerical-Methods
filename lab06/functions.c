#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "mgmres.c"

/*********** Definicje funkcji ***********/

int calc_j(int, int);
int calc_i(int, int);
double calc_epsl(int, int, int ,int);
double calc_ro1(double, double, double, double, double);
double calc_ro2(double, double, double, double, double);
int BC_Dirichlet(double, double, double, double *, int *, int *, double *,int nx, int, int, int, int, int, int, int);
int check_Dirichlet(FILE * , FILE * , double, double, double, double *, int *, int *, double *,int nx, int, int, int, int, int, int, int);
void check_poisson(FILE *, FILE *, double, double, double, int, int, int, int, int, int, int, int);
void algebrae_poisson(FILE * , double, double, double, int, int, int, int, int, int, int, int);



/*********** Implementacja ***********/

//Funkcja obliczajaca j wg. wzoru (12)
int calc_j(int nx, int l){

    return (l/(nx+1));

}

//Funkcja obliczajaca j wg. wzoru (13)
int calc_i(int nx, int l){

    return l-calc_j(nx,l)*(nx+1);

}

//Funkcja obliczajaca epsilon_l wg. wzoru (21)
double calc_epsl(int eps1, int eps2, int nx, int l){

    if(calc_i(nx,l) <= nx/2)
        return eps1;
    else
        return eps2;

}

//Funkcja obliczajca ro(1)(x,y) wg. wzoru (25)
double calc_ro1(double x, double y, double xmax, double ymax, double sigma){

    return exp(-1*pow(x - 0.25*xmax, 2)/pow(sigma, 2) - pow(y - 0.5*ymax, 2)/pow(sigma, 2));

}

//Funkcja obliczajaca ro(2)(x,y) wg. wzoru (26)
double calc_ro2(double x, double y, double xmax, double ymax, double sigma){
    
    return -1*exp(-1*pow(x - 0.75*xmax, 2)/pow(sigma, 2) - pow(y - 0.5*ymax, 2)/pow(sigma, 2));

}

//Funkcja rozwiacujaca algebraicznie row. Poissona
void algebrae_poisson(FILE * map, double delta, double xmax, double ymax, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4){

    
    int N = (nx+1)*(ny+1);   
    
    double *a = malloc((5*N) * sizeof(double));
    int *ja = malloc((5*N) * sizeof(int));
    int *ia = malloc((N + 1) * sizeof(int)); 
    double *b = malloc(N * sizeof(double));
    double *V = malloc(N * sizeof(double));
   

    int nz_num = BC_Dirichlet(delta, xmax, ymax, a, ia, ja, b, nx, ny, eps1, eps2, V1, V2, V3, V4);

    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10,-8);
    double tol_rel = pow(10,-8);

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);


    double temp = 0.0;

    //Zapisywanie danych do stworzenia mapy potencjalu
	for(int l = 0; l < N; ++l){

		if( delta*calc_j(nx,l) > temp)
			fprintf(map,"\n");

		fprintf(map,"%f %f %f \n", delta*calc_j(nx,l), delta*calc_i(nx,l), V[l]);
		temp = delta*calc_j(nx,l);
	}

    free(ja);
    free(ia);
    free(a);
    free(b);
    free(V);
}

//Warunki Brzegowe Dirichleta wg algorytmu przedstawionego w sekcji (3)
int BC_Dirichlet(double delta, double xmax, double ymax, double *a, int *ia, int *ja, double *b,int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4){

    
	int k = -1;

    int N = (nx+1)*(ny+1);

	for(int l = 0; l < N; ++l) {

		int brzeg = 0;  
		double vb = 0;  

		if(calc_i(nx, l) == 0){
			brzeg = 1;
			vb = (double)V1;
		}

		else if(calc_i(nx, l) == nx){
			brzeg = 1;
			vb = (double)V3;
		}

		else if(calc_j(nx, l) == ny){
			brzeg = 1;
			vb = (double)V2;
		}

		else if(calc_j(nx, l) == 0){
			brzeg = 1;
			vb = (double)V4;
		}

        b[l] = (-1)*(calc_ro1(delta*calc_i(nx,l), delta*calc_j(nx,l),xmax,ymax,xmax/10) + calc_ro2(delta*calc_i(nx,l), delta*calc_j(nx,l),xmax,ymax,xmax/10)); 


		if(brzeg == 1)
			b[l] = vb;


		ia[l] = -1;

		if(l - nx - 1 > 0 && brzeg == 0){
			k++;
			if(ia[l] < 0)
                ia[l] = k;

			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l - nx - 1;
		}

		if(l-1 > 0 && brzeg == 0) {
			k++;
			if(ia[l] < 0)
                ia[l] = k;
			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l-1;
		}

		k++;
		if(ia[l] < 0)
            ia[l] = k;

		if(brzeg == 0)
			a[k] = -(2*calc_epsl(nx,l,eps1,eps2)+calc_epsl(nx,l+1,eps1,eps2)+calc_epsl(nx,l+nx+1,eps1,eps2))/(delta*delta);
		else
			a[k] = 1;

		ja[k] = l;

		if(l < N && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+1,eps1,eps2)/(delta*delta);
			ja[k] = l + 1;
		}

		if(l < N-nx-1 && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+nx+1,eps1,eps2)/(delta*delta);
			ja[k] = l + nx + 1;
		}

    }

	ia[N] = k+1;
    return k+1;
}

//Funkcja sluzaca do sprawdzenia poprawnosc wypelnienia macierzy i wektora
void check_poisson(FILE * matrix, FILE * vector, double delta, double xmax, double ymax, int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4){
    int N = (nx+1)*(ny+1);
    double a[5*N];
    int ja[5*N];
    int ia[N+1];
    double b[N];
    double V[N];

    int nz_num = check_Dirichlet( matrix, vector, delta, xmax, ymax, a, ia, ja, b, nx, ny, eps1, eps2, V1, V2, V3, V4);

    int itr_max = 500;
    int mr = 500;
    double tol_abs = pow(10,-8);
    double tol_rel = pow(10,-8);

    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);
}

//Funkcja sluzaca do sprawdzenia poprawnosc wypelnienia macierzy i wektora wpisujaca dane do plikow
int check_Dirichlet(FILE * matrix, FILE * vector, double delta, double xmax, double ymax, double *a, int *ia, int *ja, double *b,int nx, int ny, int eps1, int eps2, int V1, int V2, int V3, int V4){


	int k = -1;

    int N = (nx+1)*(ny+1);

	for(int l = 0; l < N; ++l) {

		int brzeg = 0;  
		double vb = 0;  

		if(calc_i(nx, l) == 0){
			brzeg = 1;
			vb = (double)V1;
		}

		else if(calc_i(nx, l) == nx){
			brzeg = 1;
			vb = (double)V3;
		}

		else if(calc_j(nx, l) == ny){
			brzeg = 1;
			vb = (double)V2;
		}

		else if(calc_j(nx, l) == 0){
			brzeg = 1;
			vb = (double)V4;
		}

        b[l] = (-1)*(calc_ro1(delta*calc_i(nx,l), delta*calc_j(nx,l),xmax,ymax,xmax/10) + calc_ro2(delta*calc_i(nx,l), delta*calc_j(nx,l),xmax,ymax,xmax/10)); 


		if(brzeg == 1)
			b[l] = vb;


		ia[l] = -1;

		if(l - nx - 1 > 0 && brzeg == 0){
			k++;
			if(ia[l] < 0)
                ia[l] = k;

			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l - nx - 1;
		}

		if(l-1 > 0 && brzeg == 0) {
			k++;
			if(ia[l] < 0)
                ia[l] = k;
			a[k] = calc_epsl(nx,l,eps1,eps2)/(delta*delta);
			ja[k] = l-1;
		}

		k++;
		if(ia[l] < 0)
            ia[l] = k;

		if(brzeg == 0)
			a[k] = -(2*calc_epsl(nx,l,eps1,eps2)+calc_epsl(nx,l+1,eps1,eps2)+calc_epsl(nx,l+nx+1,eps1,eps2))/(delta*delta);
		else
			a[k] = 1;

		ja[k] = l;

		if(l < N && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+1,eps1,eps2)/(delta*delta);
			ja[k] = l + 1;
		}

		if(l < N-nx-1 && brzeg == 0){
			k++;
			a[k] = calc_epsl(nx,l+nx+1,eps1,eps2)/(delta*delta);
			ja[k] = l + nx + 1;
		}
        if(l%5 == 0 && l != 0)
            fprintf(vector, "\n");
        fprintf(vector, "%d %d %d %f \n", l, calc_i(nx,l), calc_j(nx,l), b[l]);

    }


	ia[N] = k+1;
    for(int i = 0; i < 5*N; i++)
        fprintf(matrix,"%d %0.f \n", i, a[i]);

    return k+1;
}