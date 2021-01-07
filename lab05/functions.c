#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*Stale*/
#define delta  0.2
#define nx  128
#define ny  128
const double xmax = nx*delta;
const double ymax = ny*delta;
const double TOL = pow(10,-8);



/*Deklaracja funkcji*/
double warunek_stopu(double V[][ny+1], int k);
void fill_V(double V[][ny+1]);
void WB_V(double V[][ny+1]);
void zageszczanie_siatki(double V[][ny+1], int k);
void poisson_multigrid();


/*Implementacja funkcji*/
double warunek_stopu(double V[][ny+1], int k){
    double S = 0.0;

    for(int i = 0; i <= nx-k; i += k){
		for(int j = 0; j <= ny-k; j += k){
			S += 0.5*pow(k*delta, 2)*(pow((V[i+k][j] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i][j+k])/(2*k*delta), 2) +
					pow((V[i][j+k] - V[i][j])/(2*k*delta) + (V[i+k][j+k] - V[i+k][j])/(2*k*delta), 2));
		}
	}
    return S;
}


void poisson_multigrid(){

    FILE * map = fopen("maps.dat" , "w");
    FILE * integral = fopen("integral.dat", "w");

    double s[2] = {0.0, 0.0};
    int it= 0;
    double V[nx+1][ny+1];
    
    //Wypelnienie V zerami
    fill_V(V);

    //Warunki brzegowe
    WB_V(V);

    for(int k = 16; k>0; k=k/2){

        s[1] = warunek_stopu(V,k);

        do{
            //Dyskretyzacja siatki
            for(int i=k; i<=nx-k; i+=k)
                for(int j=k; j<= ny-k; j+=k)
                    V[i][j] = 0.25*(V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);

            s[0] = s[1];
            s[1] = warunek_stopu(V, k);
            fprintf(integral, "%d %d %f \n", k, it, s[1]);
            it++;

        }while(fabs((s[1] - s[0])/s[0]) > TOL);

        fprintf(integral, "\n \n");

        //mapowanie
        for(int i=0;i<nx+1;i+=k){
            for(int j=0;j<ny+1;j+=k){
                fprintf(map, "%f %f %f \n", delta*i, delta*j, V[i][j]);
            }
            fprintf(map,"\n");
        }
        fprintf(map, "\n \n");

        //zageszczenie siatki
        zageszczanie_siatki(V, k);

    }

    fclose(map);
    fclose(integral);
}

void fill_V(double V[][ny+1]){

    for(int i = 0 ; i < nx+1; i++)
        for(int j = 0; j < ny+1; j++)
            V[i][j] = 0.0;

}

void WB_V(double V[][ny+1]){

    for(int y=0; y<ny+1;y++)
        V[0][y] = (1.0) * sin(M_PI * delta*y / ymax);
    
    for(int x=0; x<nx+1;x++)
        V[x][ny] = (-1.0) * sin(2*M_PI * delta*x / xmax);

    for(int y=0;y<ny+1;y++)
        V[nx][y] = (1.0) * sin(M_PI * delta*y / ymax);

    for(int x=0;x<nx+1;x++)
        V[x][0] =  (1.0) * sin(2*M_PI* delta*x/xmax);

}

void zageszczanie_siatki(double V[][ny+1],int k){

    for(int i = 0; i <= nx-k; i += k){
		for(int j = 0; j <= ny-k; j += k){

            V[i + k/2][j + k/2] = 0.25*(V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);

            if(i!=nx-k)
                V[i + k][j + k/2] = 0.5*(V[i+k][j] + V[i+k][j+k]);
            if(j!=ny-k)
                V[i + k/2][j + k] = 0.5*(V[i][j+k] + V[i+k][j+k]);
            if(j!=0)
                V[i + k/2][j] = 0.5*(V[i][j] + V[i+k][j]);
            if(i!=0)
                V[i][j + k/2] = 0.5*(V[i][j] + V[i][j+k]);
        }
    }
}