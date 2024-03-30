#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
using namespace std;

//Función a la que es igual la segunda derivada de R: D2R=-DR/r-R
double D2R(double r, double R, double DR){
  return -DR/r-R;
}

int main() {
  //Abrimos el archivo .dat
  ofstream outfile;
  outfile.open("datos_a.dat");

  //Elegimos dr
  double dr=0.01;
  
  //Dominio: [r0, rf]
  double r0=0.01;
  double rf=10.0;
  
  //Condiciones iniciales (aproximadas)
  double R0=1.0; 
  double DR0=0.0;
    
  //Declaramos las variables
  double r, R, DR;  //R(r), DR(r)

  /*En este loop inicializamos las variables, y evaluamos los diferentes valores de r mientras
  resolvemos la ecuación diferencial aplicando RK4*/
  for(r=r0, R=R0, DR=DR0; r<rf;){   
    //Mandamos los datos (r, R(r), DR(r))
    outfile << r << " " << R << " " << DR << endl;

    //RK4 para sistemas acoplados
    //RK4 para la primera derivada        //RK4 para la segunda derivada
    double K1_R=dr*DR;                    double K1_DR=dr*D2R(r, R, DR);
    double K2_R=dr*(DR+K1_DR/2);          double K2_DR=dr*D2R(r+dr/2, R+K1_R/2, DR+K1_DR/2);
    double K3_R=dr*(DR+K2_DR/2);          double K3_DR=dr*D2R(r+dr/2, R+K2_R/2, DR+K2_DR/2);
    double K4_R=dr*(DR+K3_DR);            double K4_DR=dr*D2R(r+dr, R+K3_R, DR+K3_DR);
    R=R+(K1_R+2*K2_R+2*K3_R+K4_R)/6.0;    DR=DR+(K1_DR+2*K2_DR+2*K3_DR+K4_DR)/6.0;   
    
    //Variación de r
    r=r+dr;   
  }   
  
  //Cerramos el archivo .dat
  outfile.close();
  
  return 0; 
}

//La gráfica fue hecha con Gnuplot
