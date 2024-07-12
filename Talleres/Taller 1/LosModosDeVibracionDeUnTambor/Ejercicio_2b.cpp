#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
using namespace std;

//Función a la que es igual la segunda derivada de R: D2R=-DR/r-(lambda)^2*R
double D2R(double lambda, double r, double R, double DR){
  return -DR/r-lambda*lambda*R;
}

int main() {
  //Abrimos el archivo .dat
  ofstream outfile;
  outfile.open("datos_b.dat");

  //Elegimos un diferencial para lambda (dl) y para r (dr)
  double dl=0.01;
  double dr=0.01;
  
  //Dominio de lambda: [l0, lf]
  double l0=0.1;
  double lf=15.0;
  
  //Dominio de r: [r0, rf]
  double r0=0.01;
  double rf=1.0;  //rf es uno porque solo nos interesa la función R(r=1, lambda) 
  
  //Condiciones iniciales (aproximadas)
  double R0=1.0;
  double DR0=0.0;
  
  //Declaramos las variables    
  double lambda, r, R, DR;  //R(r, lambda), DR(r, lambda)

  //Este loop evalúa los diferentes valores de lambda
  for(lambda=l0; lambda<lf;){
    /*En este loop inicializamos las variables, y evaluamos los diferentes valores de r mientras
    resolvemos la ecuación diferencial aplicando RK4*/
    for(r=r0, R=R0, DR=DR0; r<rf;){
      
      //RK4 para sistemas acoplados
      //RK4 para la primera derivada        //RK4 para la segunda derivada
      double K1_R=dr*DR;                    double K1_DR=dr*D2R(lambda, r, R, DR);
      double K2_R=dr*(DR+K1_DR/2);          double K2_DR=dr*D2R(lambda, r+dr/2, R+K1_R/2, DR+K1_DR/2);
      double K3_R=dr*(DR+K2_DR/2);          double K3_DR=dr*D2R(lambda, r+dr/2, R+K2_R/2, DR+K2_DR/2);
      double K4_R=dr*(DR+K3_DR);            double K4_DR=dr*D2R(lambda, r+dr, R+K3_R, DR+K3_DR);
      R=R+(K1_R+2*K2_R+2*K3_R+K4_R)/6.0;    DR=DR+(K1_DR+2*K2_DR+2*K3_DR+K4_DR)/6.0;   
      r=r+dr;   
    }

    //Mandamos los datos (lambda, R(r=1, lambda))	 
    outfile << lambda << " " << R << endl;
    
    //Variación de lambda
    lambda=lambda+dl;
  }
  
  //Cerramos el archivo .dat
  outfile.close();
  
  return 0; 
}
//Los puntos donde la gráfica se hace cero son aproximadamente: 2.3, 5.6, 8.8, 11.9, 15
