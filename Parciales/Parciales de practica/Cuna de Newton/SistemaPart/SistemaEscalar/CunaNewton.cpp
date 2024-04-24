/*
Cuna de Newton

Algoritmo de PERFL para mover varios péndulos simples bajo el efecto de la gravedad y 
la fuerza de Hertz para modelar la colisión entre ellos.

Obs:
- El paso de tiempo dt debe de ser mas pequeño a medida la constante de Hertz K sea mas grande.
- Se debe de tener cuidado con el valor de la constante de Hertz K, ya que si es muy grande
  los péndulos pueden colapsar.
*/


#include <iostream>
#include <cmath>
#include "../../../vector.h"
using namespace std;

//--------------- Constantes globales ------------
const double g=9.8, L=12;
double K=1.0e9;

//Número de moleculass
const int N = 3;

//Dimensiones de la caja
const double Lx=0, Ly=0; 



//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);

//--------------- Declarar las clases-----------
class Cuerpo;
class Colisionador;

//--------- Declarar las interfases de las clases---------
class Cuerpo{
private:
  double xs,ys;//punto de apoyo del pendulo
  double theta, omega, tau;
  double m,R;
public:
  void Inicie(double xs, double ys, double theta0, double omega0, 
              double m0, double R0);
  void Mueva_Theta(double dt,double coeficiente);
  void Mueva_Omega(double dt,double coeficiente);
  void Dibujese(void);
  void BorreTorque(void){tau = 0 ;};// Inline
  void SumeTorque(double dtau){tau+=dtau;};// Inline
  double GetX(void){return xs + L*sin(theta);};
  double GetY(void){return ys - L*cos(theta);};
  friend class Colisionador;  
};
class Colisionador{
private:
public:
  void CalculeTodosLosTorques(Cuerpo * pendulos);
  void CalculeTorqueEntre(Cuerpo & moleculas1,Cuerpo & moleculas2);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0, double y0, double theta0, double omega0, 
            double m0, double R0){
  xs=x0; ys=y0;
  theta=theta0; omega=omega0; 
  m=m0; R=R0;
}

void Cuerpo::Mueva_Theta(double dt,double coeficiente){
  theta+=omega*(coeficiente*dt);
}

void Cuerpo::Mueva_Omega(double dt,double coeficiente){
  double I_L = 1/3.0*m*L*L, I_R = 0.5*m*R*R;
  double I = I_L + I_R; 
  double alpha = tau / I;
  omega += alpha * (coeficiente * dt);
}

void Cuerpo::Dibujese(void){
  double x = L*sin(theta);
  double y = -L*cos(theta);
  cout<<" , "<<x + xs<<"+"<<R<<"*cos(t),"<<y + ys<<"+"<<R<<"*sin(t) , "
  <<xs<<"+"<<x/7.0<<"*t,"<<ys<<"+"<<y/7.0<<"*t";
}

//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodosLosTorques(Cuerpo * pendulos){
  int i,j;

  //Borro las fuerzas de todos los moleculass
  for(i=0;i<N;i++)
    pendulos[i].BorreTorque();

  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++){

    double theta = pendulos[i].theta;
    double m = pendulos[i].m;

    //Fuerza de la gravedad
    double tau = -m * g * L * sin(theta);
    pendulos[i].SumeTorque(tau);

    for(j=i+1;j<N;j++)
  		CalculeTorqueEntre(pendulos[i],pendulos[j]);
  }

  // for(i=0;i<N-1;i++)
  //   for(j=i+1;j<N;j++)
  //     CalculeTorqueEntre(pendulos[i],pendulos[j]);
}

void Colisionador::CalculeTorqueEntre(Cuerpo & pendulos1, Cuerpo & pendulos2){
    
    //distancia relativa entre los pendulos
    double d = pendulos2.GetX() - pendulos1.GetX();

    //distancia de interpentracion
    double s = pendulos1.R + pendulos2.R - d;
    

    if(s>0){
      double F = K*pow(s,1.5);
      double tau = F*L;
      pendulos1.SumeTorque(-tau);
      pendulos2.SumeTorque(tau);
      
    }
}


//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(double dt){
  cout<<"set terminal gif animate delay 25"<<endl; 
  cout<<"set output 'CunaNewton_dt="<<dt<<".gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-20:"<<Lx + 20<<"]"<<endl;
  cout<<"set yrange[-20:"<<Ly + 20<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(double t){
  cout << "set label \'t=" << t << "s \' at graph 0.5,0.9 center" << endl;
  cout<<"plot 0,0 ";

    // cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    // cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    // cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    // cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
    
}

void TermineCuadro(void){
    cout<<endl;
    cout<<"unset label"<<endl;
}

int main(int argc, char *argv[]){
  
  double dt = atof(argv[1]);
  
  Cuerpo pendulos[N];
  Colisionador Newton;

  //Parametros de la simulación
  double m0=0.1; double R0=1.5;
  double theta01=-M_PI/12,theta0i=0;
  double T = M_PI*(L/g);
  
  double x0=0,y0=L,omega0=0;
  double dx = 2*R0;//Separación entre pendulos 
 
  //Variables auxiliares para correr la simulacion
  int i, Ncuadros=100; 
  double t,tdibujo,tmax=5*T,tcuadro=tmax/Ncuadros; 

  InicieAnimacion(dt);
  
  //Inicie l0s pendulos
  //Inicializacion del primer pendulo
  pendulos[0].Inicie(x0,y0,3*theta01,omega0,m0,R0);
  
  //Inicializacion de los demas pendulos
  for(i=1;i<N;i++){
    x0=(i+1)*dx;
    pendulos[i].Inicie(x0,y0,theta0i,omega0,m0,R0);
  }
 
  
 
  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    
    if(tdibujo>=tcuadro){
      
      InicieCuadro(t);
      for(i=0;i<N;i++) pendulos[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }
  
    
    for(i=0;i<N;i++) pendulos[i].Mueva_Theta(dt,xi);    
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++) pendulos[i].Mueva_Omega(dt,Um2lambdau2);

    for(i=0;i<N;i++) pendulos[i].Mueva_Theta(dt,chi);
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++) pendulos[i].Mueva_Omega(dt,lambda);
    
    for(i=0;i<N;i++) pendulos[i].Mueva_Theta(dt,Um2chiplusxi);
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++)pendulos[i].Mueva_Omega(dt,lambda);

    for(i=0;i<N;i++) pendulos[i].Mueva_Theta(dt,chi);
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++)pendulos[i].Mueva_Omega(dt,Um2lambdau2);

    for(i=0;i<N;i++) pendulos[i].Mueva_Theta(dt,xi);
    
    clog<<"\nPorcentaje de avance: "<<(t/tmax)*100<<"%"<<endl;
  }
  return 0;
}
