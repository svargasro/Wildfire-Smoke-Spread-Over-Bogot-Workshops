#include <iostream>
#include <cmath>
#include "../../vector.h"
#include "../../Random64.h"
using namespace std;

//--------------- Constantes globales ------------
const double g=9.8, L=12;
const double K=1.0e4;

//Número de moleculass
const int N = 1;

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
  vector3D r, tau;
  double xs,ys;//punto de apoyo del pendulo
  double theta, omega;
  double m,R;
public:
  void Inicie(double xs, double ys, double theta0, double omega0, 
              double m0, double R0);
  void Mueva_theta(double dt,double coeficiente);
  void Mueva_omega(double dt,double coeficiente);
  void Dibujese(void);
  void BorreTorque(void){tau.load(0,0,0);};// Inline
  // void SumeTorque(double dtau){tau+=dtau;};
  void SumeTorque(vector3D dtau){tau+=dtau;};
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
  r.load(x0 + L*sin(theta),y0 - L*cos(theta),0);
}
#include <fstream>

void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(coeficiente*dt);

  ofstream archivo1("CunaNewtonTheta.txt", ios::app); // Open the file in append mode
  archivo1 << theta << endl; // Write the value of tau to the file
  archivo1.close(); // Close the file
  // clog<<"Theta = "<<theta<<endl;

  ofstream archivo2("CunaNewtonOmega.txt", ios::app); // Open the file in append mode
  archivo2 << omega << endl; // Write the value of tau to the file
  archivo2.close(); // Close the file
  // clog<<"Omega = "<<omega<<endl;
}
void Cuerpo::Mueva_omega(double dt,double coeficiente){
  double I_L = 1/3.0*m*L*L, I_R = 0.5*m*R*R;
  double I = I_L + I_R; 
  double alpha = tau.norm()/I;
  
  omega+=alpha*(coeficiente*dt);

  ofstream archivo3("CunaNewtonAlpha.txt", ios::app); // Open the file in append mode
  archivo3 << alpha << endl; // Write the value of tau to the file
  archivo3.close(); // Close the file
  // clog<<"Alpha = "<<alpha<<endl;
  
}
void Cuerpo::Dibujese(void){
  double x = L*sin(theta);
  double y = -L*cos(theta);
  cout<<" , "<<x + xs<<"+"<<R<<"*cos(t),"<<y + ys<<"+"<<R<<"*sin(t) , "
  <<xs<<"+"<<x/7.0<<"*t,"<<ys<<"+"<<y/7.0<<"*t";
}
//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodosLosTorques(Cuerpo * pendulos){
  int i;


  //Borro las fuerzas de todos los moleculass
  for(i=0;i<N;i++)
    pendulos[i].BorreTorque();

  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++){
    pendulos[i].BorreTorque();

    // double theta = pendulos[i].theta
    double m = pendulos[i].m;

    //Fuerza de la gravedad
    // double Fg = -m*g;
    // double xs = pendulos[i].xs,ys = pendulos[i].ys;
    // vector3D rL,F;
    // rL.load(L*sin(theta) + xs,-L*cos(theta) + ys,0);
    // F.load(0,Fg,0);
    //pendulos[i].SumeTorque(rL^F);
    
    vector3D F,r = pendulos[i].r;
    F.load(0,-m*g,0);
    pendulos[i].SumeTorque(r^F);
    
    // for(j=0;j<i;j++)
    //   CalculeTorqueEntre(pendulos[i],pendulos[j]);
  }

    
}
void Colisionador::CalculeTorqueEntre(Cuerpo & moleculas1, Cuerpo & moleculas2){
  
  //Calcular el vector normal
  vector3D r21=moleculas2.r-moleculas1.r; 
  double d=r21.norm();


  //Determinar si hay colisión
  double s= moleculas1.R + moleculas2.R - d;
  
  if(s>0){
    double mag= K*pow(s,1.5);

    //Calcular vector normal
    vector3D n = r21*(1.0/d);

    //Calcular la fuerza
    vector3D F=n*mag;

    //Calcular el torque
    vector3D r1 = n*moleculas1.R;
    vector3D r2 = n*moleculas2.R;
    vector3D rcm = r1 - r2;
    vector3D tau = rcm^F;
    
    moleculas2.SumeTorque(tau);  moleculas1.SumeTorque(tau*(-1));

  }
}


//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  cout<<"set terminal gif animate delay 20"<<endl; 
  cout<<"set output 'CunaNewtonVectorial.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-20:"<<Lx + 20<<"]"<<endl;
  cout<<"set yrange[-20:"<<Ly + 20<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(double tiempo){
    cout << "set label \'t=" << tiempo << "s \' at graph 0.5,0.9 center" << endl;
    cout<<"plot 0,0 ";;
    // cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    // cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    // cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    // cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
    
}
void TermineCuadro(void){
    cout<<endl;
    cout<<"unset label" << endl;
}

int main(){
  
  Cuerpo pendulos[N];
  Colisionador Newton;

  //Parametros de la simulación
  double m0=0.1; double R0=1.5;
  double theta01=-M_PI/12,theta0i=0;
  double T = 2*M_PI*sqrt(L/g);
  
  double x0=0,y0=L,omega0=0;
  double dx = 2*R0;//Separación entre pendulos 
 
  //Variables auxiliares para correr la simulacion
  int i, Ncuadros=100; 
  double t,tdibujo,dt=1e-3,tmax=3*T,tcuadro=tmax/Ncuadros;  

  InicieAnimacion();
  
  //Inicie las pendulos
  for(i=0;i<N;i++){
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
  
    
    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,xi);    
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++) pendulos[i].Mueva_omega(dt,Um2lambdau2);

    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,chi);
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++) pendulos[i].Mueva_omega(dt,lambda);
    
    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,Um2chiplusxi);
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++)pendulos[i].Mueva_omega(dt,lambda);

    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,chi);
    Newton.CalculeTodosLosTorques(pendulos); 
    for(i=0;i<N;i++)pendulos[i].Mueva_omega(dt,Um2lambdau2);

    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,xi);
    
    clog<<"\nt = "<<t<<endl;
  }
  return 0;
}
