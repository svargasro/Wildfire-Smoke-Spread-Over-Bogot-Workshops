#include <iostream>
#include <cmath>
#include <random>
#include "../vector.h"

//Constantes del problema físico
double Lx=60, Ly=60; //Tamaño de la caja
const int Nx=5, Ny=5; 
const int N=Nx*Ny, Ntot=N+1; //Numero de partículas Nx*Ny y una más que servira de pared
const double e=1.0, r_0 = 10.0; //Constantes del potencial de Lennard-Jones
const double Rpared=50; //Radio de la pared circular
const double Mpared=0; //Masa de la pared circular (no tiene ningún efecto en la simulación por lo que puede tomar cualquier valor)

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
  vector3D r,V,F; double m,R; double theta,omega,tau,I;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0); tau=0;};// Inline
  void SumeFuerza(vector3D dF){F+=dF;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  vector3D Getr(void){return r;}; // Inline
  friend class Colisionador;
};
class Colisionador{
public:
  void CalculeTodasLasFuerzas(Cuerpo * Grano);
  void CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
  theta=theta0; omega=omega0; I=2.0/5.0*m*R*R;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt); theta+=omega*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m); omega+=tau*(coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
      //<<" , "<<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t"; //Dibuja una linea en cada grano para apreciar la rotación 
}
//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * Grano){
  int i,j;
  //Borro las fuerzas de todos los granos
  for(i=0;i<N+1;i++)
    Grano[i].BorreFuerza();

  //Fuerza de las paredes  
  vector3D Fpared, n;
  double r, r_centro, s;  
  int K = 10000;
  for(i=0;i<N;i++){
    r_centro = sqrt(pow(Grano[i].Getx(),2)+pow(Grano[i].Gety(),2)) + Grano[i].R;
    s = r_centro - Rpared;
    if(s>0){
      r = Grano[i].r.norm();
      n = Grano[i].r*(1.0/r);
      Fpared = (-K*s)*n;
      Grano[i].SumeFuerza(Fpared);
    }
  }
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=i+1;j<N;j++){
      CalculeFuerzaEntre(Grano[i],Grano[j]);};
}
void Colisionador::CalculeFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2){
  //Calcular el vector normal
  vector3D r21=Grano2.r-Grano1.r; 
  double r=r21.norm();
  //Calcular la fuerza
  double aux = (12*e/r)*(pow(r_0/r,12)-pow(r_0/r,6));
  vector3D n = r21*(1.0/r);  
  vector3D F1=n*aux;  
  //Sumar las fuerzas
  Grano1.SumeFuerza(-1*F1);  
  Grano2.SumeFuerza(F1);
}

//------------Funciones Globales ----------------
// double R_giro(){
//   vector3D r_cm;
//   for (int i=0; i<N; i++){
//     r_cm += Grano[i].r;
//   }


//   r_cm = r_cm*(1/N);
//   return R_giro
// }

int main(){
  Cuerpo Grano[Ntot];
  Colisionador LJ;
  int i,ix,iy;
  //Parametros de la simulación
  double m0=1.0; double R0=2.0; //Masas y radios de las particulas
  double KbT=10; //Constante de Boltzmann 
  double dx=10.0,dy=10.0; //Separación entre partículas
  double theta; double V0=sqrt(KbT/m0), omega0=3*M_PI;//omega es la velocidad angular inicial;
  double x0,y0,Vx0,Vy0;
  //Variables auxiliares para correr la simulacion
  int Ncuadros=100; double t,tdibujo,dt=5e-4,tmax=100,tcuadro=tmax/Ncuadros; 
  int timp;
  int i_prom = 0;
  double R_prom = 0.0;
  
  //INICIO
  //Inicializar la pared
  //-------------(x0,y0,Vx0,Vy0,theta0,omega0,    m0,    R0)
  Grano[N].Inicie( 0, 0,  0,  0,     0,     0,Mpared,Rpared); //Pared circular

  // Inicializar los granos
  for(ix=0;ix<Nx;ix++)
    for(iy=0;iy<Ny;iy++){
      std::mt19937_64 rng{std::random_device{}()}; // Inicializa el generador con una semilla aleatoria
      std::uniform_real_distribution<double> dist{0.0, 1.0};
      double random_number = dist(rng);
      theta=2*M_PI*dist(rng); //Direccion de la velocidad inicial	aleatoria
      x0=-20.0+ix*dx; y0=-20.0+iy*dy; Vx0=V0*cos(theta); Vy0=V0*sin(theta); //Posiciones y velocidades iniciales
      //--------------------(x0,y0,Vx0,Vy0,theta0,omega0,m0,R0)
      Grano[iy*Nx+ix].Inicie(x0,y0,Vx0,Vy0,     0,omega0,m0,R0);	
    }
      
  //CORRO
  for(t=timp=0;t<tmax;t+=dt,timp++){
    if(timp%100==0){      
      // for(i=0;i<N;i++) Grano[i].Dibujese(); 
      // tdibujo=0;
      vector3D r_cm;
      r_cm.load(0.0,0.0,0.0);
      for(i=0;i<N;i++){
        r_cm = (r_cm + Grano[i].Getr());  
      }   
      r_cm = r_cm*(1.0/N);
      double R_giro = 0;
      vector3D disr_cm;
      for(i=0;i<N;i++){
        disr_cm = Grano[i].Getr()-r_cm; 
        R_giro += disr_cm.norm2();
      }
      R_giro = sqrt(R_giro*(1.0/N));
      std::cout<<t<<"\t"<<R_giro<<std::endl;
      timp=0;
      if(t > 20){
        i_prom++;
        R_prom += R_giro;
      }
    }
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);    
    LJ.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    LJ.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++) Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,Um2chiplusxi);
    LJ.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,chi);
    LJ.CalculeTodasLasFuerzas(Grano); for(i=0;i<N;i++)Grano[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Grano[i].Mueva_r(dt,xi);     
  }
  R_prom = R_prom/i_prom;
  std::cout<<KbT<<std::endl;
  std::cout<<R_prom<<std::endl;
  return 0;
}



