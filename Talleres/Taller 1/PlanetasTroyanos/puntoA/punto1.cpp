#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"

using std::cout;
using std::sqrt;
using std::endl;
using std::fabs;

//Constantes del problema físico
const int N=2; //Número de cuerpos.
const double G=1.0;

//Constantes del algoritmo de integración: (PEFRL)
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);


//--------------- Declarar las clases-----------
class Cuerpo;
class Interaction;



//--------- Declarar las interfases de las clases---------
class Cuerpo{
private:
  vector3D r,V,F; double m;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0);
  void BorreFuerza(void){F.load(0,0,0);}; // Inline
  void SumeFuerza(vector3D dF){F+=dF;}; // Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  double Getx(void){return r.x();}; // Inline
  double Gety(void){return r.y();}; // Inline
  friend class Interaction;
};

class Interaction{
private:
public:
  void CalculeTodasLasFuerzas(Cuerpo * Planetas);
  void CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);
};

//-------Implementar las funciones de las clases------


//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double z0,
                    double Vx0,double Vy0,double Vz0,double m0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0;
}

void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}

void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}

//------- Funciones de la clase Interaction --------
void Interaction::CalculeTodasLasFuerzas(Cuerpo * Planeta){
  int i,j;
  //Borro las fuerzas de todos los planetas
  for(i=0;i<N;i++)
    Planeta[i].BorreFuerza();
  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++)
    for(j=0;j<i;j++)
      CalculeFuerzaEntre(Planeta[i],Planeta[j]);
}
void Interaction::CalculeFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  double m1=Planeta1.m, m2=Planeta2.m;
  vector3D r21=Planeta2.r-Planeta1.r;
  double r2=r21.norm2();
  double aux=G*m2*m1*pow(r2,-1.5);
  vector3D F1=r21*aux;
  Planeta1.SumeFuerza(F1);  Planeta2.SumeFuerza(F1*(-1));
}
//----------- Funciones Globales -----------

int main(){
  double r=1000.0; //Distancia entre Júpiter y el Sol.
  double m0= 1047.0;  //Masa del sol.
  double m1=1; //Masa de Júpiter.
  double M=m0+m1;
  double x0=-m1*r/M,x1=m0*r/M; //Condición inicial en la posición para movimiento circular.
  double omega=sqrt(G*M/(r*r*r)); //Condición inicial en la velocidad angular para movimiento circular.
  double T=2*M_PI/omega;
  double V0=omega*x0, V1=omega*x1; //Condiciones iniciales en la velocidad para movimiento circular.
  int numOrbitas = 20; //Número total de órbitas.
  double t, ttotal=numOrbitas*T;

  double dt=0.1; //Paso de tiempo
  Cuerpo Planeta[N];
  Interaction Newton;
  int i;

  //Iniciializacón.
  //---------------(x0,y0,z0,Vx0,Vy0,Vz0,m0)
  Planeta[0].Inicie(x0, 0, 0,  0, V0,  0,m0);
  Planeta[1].Inicie(x1, 0, 0,  0, V1,  0,m1);



  //Manejo de ceración de archivo.
  std::ofstream fout;
  fout.open("data.txt");


  //Se corren las iteraciones.
  for(t=0;t<ttotal;t+=dt){

    fout<<Planeta[0].Getx()<<" "<<Planeta[0].Gety()<<" "<<Planeta[1].Getx()<<" "<<Planeta[1].Gety()<<endl;

    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,xi);
    Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++) Planeta[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,Um2chiplusxi);
    Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++)Planeta[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Planeta); for(i=0;i<N;i++)Planeta[i].Mueva_V(dt,Um2lambdau2);
    for(i=0;i<N;i++) Planeta[i].Mueva_r(dt,xi);


  }

  fout.close();
  return 0;
}
