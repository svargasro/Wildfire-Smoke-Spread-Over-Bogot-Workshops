#include <iostream>
#include <cmath>
#include <fstream>
#include "vector.h"

using std::cout;
using std::sqrt;
using std::endl;
using std::fabs;
using std::sin;
using std::cos;

//Constantes del problema físico
const int N=3;
const double G=1.0;

//Constantes del algoritmo de integración
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

  //Variables tercer cuerpo.
  double x2 = x1*cos(M_PI/3.0); //Posición inicial en x
  double y2 = x1*sin(M_PI/3.0); //Posición inicial en y
  double V2x = -omega*y2; //Velocidad inicial en x

  V2x = V2x + 5*V2x/1000.0; //Perturbación de 5 partes por mil en la velocidad en x.

  double V2y = omega*x2; //Velocidad inicial en y

  double m2= 0.005; //Masa del tercer cuerpo

  int numOrbitas = 20; //Número total de órbitas.

  double t, ttotal=numOrbitas*T;

  double dt=1; //Paso de tiempo
  Cuerpo Planeta[N];
  Interaction Newton;
  int i;

  //INICIO
  //---------------(x0,y0,z0,Vx0,   Vy0,Vz0,m0)
  Planeta[0].Inicie(x0, 0, 0,  0, V0,  0,m0);
  Planeta[1].Inicie(x1, 0, 0,  0, V1,  0,m1);
  Planeta[2].Inicie(x2, y2, 0,  V2x, V2y,  0,m2);

  std::ofstream fout;
  fout.open("data.txt");

  //Se calcula el movimiento.
  for(t=0;t<ttotal;t+=dt){

    //Se declarar variables para las coordenadas del Sol y Júpiter en el sistema sin rotar puesto que serán utilizadas varias veces.
    double x_s = Planeta[0].Getx();
    double y_s = Planeta[0].Gety();
    double x_j = Planeta[1].Getx();
    double y_j = Planeta[1].Gety();

    double x_t = Planeta[2].Getx();
    double y_t = Planeta[2].Gety();


    double omegaT= omega*t; //Ángulo que se recorre por unidad de tiempo.

    /*Se cambia al sistema rotante (primado), que satisface el cambio de base dado por la matriz de rotación con theta, el ángulo de giro,
      igual a omega*T. */

    double x_sp = x_s*cos(omegaT) + y_s*sin(omegaT);
    double y_sp = -x_s*sin(omegaT) + y_s*cos(omegaT);

    double x_jp = x_j*cos(omegaT) + y_j*sin(omegaT);
    double y_jp = -x_j*sin(omegaT) + y_j*cos(omegaT);

    //Se agrega el cálculo de las coordenadas del planeta troyano en el sistema primado.
    double x_tp = x_t*cos(omegaT) + y_t*sin(omegaT);
    double y_tp = -x_t*sin(omegaT) + y_t*cos(omegaT);


    fout<<x_sp<<" "<<y_sp<<" "<<x_jp<<" "<<y_jp<<" "<<" "<<x_tp<<" "<<y_tp<<" "<<t<<endl; //Se pasan las coordenadas del sistema primado del Sol, Júpiter y del planeta troyano, además del tiempo.

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
