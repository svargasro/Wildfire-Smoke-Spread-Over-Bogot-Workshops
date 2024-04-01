#include <vector>
#include <cmath> 
#include <utility> 
#include <tuple>
#include <iostream>

//Constantes del problema físico
const int KbT = 1;
const double epsilon = 1.0;
const double r0 = 10.0;

//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);

class Molecula;

//Función que calcula las componentes de la fuerza a partir de la posición
std::pair<double, double> CalculeFuerza(double x, double y){
  double r = std::sqrt(x*x + y*y);
  double Fr = ((12*epsilon)/r)*(std::pow(r0/r,12) - std::pow(r0/r,6));
  double theta = std::atan2(y,x);
  double Fx = Fr*std::cos(theta);
  double Fy = Fr*std::sin(theta);
  return std::make_pair(Fx,Fy); //Crea un contenedor con las componentes de la fuerza
};

class Molecula{
private:
    double x, y, Vx, Vy, Fx, Fy;
    double m, R;
public:
    void Inicie(double x0,double y0,double Vx0,double Vy0,double m,double R);
    void Mueva_r(double dt,double coeficiente);
    void Mueva_V(double dt,double coeficiente, double Fx, double Fy);
    double Getx(void){return x;}; // Inline
    double Gety(void){return y;}; // Inline
};

void Molecula::Inicie(double x0,double y0,double Vx0,double Vy0,double m,double R){
  x=x0; y=y0; Vx=Vx0; Vy=Vy0; Fx=0; Fy=0; this->m=m; this->R=R;
}

void Molecula::Mueva_r(double dt,double coeficiente){
  x+=Vx*(coeficiente*dt);
  y+=Vy*(coeficiente*dt);
}
void Molecula::Mueva_V(double dt,double coeficiente, double Fx, double Fy){
  Vx+=Fx*(coeficiente*dt);
  Vy+=Fy*(coeficiente*dt);
}

int main (){
    double dt=1e-2; //dt ajustado para que la energía se conserve
    int t=100; //t max
    double Fx, Fy;
    Molecula Pepe; //Molecula de ejemplo
    double x0=10.0, y0=0.0, m=1.0, R=3.0;
    double Vx0=sqrt(KbT/m), Vy0=0.0; 
    Pepe.Inicie(x0,y0,Vx0,Vy0,m,R); //Atributos iniciales de la molécula

    for (double i = 0; i < t; i=i+dt){
      //Implementación Algoritmo Omelyan Position-Extended-Forest-Ruth-Like (PEFRL) 4th Order (2002)
      Pepe.Mueva_r(dt,xi);
      std::tie(Fx,Fy)=CalculeFuerza(Pepe.Getx(), Pepe.Gety()); //std::tie desempaqueta el resultado de la función CalculeFuerza en dos variable  
      Pepe.Mueva_V(dt,Um2lambdau2,Fx,Fy);

      Pepe.Mueva_r(dt,chi);
      std::tie(Fx,Fy)=CalculeFuerza(Pepe.Getx(), Pepe.Gety());
      Pepe .Mueva_V(dt,lambda,Fx,Fy);

      Pepe.Mueva_r(dt,Um2chiplusxi);    
      std::tie(Fx,Fy)=CalculeFuerza(Pepe.Getx(), Pepe.Gety());
      Pepe.Mueva_V(dt,lambda,Fx,Fy);

      Pepe.Mueva_r(dt,chi);
      std::tie(Fx,Fy)=CalculeFuerza(Pepe.Getx(), Pepe.Gety());
      Pepe.Mueva_V(dt,Um2lambdau2,Fx,Fy);
      
      Pepe.Mueva_r(dt,xi);  

      std::cout<<i<<"\t"<<Pepe.Getx()<<std::endl;
    }
    return 0;
}