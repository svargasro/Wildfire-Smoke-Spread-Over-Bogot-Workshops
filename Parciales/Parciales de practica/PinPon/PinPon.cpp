#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "../vector.h"
#include "../Random64.h"
using namespace std;

//Constantes del problema físico
const double g=9.81;
//Longitud de la raqueta
const double L=10;

//Número de moleculass
const double Lx=12, Ly=50; 
const int N = 1;

//Elasticidad de la colisión y 
const double K=1.0e4;


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
//Clase Cuerpo
class Cuerpo{
private:
  vector3D r,V,F; double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};// Inline
  void SumeFuerza(vector3D dF){F+=dF;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void Dibujese(void);
  double GetX(void){return r.x();}
  double GetY(void){return r.y();}
  friend class Colisionador;
};

//Clase Raqueta
class Raqueta{
private:
  vector3D r,V;
  //Parámetros de la raqueta
  double x0,y0,z0,Vx0,Vy0,Vz0;
  //Amplitud, frecuencia angular, longitud
  double A,omega;
  friend class Colisionador;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
              double omega0, double A0);
  void Muevase(double dt);
  void Dibujese(void);
};

//Clase Colisionador
class Colisionador{
private:
public:
  void CalculeTodasLasFuerzas(Cuerpo * moleculas, Raqueta & raqueta);
  void CalculeFuerzaEntre(Cuerpo & moleculas, Raqueta & raqueta);
};



//-------Implementar las funciones de las clases------

//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.load(x0,y0,z0);  V.load(Vx0,Vy0,Vz0); m=m0; R=R0;
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  r+=V*(coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt,double coeficiente){
  V+=F*(coeficiente*dt/m);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//------- Funciones de la clase Raqueta --------
void Raqueta::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,
              double omega0, double A0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0); 
  omega=omega0; A=A0;
}
void Raqueta::Muevase(double t){
  r.load(0,A*sin(omega*t),0);

  // V.load(0,-A*omega*cos(omega*dt),0);
  // r+=V*dt;
}
void Raqueta::Dibujese(void){
  cout << ", " << r.x() + L<< "*cos(t)," << r.y();
}

//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * moleculas, Raqueta & raqueta){
  int i,j;
  //Borro las fuerzas de todos los moleculass  
  for(i=0;i<N;i++)
    moleculas[i].BorreFuerza();

  //Fuerza de gravedad
  vector3D Fg;
  
  for(i=0;i<N;i++){
    Fg.load(0,-moleculas[i].m*g,0);
    moleculas[i].SumeFuerza(Fg);
  }
  
  //Recorro por parejas y calculo la fuerza de colision
  for(i=0;i<N;i++) CalculeFuerzaEntre(moleculas[i], raqueta);
     
}

void Colisionador::CalculeFuerzaEntre(Cuerpo & moleculas1, Raqueta & raqueta){
  
  //Calcular el vector normal
  vector3D r21=raqueta.r-moleculas1.r; 
  double d=r21.norm();


  //Determinar si hay colisión
  double s= moleculas1.R  - d;
  
  if(s>0){
    //Fuerza de la colisión
    double m = moleculas1.m;
    double Vy = moleculas1.V.y() - raqueta.V.y();
    double F= K*pow(s,1.5) - gamma_coef*m*Vy*sqrt(s); 

    if (F<0) F=0; //Evitar la interpenetración (fuerza negativa)
    else {

      //Calcular vector normal
    vector3D n = r21*(1.0/d);

    //Calcular la fuerza
    vector3D F1=n*F;

    //Sumar las fuerzas
    moleculas1.SumeFuerza(F1*(-1));

    }
  }
}


//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'PinPon.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-"<<Lx <<":"<<Lx <<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly + 10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}

void TermineCuadro(void){
    cout<<endl;
}

//---Funciones del taller---


ofstream puntoA(Cuerpo & molecula,double dt){
  string filename = "PinPonA" + to_string(dt) + ".txt";
  ofstream archivo(filename, ios::app); // Open the file in append mode
  numeral1();
  numeral2();
  return archivo;
}
void numeral1(void){
  gamma_coef = 0;
}
void numeral2(void){
  gamma_coef = 10;
}
ofstream puntoB(double t, Cuerpo & molecula){
  string filename = "PinPonB.txt";
  ofstream archivo(filename, ios::app); // Open the file in append mode

  A=1; T=5; tmax = 200;
}
// void puntoC(void){
//   string filename = "PinPonB.txt";
//   ofstream archivo(filename, ios::app); // Open the file in append mode

//   double t,tdibujo,dt=1e-3,tmax=100,tcuadro=tmax/Ncuadros; 
  
//   double T_inicial = 0.1, T_final = 10, incremento_T = 0.1;
//   A=1;
//   double umbral_energia = 0.1;
//   //Variar el periodo para encontrar T_critico
//   double T_c = T_critico(A, T_inicial, T_final, incremento_T, umbral_energia, tmax, dt); 

// }

// //Funcion para encontrar T_critico
// double T_critico(double A, double T_inicial, double T_final, double incremento_T, double umbral_energia, 
//                 double tiempo_maximo, double delta_tiempo){
//   for (double T_candidato = T_inicial; T_candidato <= T_final; T_candidato += incremento_T) {
//         if (correr_simulacion(T_candidato, A, umbral_energia, tiempo_maximo, delta_tiempo)) {
//             return T_candidato;  // Encontramos el periodo crítico T_c
//         }
//     }
//   }

// // Función para correr la simulación
// bool correr_simulacion(double T, double A, double umbral_energia, double tiempo_maximo, double delta_tiempo) {
//     double energia_inicial = calcular_energia(A);
//     double tiempo = 0;
//     double energia_actual = energia_inicial;

//     //Parametros de la simulación
//     double m0=1; double R0=2;

//     //Variables auxiliares para la condición inicial
//     //Condiciones iniciales cuerpo
//     double x0 = 0,y0 = 30,Vx0 = 0,Vy0 = 0;
//     //Condiciones iniciales raqueta
//     double x0r = 0,y0r = 0,Vx0r = 0,Vy0r = 0;
//     double omega0 = 2*M_PI/T;


//     while (tiempo < tiempo_maximo) {

//         // Correr la simulación
//         simulacion(x0, y0, Vx0, Vy0, x0r, y0r, Vx0r, Vy0r, omega0, m0, R0);
        
//         //Mirar como acceder a los objetos con los que se va a trabajar
//         //Parar mirar que criterio establecer para saber si el sistema es 
//         //caotico

//         energia_actual = calcular_energia_actual();

//         if (es_caotico(energia_inicial, energia_actual, umbral_energia)) {
//             return true;  // El sistema es caótico para este periodo T
//         }

//         tiempo += delta_tiempo;
//     }

//     return false;  // El sistema no fue caótico
// }

void graficar(Cuerpo * molecula, Raqueta & raqueta){
 if(tdibujo>=tcuadro){
      
      InicieCuadro();
      for(i=0;i<N;i++) molecula[i].Dibujese();
      raqueta.Dibujese();
      TermineCuadro();
    }
}

void paso(double dt, Cuerpo * moleculas, Raqueta & raqueta, Colisionador & Newton){
  int i;    
  for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,xi);
    raqueta.Muevase(t);     
    Newton.CalculeTodasLasFuerzas(moleculas, raqueta); 
    for(i=0;i<N;i++) moleculas[i].Mueva_V(dt,Um2lambdau2);

    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,chi);
    raqueta.Muevase(t); 
    Newton.CalculeTodasLasFuerzas(moleculas, raqueta);
    for(i=0;i<N;i++) moleculas[i].Mueva_V(dt,lambda);
    
    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,Um2chiplusxi);
    raqueta.Muevase(t); 
    Newton.CalculeTodasLasFuerzas(moleculas, raqueta);
    for(i=0;i<N;i++)moleculas[i].Mueva_V(dt,lambda);

    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,chi);
    raqueta.Muevase(t); 
    Newton.CalculeTodasLasFuerzas(moleculas, raqueta);
    for(i=0;i<N;i++)moleculas[i].Mueva_V(dt,Um2lambdau2);

    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,xi);
    raqueta.Muevase(t); 

}

double simulacion(Cuerpo * moleculas, Raqueta & raqueta, Colisionador & Newton){

  int Ncuadros=100, i; 
  double t,tdibujo,dt=1e-3,tmax=14,tcuadro=tmax/Ncuadros; 

  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    
    graficar(moleculas, raqueta);
    
    paso(dt, moleculas, raqueta, Newton);
    
  }
}

// // Función para calcular la energía del sistema
// double calcular_energia(double A) {
//     //Energía de la partícula
//     double E = 0.5 * m * V.norm2();
//     return E;
// }

// // Función para mirar si el sistema es caótico
// bool es_caotico(double energia_inicial, double energia_actual, double umbral_energia) {
//     // Un sistema se considera caótico si hay una gran variación en la energía
//     // Esta es una simplificación y puede necesitar ajuste según tu sistema específico
//     double variacion_energia = std::abs(energia_actual - energia_inicial);
//     return variacion_energia > umbral_energia;
// }

// //Fucniiones para guardar datos
// void datos_yparticula(double t, Cuerpo & molecula, ofstream &archivo){
//   archivo << t<< "\t"<< molecula.GetY() << endl; // Write the value of tau to the file
// }



//A=1, Tc=2

double T_critico(Cuerpo * moleculas, Raqueta & raqueta, Colisionador & Newton){
  double dt = 1e-3;
  double previus_position = moleculas[0].GetY();
  double current_position;

  for (int i = 0; i < 100; i++) {
    paso(dt, moleculas, raqueta, Newton);
    current_position = moleculas[0].GetY();
    if (current_position < previus_position) {
      clog<<"T_critico: "<<T<<endl;
    }
    previus_position = current_position;
  }
}

double gamma_coef = 10; 
double T = 10, A = 1;


int main(int argc, char ** argv){
// int main(){
  
  A = atof(argv[1]);
  T = atof(argv[2]); 
  

  //Parametros de la simulación
  double m0=1; double R0=2;

  //Variables auxiliares para la condición inicial
  //Condiciones iniciales cuerpo
  double x0 = 0,y0 = 30,Vx0 = 0,Vy0 = 0;
  //Condiciones iniciales raqueta
  double x0r = 0, y0r = 0, Vx0r = 0, Vy0r = 0;
  double omega0 = 2*M_PI/T;

  

  InicieAnimacion();

  Cuerpo moleculas[N];
  Raqueta raqueta;
  Colisionador Newton;

  int i;
  
  for (i=0;i<N;i++) moleculas[i].Inicie(x0,y0, 0,Vx0,Vy0,  0,m0,R0);

  //Inicializar de la raqueta
  //------------(x0,y0,z0,Vx0,Vy0,Vz0,omega0)
  raqueta.Inicie(x0r,y0r,0,Vx0r,Vy0r,0,omega0,A);

  simulacion(moleculas, raqueta, Newton);

  //INICIO
  // for (i=0;i<N;i++) moleculas[i].Inicie(x0,y0, 0,Vx0,Vy0,  0,m0,R0);

  // //Inicializar de la raqueta
  // //------------(x0,y0,z0,Vx0,Vy0,Vz0,omega0)
  // raqueta.Inicie(x0r,y0r,0,Vx0r,Vy0r,0,omega0,A);
 
  //CORRO
  // for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

  //   if(tdibujo>=tcuadro){
      
  //     InicieCuadro();
  //     for(i=0;i<N;i++) moleculas[i].Dibujese();
  //     raqueta.Dibujese();
  //     TermineCuadro();
      
  //     // puntoA();
  //     // puntoB();
  //     // puntoC();

  //     tdibujo=0;
  //   }
    
    
  //   for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,xi);
  //   raqueta.Muevase(t);     
  //   Newton.CalculeTodasLasFuerzas(moleculas, raqueta); 
  //   for(i=0;i<N;i++) moleculas[i].Mueva_V(dt,Um2lambdau2);

  //   for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,chi);
  //   raqueta.Muevase(t); 
  //   Newton.CalculeTodasLasFuerzas(moleculas, raqueta);
  //   for(i=0;i<N;i++) moleculas[i].Mueva_V(dt,lambda);
    
  //   for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,Um2chiplusxi);
  //   raqueta.Muevase(t); 
  //   Newton.CalculeTodasLasFuerzas(moleculas, raqueta);
  //   for(i=0;i<N;i++)moleculas[i].Mueva_V(dt,lambda);

  //   for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,chi);
  //   raqueta.Muevase(t); 
  //   Newton.CalculeTodasLasFuerzas(moleculas, raqueta);
  //   for(i=0;i<N;i++)moleculas[i].Mueva_V(dt,Um2lambdau2);

  //   for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,xi);
  //   raqueta.Muevase(t); 
    
  // }
  return 0;

}
