#include "functions.hpp"

int main(){

//-------------------Parametros de simulacion-------------------
//Parametros de algorimto de integracion
double xi=0.1786178958448091;
double lambda=-0.2123418310626054;
double chi=-0.06626458266981849;
double Um2lambdau2=(1-2*lambda)/2;
double Um2chiplusxi=1-2*(chi+xi);

//Parametros del problema fisico
int kT = 1;
double epsilon = 1.0;
double r0 = 10.0;
const int Nx=1; 
const int Ny=1; 
double Lx=60, Ly=60; 
double L = sqrt(Lx*Lx/ + Ly*Ly);
const int N = Nx*Ny;

//-------------------Condiciones iniciales-------------------
Cuerpo Moleculas[N]; //Moleculas[i] de ejemplo
Colisionador Newton; //Colisionador de ejemplo
Crandom ran64(1);
double x0=10.0, y0=0.0, m0=1.0, R0=3.0;
double v0=sqrt(kT/m0), theta, vx0, vy0;

// Paso de tiempo y tiempo maximo
int Ncuadros=1000; 
double t,tdibujo,dt=1e-3,tmax=L/v0,tcuadro=tmax/Ncuadros; 



for (int i = 0; i < N; i++) {
    theta=2*M_PI*ran64.r();
    vx0=v0*cos(theta); vy0=v0*sin(theta);
    Moleculas[i].Inicie(x0,y0,vx0,vy0,m0,R0);
    } //Atributos iniciales de la molécula

//-------------------Simulacion-------------------
//INICIO
// Inicio de la animacion
InicieAnimacion();

// Inicio de la simulacion
for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>=tcuadro){
      
      InicieCuadro();
      for(int i=0;i<N;i++) Moleculas[i].Dibujese();
      TermineCuadro();
      
      tdibujo=0;
    }
    
    
    for(int i = 0; i < N; i++) Moleculas[i].Mueva_r(dt,xi);
    Newton.CalculeTodasLasFuerzas(Moleculas, epsilon, r0, N);
    for(int i = 0; i < N; i++) Moleculas[i].Mueva_v(dt,Um2lambdau2);


    for(int i = 0; i < N; i++) Moleculas[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Moleculas, epsilon, r0, N);
    for(int i = 0; i < N; i++) Moleculas[i].Mueva_v(dt,lambda);

    for(int i = 0; i < N; i++) Moleculas[i].Mueva_r(dt,Um2chiplusxi);
    Newton.CalculeTodasLasFuerzas(Moleculas, epsilon, r0, N);
    for(int i = 0; i < N; i++) Moleculas[i].Mueva_v(dt,lambda);

    for (int i = 0; i < N; i++) Moleculas[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(Moleculas, epsilon, r0, N);
    for(int i = 0; i < N; i++) Moleculas[i].Mueva_v(dt,Um2lambdau2);

    for(int i = 0; i < N; i++) Moleculas[i].Mueva_r(dt,xi);

    }

  return 0;
}


//Difenrecias entre implementacion con punteros y referencias:
// -Sintactica:
// En la implementacion con punteros se tiene que acceder a los 
//atributos de la clase con el operador ->, mientras que en la 
//implementacion con referencias se accede con el operador .
//-Rendimiento:
// Al utilizar punteros en lugar de copias directas de datos, se 
// evita la sobrecarga asociada con la copia de objetos grandes, lo 
// que puede resultar en un rendimiento mejorado en términos de tiempo 
// de ejecución y consumo de memoria, especialmente si los objetos son
// grandes o si se están creando muchos de ellos.

// Sin embargo, también se deben considerar las implicaciones de la gestión 
// de la memoria al utilizar punteros. Si la clase Moleculas[i] no gestiona correctamente 
// la memoria (por ejemplo, no se eliminan adecuadamente los objetos creados dinámicamente), puede resultar en fugas de memoria y problemas de rendimiento.
