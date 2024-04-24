#include <iostream>
#include <chrono>
#include <ctime>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include "../vector.h"

using namespace std;

//--------------- Constantes globales ------------
const double g=9.8, L=12;


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
  void Mueva_theta(double dt,double coeficiente);
  void Mueva_omega(double dt,double coeficiente);
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
  void CalculeTodosLosTorques(Cuerpo * pendulos, double K);
  void CalculeTorqueEntre(Cuerpo & moleculas1,Cuerpo & moleculas2, double K);
};


//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0, double y0, double theta0, double omega0, 
            double m0, double R0){
  xs=x0; ys=y0;
  theta=theta0; omega=omega0; 
  m=m0; R=R0;
}

void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(coeficiente*dt);
  
}

void Cuerpo::Mueva_omega(double dt,double coeficiente){
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
void Colisionador::CalculeTodosLosTorques(Cuerpo * pendulos, double K){
  int i, j;

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
  		CalculeTorqueEntre(pendulos[i],pendulos[j], K);
  }
}

void Colisionador::CalculeTorqueEntre(Cuerpo & pendulos1, Cuerpo & pendulos2, double K){
    
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
void InicieAnimacion(double K){
  double base = K / 1e10;

  ostringstream oss;
  oss << std::scientific << std::setprecision(1) << base << "e+10";

  clog<<"'Pendulos_K="<<oss.str()<<".gif'"<<endl;

  cout<<"set terminal gif animate delay 10"<<endl; 
  cout<<"set output 'Pendulos_K="<<oss.str()<<".gif'"<<endl;
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
}

void TermineCuadro(void){
    cout<<endl;
    cout<<"unset label"<<endl;
}

void Datos(Cuerpo & pendulo, double K){ 

  double base = K / 1e10;

  ostringstream oss;
  oss << std::scientific << std::setprecision(1) << base << "e+10";
  string filename = "PenduloXvsY_K=" + oss.str()  + ".txt";

  ofstream archivo(filename, ios::app); // Open the file in append mode
  archivo << pendulo.GetX() << " " << pendulo.GetY() << endl; // Write the value of tau to the file
  archivo.close(); // Close the file
}

void graficar(double K){

  double base = K / 1e10;

  ostringstream oss;
  oss << std::scientific << std::setprecision(1) << base << "e+10";

  ofstream gp("graph.gp", ios::out);
  if (!gp) {
    cout << "Error al abrir el archivo" << endl;
    return;
  }
  gp << "n=" << N << endl;
  gp << "set terminal png size 800,600" << endl;
  gp << "set output 'pendulos.png'" << endl;
  gp << "set title 'Trayectoria de los péndulos'" << endl;
  gp << "set xlabel 'Posición X'" << endl;
  gp << "set ylabel 'Posición Y'" << endl;
  gp << "plot for [i=0:n-1] 'PenduloXvsY_K="<<oss.str()<<".txt' every n::i using 1:2 with linespoints pt 7 ps 1 title sprintf('Péndulo %d', i+1)" << endl;
  gp << "set output" << endl;
  gp.close();
}

//------------------ Funciones rendimiento ------------------
void time(double & wsum , double & wsum2 , double & csum , double & csum2,
          double & wtime , double & ctime )
{
    
    wsum  += wtime;
    wsum2 += wtime*wtime;
    csum  += ctime;
    csum2 += ctime*ctime;
    
    
}

void stats(double & wsum , double & wsum2 , double & csum , double & csum2,
          double & wtime , double & ctime, double reps,
          double & mean_wtime, double & sigma_wtime,
          double & mean_ctime, double & sigma_ctime, double K)
{
  mean_wtime = wsum/reps;
  sigma_wtime = std::sqrt(reps*(wsum2/reps - mean_wtime*mean_wtime)/(reps-1));
  mean_ctime = csum/reps;
  sigma_ctime = std::sqrt(reps*(csum2/reps - mean_ctime*mean_ctime)/(reps-1));
  
  double base = K / 1e10;

  ostringstream oss;
  oss << std::scientific << std::setprecision(1) << base << "e+10";
  string filename = "Benchmark_CunaNewton_K=" + oss.str() + ".txt";
  
  ofstream archivo(filename, ios::app); // Open the file in append mode
  archivo<<"\n- Simulación con K = "<<K<<endl;
  archivo<<"\nmean_wtime = "<<mean_wtime<<"s  sigma_wtime = "<<sigma_wtime<<"s"<<endl;
  archivo<<"mean_ctime = "<<mean_ctime<<"s sigma_ctime = "<<sigma_ctime<<"s"<<endl;

  archivo.close();
}

int main(int argc, char *argv[]){
  
  double K = atof(argv[1])*1e10;
  clog<<"\n\n\nK = "<<K<<"\n\n\n";

  Cuerpo pendulos[N];
  Colisionador Newton;

  //Parametros de la simulación
  double m0=0.1; double R0=1.5;
  double theta01=-M_PI/12;
  double theta0i=0;
  double T = M_PI*(L/g);
  
  double x0=0,y0=L,omega0=0;
  double dx = 2*R0;//Separación entre pendulos 
 
  //Variables auxiliares para correr la simulacion
  int i, Ncuadros=100; 
  double t,tdibujo,dt=1e-5,tmax=5*T,tcuadro=tmax/Ncuadros; 
  clog<<"T = "<<T<<"s tmax = "<<tmax<<"s"<<endl;

  
  //Variables para medir el rendimiento
  double mean_wtime, sigma_wtime;
  double mean_ctime, sigma_ctime;
  double wsum = 0, wsum2 = 0, csum = 0, csum2 = 0;
  double wtime, ctime;
  double reps = tmax/dt;
  auto start = std::chrono::system_clock::now(); // measures wall time
  std::clock_t c1 = std::clock();

  InicieAnimacion(K);
  graficar(K);
  
  //Inicie l0s pendulos
  //Inicializacion del primer pendulo
  pendulos[0].Inicie(x0,y0,theta01,omega0,m0,R0);
  
  //Inicializacion de los demas pendulos
  for(i=1;i<N;i++){
    x0=i*dx;
    pendulos[i].Inicie(x0,y0,theta0i,omega0,m0,R0);
  }
  
 
  //CORRO
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    // if(tdibujo>=tcuadro){
      
    //   InicieCuadro(t);
    //   for(i=0;i<N;i++) {
    //     pendulos[i].Dibujese();
    //     Datos(pendulos[i], K); 
    //     }
    //   TermineCuadro();
      
    //   tdibujo=0;
    // }
  
    
    // for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,xi);    
    // Newton.CalculeTodosLosTorques(pendulos,K); 
    // for(i=0;i<N;i++) pendulos[i].Mueva_omega(dt,Um2lambdau2);

    // for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,chi);
    // Newton.CalculeTodosLosTorques(pendulos,K); 
    // for(i=0;i<N;i++) pendulos[i].Mueva_omega(dt,lambda);
    
    // for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,Um2chiplusxi);
    // Newton.CalculeTodosLosTorques(pendulos,K); 
    // for(i=0;i<N;i++)pendulos[i].Mueva_omega(dt,lambda);

    // for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,chi);
    // Newton.CalculeTodosLosTorques(pendulos,K); 
    // for(i=0;i<N;i++)pendulos[i].Mueva_omega(dt,Um2lambdau2);

    // for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,xi);

    clog<<"\nPorcentaje de avance: "<<(t/tmax)*100<<"%"<<endl;

    auto end = std::chrono::system_clock::now(); // wall time
    std::clock_t c2 = std::clock(); // cpu time

    std::chrono::duration<double> elapsed_seconds = end-start;
    
    ctime = 1.0*(c2-c1)/CLOCKS_PER_SEC;
    wtime = elapsed_seconds.count();

    time(wsum, wsum2, csum, csum2, wtime, ctime);

  }

  stats(wsum, wsum2, csum, csum2, wtime, ctime, reps, mean_wtime, sigma_wtime, mean_ctime, sigma_ctime, K);

  return 0;
}
