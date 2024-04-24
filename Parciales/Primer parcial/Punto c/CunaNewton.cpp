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

//Constantes globales

const int N=3;
//
const double g=980.665, L=12;
//Gravedad (cm/s^2), Longitud del pendulo (cm)

//Dimensiones de la caja
const double Lx=0, Ly=0; 

//Constantes del algoritmo de integración
const double xi=0.1786178958448091;
const double lambda=-0.2123418310626054;
const double chi=-0.06626458266981849;
const double Um2lambdau2=(1-2*lambda)/2;
const double Um2chiplusxi=1-2*(chi+xi);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  double xs,ys;//punto de apoyo del pendulo
  double theta, omega, tau;
  double m,R, I;
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
  double GetTau(void){return tau;};
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
  I = m0*L*L;
}

void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(coeficiente*dt);
  
}

void Cuerpo::Mueva_omega(double dt,double coeficiente){
  
  omega += (tau / I) * (coeficiente * dt);
  
}
void Cuerpo::Dibujese(void){

  cout<<" , "<<GetX()<<"+"<<R<<"*cos(t),"<<GetY()<<"+"<<R<<"*sin(t) , ";
  cout<<xs<<"+"<<(GetX() - xs)/7.0<<"*t,"<<ys<<"+"<<(GetY()- ys)/7.0<<"*t";
  
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
void InicieAnimacion(string K){
  
  cout<<"set terminal gif animate delay 10"<<endl; 
  cout<<"set output 'Pendulos_K="<<K<<".gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-20:"<<Lx + 20<<"]"<<endl;
  cout<<"set yrange[-20:"<<Ly + 20<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl; 
  cout<<"set title 'Cuna de Newton para K="<<K<<"'"<<endl; 
}

void InicieCuadro(double t){
  cout <<"set label \'t=" << t << "s \' at graph 0.5,0.9 center" << endl;
  cout<<"plot 0,0 ";
}

void TermineCuadro(void){
    cout<<endl;
    cout<<"unset label"<<endl;
}

void Datos(Cuerpo & pendulo, double K, ofstream & archivo){   
  archivo << pendulo.GetX() << " " << pendulo.GetY() << endl; // Write the value of tau to the file
}

void graficar(string K, ofstream & gp){
  
  if (!gp) {
    cout << "Error al abrir el archivo" << endl;
    return;
  }
  gp << "n=" << N << endl;
  gp << "set terminal png size 800,600" << endl;
  gp << "set output 'Pendulos_K="<<K<<".png'" << endl;
  gp << "set title 'Trayectoria de los péndulos'" << endl;
  gp << "set xlabel 'Posición X'" << endl;
  gp << "set ylabel 'Posición Y'" << endl;
  gp << "plot for [i=0:n-1] 'PenduloXvsY_K="<<K<<".txt' every n::i using 1:2 with linespoints pt 7 ps 1 title sprintf('Péndulo %d', i+1)" << endl;
  gp << "set output" << endl;
  gp.close();
}


//-----------------Funciones numeral c -----------------------
void DatosTorque(double Tau, double t, ofstream & archivoTauK) {
  archivoTauK<<t<<" "<<Tau <<endl;
}

// void graficarTau(string K, ofstream & gp, bool isFirstTime ){
  
//   if (isFirstTime) {
//     gp << "set terminal png size 800,600 enhanced font 'Verdana,10'" << endl;
//     gp << "set output 'CunaNewtonTau_AllK.png'" << endl;
//     gp << "set title 'Torque vs tiempo para varios K'" << endl;
//     gp << "set xlabel 'Tiempo'" << endl;
//     gp << "set ylabel 'Torque'" << endl;
//     gp << "set grid" << endl;
//     gp <<"style data points"<<endl;
//     gp << "plot 'CunaNewtonTau_K="<<K<<".txt' using 1:2 with linespoints pt 7 ps 1 title 'Torque vs tiempo'" << endl;
//     isFirstTime = false;
//   }
//   if (!gp) {
//     cout << "Error al abrir el archivo" << endl;
//     return;
//   }
//   else{
//     gp << "plot 'CunaNewtonTau_K="<<K<<".txt' using 1:2 with linespoints pt 7 ps 1 title 'Torque vs tiempo'" << endl;
//   }
// }
// Error en la logica de programación, no se puede hacer un grafico por cada iteración del bucle, 
// no se puede tener información sobre cuantas veces entra debido a que al llamarse el archivo
// parallelmente se pierde la información de la variable isFirstTime, por lo que se debe hacer
// un grafico por cada K, y no por cada iteración del bucle, pero esto no es lo que se quiere
// por tanto se va a crear un archivo gp que ya contenga las intrucciiones para graficar todos los K

//-----------------Funciones numeral d -----------------------
void DatosTorqueMax(double tMax, double TauMax, double K, ofstream & archivo){
   // Open the file in append mode
  archivo << tMax << " " << TauMax << endl; // Write the data to the file
}

void TorqueMax(double Tau, double & TauMax){

  if(Tau > TauMax){
    TauMax = Tau;
  }
  
}
void TimeOscilationMax(double Tau, double t, double & tMax, bool & firstTime){
  
  double t_start;
  double t_end;
  double pre = 1e6;
  // clog<<Tau<<endl;
  if(abs(Tau) > pre){
      if (firstTime) {
        t_start = t;
        firstTime = false;
        }
      else
        t_end = t;
  }
  tMax = t_end - t_start;

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


  string filename = "Benchmark_CunaNewton.txt";
  ofstream archivo(filename, ios::app); // Open the file in append mode
  archivo<<"\n- Simulación con K = "<<K<<endl;
  archivo<<"\nmean_wtime = "<<mean_wtime<<"s  sigma_wtime = "<<sigma_wtime<<"s"<<endl;
  archivo<<"mean_ctime = "<<mean_ctime<<"s sigma_ctime = "<<sigma_ctime<<"s"<<endl;

  archivo.close();
}

 
  
  
int main(int argc, char *argv[]){
  
  float K = atof(argv[1])*1e9;
  clog<<"\n\n\nK = "<<K<<"\n\n\n";

  Cuerpo pendulos[N];
  Colisionador Newton;
  
  //Archivos para guardar los datos
  ostringstream ossk;
  ossk << scientific << setprecision(1) << static_cast<double>(K);
 
  ofstream archivoXvsY("PenduloXvsY_K=" + ossk.str() + ".txt", ios::app); //Archivo para guardar las posiciones de los pendulos
  ofstream archivoTauK("CunaNewtonTau_K=" + ossk.str() + ".txt", ios::app); // Archivo para guardar el torque de un pendulo en funcion del tiempo
  ofstream gp_pendulos("PenduloXvsY_K=" + ossk.str() + "_graph.gp", ios::app); // Archivo para graficar las posiciones de los pendulos
  // ofstream gp_tau("CunaNewtonTau_AllK.gp", ios::app); // Archivo para graficar el torque de un pendulo en funcion del tiempo
  ofstream archivotMaxTauMax("Punto(tMax,TauMax)_K=" + ossk.str() + ".txt", ios::app); // Archivo para guardar el punto de torque maximo y tiempo
  ofstream benchmark("Benchmark_CunaNewton.txt", ios::app); // Open the file in append mode

 
  //Variables auxiliares para correr la simulacion
  double m0=100, R0=1.5;
  double theta01=-M_PI/12;
  double theta0i=0;
  double T=2*M_PI*sqrt(L/g);
  double ajuste = 0.2; //Como calcular este ajuste?
  double T_ajuste = ajuste*T;
  double T_total = T + T_ajuste;
  // double t,tstart = 0.174,tmax=0.184,dt=1e-5;
  int i, Ncuadros=100;
  double t,tstart = 0.174,tdibujo,dt=1e-5,tmax=2*T_total ,tcuadro=tmax/Ncuadros; 

  double x0=0,y0=L,omega0=0;
  double dx = 2*R0;//Separación entre pendulos 

  double tMax = 0, TauMax = 0;

  
  //Variables para medir el rendimiento
  double mean_wtime, sigma_wtime;
  double mean_ctime, sigma_ctime;
  double wsum = 0, wsum2 = 0, csum = 0, csum2 = 0;
  double wtime, ctime;
  double reps = tmax/dt;
  auto start = std::chrono::system_clock::now(); // measures wall time
  std::clock_t c1 = std::clock();
  bool firstTime = true;
      
  
  InicieAnimacion(ossk.str());
  //No usar salida estandar para los datos de la animación

  graficar(ossk.str(), gp_pendulos);
  // graficarTau(ossk.str(), gp_tau);

  
  //Inicie l0s pendulos
  //Inicializacion del primer pendulo
  pendulos[0].Inicie(x0,y0,theta01,omega0,m0,R0);
  
  //Inicializacion de los demas pendulos
  for(i=1;i<N;i++){
    x0=i*dx;
    pendulos[i].Inicie(x0,y0,theta0i,omega0,m0,R0);
  }
  
  //Inicie la animación
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>=tcuadro){
      
      InicieCuadro(t);
      for(i=0;i<N;i++) {
        pendulos[i].Dibujese();
        
        }
      TermineCuadro();
      
      tdibujo=0;
    }
    for (i = 0; i < N; i++) {
      Datos(pendulos[i], K, archivoXvsY);
    }
    if(t > tstart && t<tstart + 0.01) { 
          double Tau = pendulos[1].GetTau();
          
          TorqueMax(Tau, TauMax);
          TimeOscilationMax(Tau, t, tMax, firstTime);
          DatosTorque(Tau,t,archivoTauK);
          //  archivoTauK<<t<<" "<<pendulos[1].GetTau() <<"\n";
        }
    



    // archivoTauK<<t<<" "<<pendulos[1].GetTau() <<"\n";
    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,xi);    
    Newton.CalculeTodosLosTorques(pendulos,K); 
    for(i=0;i<N;i++) pendulos[i].Mueva_omega(dt,Um2lambdau2);

    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,chi);
    Newton.CalculeTodosLosTorques(pendulos,K); 
    for(i=0;i<N;i++) pendulos[i].Mueva_omega(dt,lambda);
    
    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,Um2chiplusxi);
    Newton.CalculeTodosLosTorques(pendulos,K); 
    for(i=0;i<N;i++)pendulos[i].Mueva_omega(dt,lambda);

    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,chi);
    Newton.CalculeTodosLosTorques(pendulos,K); 
    for(i=0;i<N;i++)pendulos[i].Mueva_omega(dt,Um2lambdau2);

    for(i=0;i<N;i++) pendulos[i].Mueva_theta(dt,xi);
    

      // DatosTorque(pendulos[1],t,K,TauPrevio,cambioDeSigno,TauActual,tMax,TauMax,archivoTauK);
      // DatosXyY(pendulos[1], K, archivoXvsY);
  
      clog<<"\nPorcentaje de avance: "<<(t/tmax)*100<<"%"<<endl;

      auto end = std::chrono::system_clock::now(); // wall time
      std::clock_t c2 = std::clock(); // cpu time

      std::chrono::duration<double> elapsed_seconds = end-start;
      
      ctime = 1.0*(c2-c1)/CLOCKS_PER_SEC;
      wtime = elapsed_seconds.count();

      time(wsum, wsum2, csum, csum2, wtime, ctime);

    }
    clog<<"Tiempo de simulación: "<<wtime<<"s"<<endl;
    clog<<tMax<<endl;
    DatosTorqueMax(tMax, TauMax, K, archivotMaxTauMax);
    stats(wsum, wsum2, csum, csum2, wtime, ctime, reps, mean_wtime, sigma_wtime, mean_ctime, sigma_ctime, K);
    
    //Cierre de archivos
    archivoXvsY.close();
    archivotMaxTauMax.close();
    archivoTauK.close();
    gp_pendulos.close();
    // gp_tau.close();

    
  
  
  return 0;
}
