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
const double g=980, L=12;

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
void InicieAnimacion(double K, double dt){
  
  ostringstream oss, osst;
  oss << scientific << setprecision(1) << static_cast<double>(K);
  osst << scientific << setprecision(1) << static_cast<double>(dt);

  cout<<"set terminal gif animate delay 1"<<endl; 
  cout<<"set output 'Pendulos_K="<<oss.str()<<".gif'"<<endl;
  cout<<"set title 'Pendulos con K=" << oss.str() <<" kg/(m1/2s2) y dt=" << oss.str() <<" s'" << endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-20:"<<Lx + 20<<"]"<<endl;
  cout<<"set yrange[-20:"<<Ly + 20<<"]"<<endl;
  cout << "set xlabel 'X(cm)'" << endl;
  cout << "set ylabel 'Y(cm)'" << endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}

void InicieCuadro(double t){
  cout << "set label \'t=" << t*1e2 << "s \' at graph 0.5,0.9 center" << endl;
  cout<<"plot 0,0 ";
}

void TermineCuadro(void){
    cout<<endl;
    cout<<"unset label"<<endl;
}

//-----------------Funcion numeral a -----------------------
void Datos(Cuerpo & pendulo, double K, ofstream & archivo){

  archivo << pendulo.GetX() << " " << pendulo.GetY() << endl; // Write the value of Tau to the file
   // Close the file
}

void graficar(double K){

  ostringstream oss;
  oss << std::scientific << std::setprecision(1) << static_cast<double>(K);

  ofstream gp("graph.gp", ios::out);
  if (!gp) {
    cout << "Error al abrir el archivo" << endl;
    return;
  }
  gp << "n=" << N << endl;
  gp << "set terminal png size 800,600" << endl;
  gp << "set output 'Pendulos_K="<<oss.str()<<".png'" << endl;
  gp << "set title 'Trayectoria de los péndulos_K="<<oss.str()<<"'" << endl;
  gp << "set xlabel 'Posición X'" << endl;
  gp << "set ylabel 'Posición Y'" << endl;
  gp << "plot for [i=0:n-1] 'PenduloXvsY_K="<<oss.str()<<".txt' every n::i using 1:2 with linespoints pt 7 ps 1 title sprintf('Péndulo %d', i+1)" << endl;
  gp << "set output" << endl;
  gp.close();
}



//-----------------Funcion numeral e -----------------------

void DatosTorqueRescalado(Cuerpo & pendulo, double t, double t0, double K, 
                          double a, double b, ofstream & archivo){
  
  // Calcula el tiempo y el torque reescalados
  double tiempo_reescalado = (t - t0) * pow(K, - b);
  double torque_reescalado = pendulo.GetTau() * pow(K, - a);
  
  // Write the rescaled values of time and torque to the file
  archivo << tiempo_reescalado << " " << torque_reescalado << endl;
 
}

void grafica_e(double K) {

    ostringstream ossK, ossa;
    ossK << scientific << setprecision(1) << static_cast<double>(K);

    string filename = "graph_e.gp";
    std::ofstream gp_file(filename.c_str());

    if (!gp_file.is_open()) {
        std::cerr << "No se pudo abrir el archivo para escribir el script de Gnuplot." << std::endl;
        return;
    }

    gp_file << "K=" << K << endl;
    gp_file << "set terminal png size 800,600" << std::endl;
    gp_file << "set output 'CunaNewtonTauRescalado.png'" << std::endl;
    gp_file << "set title 'Time vs Tau rescalado para K=" + ossK.str() + "'" << std::endl;
    gp_file << "set xlabel 'Time (s)'" << std::endl;
    gp_file << "set ylabel 'Tau (mx^2s^-2))'" << std::endl;
    gp_file << "set grid" << std::endl;
    gp_file << "set style data linespoints" << std::endl;

    gp_file << "plot ";
    gp_file << "'CunaNewtonTauRescalado_K=" + ossK.str() + ".txt' using 1:2 title sprintf('K = %2f e10, a = %.2f', K, a)";
  
    gp_file << std::endl;
    gp_file << "set output" << std::endl;

    gp_file.close();
    
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
  
  float escale = 1e9;
  float K = atof(argv[1])*escale;
  clog<<"\n\n\nK = "<<K<<"\n\n\n";

  double a = 0.400023 , b = -0.388982;

  Cuerpo pendulos[N];
  Colisionador Newton;
  
  //Archivos para guardar los datos
  ostringstream ossk, ossa, ossb;
  ossk << scientific << setprecision(1) << static_cast<double>(K);
  ossa << scientific << setprecision(1) << static_cast<double>(a);
  ossb << scientific << setprecision(1) << static_cast<double>(b);

  ofstream archivoXvsY; // Open the file in append mode
  ofstream archivoTauKRes;
  
  

    archivoXvsY.open("PenduloXvsY_K=" + ossk.str() + ".txt");
    archivoTauKRes.open("CunaNewtonTauRescalado_K=" + ossk.str()  + ".txt");
      
    //Variables para el algoritmo
    double m0=100, R0=1.5;
  double theta01=-M_PI/12;
  double theta0i=0;
    double T=2*M_PI*sqrt(L/g);
    double ajuste = 0.2;
    double T_ajuste = ajuste*T;
    double T_total = T + T_ajuste;
    // double t,tstart = 0.174,tmax=0.184,dt=1e-5;
    int i, Ncuadros=100;
    double t,tstart = 0.174,tdibujo,dt=1e-5,tmax=2*T_total ,tcuadro=tmax/Ncuadros; 

    double x0=0,y0=L,omega0=0;
    double dx = 2*R0;//Separación entre pendulos 
    
    //Variables para medir el rendimiento
    double mean_wtime, sigma_wtime;
    double mean_ctime, sigma_ctime;
    double wsum = 0, wsum2 = 0, csum = 0, csum2 = 0;
    double wtime, ctime;
    double reps = tmax/dt;
    auto start = std::chrono::system_clock::now(); // measures wall time
    std::clock_t c1 = std::clock();

    

  // InicieAnimacion(K);
  // graficar(K);
  
  // grafica_e(K);


  //---------------(Theta0,Omega0,m0,R0,L0,x00)
  pendulos[0].Inicie(x0,y0,theta01,omega0,m0,R0);
  
  //Inicializacion de los demas pendulos
  for(i=1;i<N;i++){
    x0=i*dx;
    pendulos[i].Inicie(x0,y0,theta0i,omega0,m0,R0);
  }
  

  
 
  for(t=0; t<tmax; t+=dt){

    // if(tdibujo>=tcuadro){
      
    //   InicieCuadro(t);
    //   for(i=0;i<N;i++) {
    //     pendulos[i].Dibujese();
    //     Datos(pendulos[i], K); 
    //     }
    //   TermineCuadro();
      
    //   tdibujo=0;
    // }

    if(t > tstart && t < tstart + 0.01) { 
          DatosTorqueRescalado(pendulos[1], t, tstart, K, a, b, archivoTauKRes);
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


    stats(wsum, wsum2, csum, csum2, wtime, ctime, reps, mean_wtime, sigma_wtime, mean_ctime, sigma_ctime, K);
    
    //Cierre de archivos
    archivoXvsY.close();
    archivoTauKRes.close();  
  
  return 0;
}
