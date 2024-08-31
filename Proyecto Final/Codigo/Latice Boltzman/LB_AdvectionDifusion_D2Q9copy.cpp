#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>  // iomap sirve para setw y setfill
#include <sstream>  // Para usar stringstream
using namespace std;

//------------------- CONSTANTES GLOBALES -------------
const int Lx=100;
const int Ly=100;

const int Q=9;

const double C=0.5;
const double C2 = C*C;
const double Cs2=0.57;

const double tau=0.55;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//------------------- FUNCIONES GLOBALES -------------
class LatticeBoltzmann{
private:
  double w[Q];      //Weights
  int Vx[Q],Vy[Q];  //Velocity vectors
  double *f, *fnew; //Distribution Functions
public:
  //----Constructor y destructor----
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
  //----Campos macroscopicos----
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  //----Funciones de equilibrio----
  double feq(double rho0,double Ux0,double Uy0,int i);
  //----Evolucion temporal----
  void Start(double A, double sigma0, double D, double ux, double uy, double x0, double y0);
  void Collision(void);
  void ImposeFields(void);
  void Advection(void);
  void Print(string NameFile, double t);
  void Printframe(double t);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Set the weights
  w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
  //Set the velocity vectors
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;

              Vx[5]=1;  Vx[6]=-1; Vx[7]=-1; Vx[8]=1;
              Vy[5]=1;  Vy[6]=1;  Vy[7]=-1; Vy[8]=-1;
  //Create the dynamic arrays
  int ArraySize=Lx*Ly*Q;
  f=new double [ArraySize];  fnew=new double [ArraySize];
}
LatticeBoltzmann::~LatticeBoltzmann(void){
    delete[] f;  delete[] fnew;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=fnew[n0]; else sum+=f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vx[i]*fnew[n0]; else sum+=Vx[i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  double sum; int i,n0;
  for(sum=0,i=0;i<Q;i++){
    n0=n(ix,iy,i);
    if(UseNew) sum+=Vy[i]*fnew[n0]; else sum+=Vy[i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*Vx[i]+Uy0*Vy[i], U2=Ux0*Ux0+Uy0*Uy0;
  return rho0 * w[i] * (1 + UdotVi / Cs2 + (UdotVi * UdotVi) / (2 * Cs2 * Cs2) - U2 / (2 * Cs2));
}
void LatticeBoltzmann::Start(double A, double sigma0, double D, double ux, double uy, double x0, double y0) {
    int ix, iy, i, n0;

    for (ix = 0; ix < Lx; ix++) {
        for (iy = 0; iy < Ly; iy++) {
            double rho = A * exp(-(pow(ix - x0, 2) + pow(iy - y0, 2)) / (2 * sigma0 * sigma0));
            for (i = 0; i < Q; i++) {
                n0 = n(ix, iy, i);
                f[n0] = feq(rho, ux, uy, i);
            }
        }
    }
}
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++){
      //compute the macroscopic fields on the cell
      rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){ //for each velocity vector
	n0=n(ix,iy,i);
	fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }
}


void LatticeBoltzmann::ImposeFields(void) {
    // Implementación para imponer un campo de velocidad constante
    double rho0, Ux0, Uy0;
    for (int ix = 0; ix < Lx; ix++) {
        for (int iy = 0; iy < Ly; iy++) {
            rho0 = rho(ix, iy, true); // Usar fnew para obtener la densidad
            Ux0 = 0.03; // Velocidad en x
            Uy0 = 0.03; // Velocidad en y
            for (int i = 0; i < Q; i++) {
                int n0 = n(ix, iy, i);
                fnew[n0] = feq(rho0, Ux0, Uy0, i);
            }
        }
    }
}


void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++) //for each cell
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //on each direction
	ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
	n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
	f[n0next]=fnew[n0]; //periodic boundaries
      }
}
void LatticeBoltzmann::Print(string NameFile, double t) {
    ofstream MyFile(NameFile); 
    double rho0, Ux0, Uy0;
    for (int ix = 0; ix < Lx; ix += 4) {
        for (int iy = 0; iy < Ly; iy += 4) {
            rho0 = rho(ix, iy, true);
            // Ux0 = Jx(ix, iy, true) / rho0; 
            // Uy0 = Jy(ix, iy, true) / rho0;
            // Imprimimos cada 4 celdas la densidad
            MyFile << ix << " " << iy << " " << rho0 << endl;
        }
        MyFile << endl;
    }
    MyFile.close();
}



void LatticeBoltzmann::Printframe(double t) {
    ofstream GnuplotScript("frame_script.gp");
    GnuplotScript << "set terminal pngcairo size 800,800 enhanced font 'Verdana,10'" << endl;
    GnuplotScript << "set output 'frames/density_" << setw(3) << setfill('0') << t << ".png'" << endl;
    GnuplotScript << "set pm3d map" << endl;
    GnuplotScript << "set size ratio -1" << endl;
    GnuplotScript << "set xrange [0:" << Lx << "]" << endl;
    GnuplotScript << "set yrange [0:" << Ly << "]" << endl;
    GnuplotScript << "set cbrange [0:*]" << endl;
    GnuplotScript << "set palette defined (0 'black', 1 'red', 2 'yellow', 3 'white', 4 'red')" << endl;
    GnuplotScript << "set title 'Densidad en t = " << t << "'" << endl; 
    GnuplotScript << "plot 'data/density_" << setw(3) << setfill('0') << t << ".dat' u 1:2:3 w image" << endl;
    GnuplotScript.close();

    // Ejecutar el script de Gnuplot
    system("gnuplot frame_script.gp");
}

//--------------- FUNCIONES GLOBALES ------------
int main(void){
  LatticeBoltzmann Aire;
  int t,tframe = 10 ,tmax=1000;  // Duración de la simulación
  double A = 1000.0;  // Máxima densidad
  double sigma0 = 6.4;  // Desviación estándar inicial
  double D = 0.05;  // Coeficiente de difusión
  double ux = 0.03, uy = 0.03;  // Velocidades
  double x0 = 33.0, y0 = 33.0;  // Posición inicial

  // Inicialización
  Aire.Start(A, sigma0, D, ux, uy, x0, y0);
  // Ejecutar
  for(t=0;t<tmax;t++){
    Aire.Collision();
    Aire.ImposeFields();
    Aire.Advection();
        if (t % tframe == 0) {
            // Imprimir los datos en archivo
            stringstream ss;
            ss << "data/density_" << setw(3) << setfill('0') << t << ".dat";
            Aire.Print(ss.str(), t);
            // Generar y guardar el frame con Gnuplot
            Aire.Printframe(t);
        }
    }

  return 0;
}