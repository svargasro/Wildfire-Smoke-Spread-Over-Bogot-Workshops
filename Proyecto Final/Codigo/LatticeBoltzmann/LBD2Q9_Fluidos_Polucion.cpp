#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//------------------- CONSTANTES GLOBALES -------------
const int Lx=256;
const int Ly=256;

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
  void Start(double rho0,double Ux0,double Uy0,double mu_x,double mu_y, double sigma_x, double sigma_y);
  void Collision(void);
  void ImposeFields(int t);
  void Advection(void);
  void Print(const char * NameFile);
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
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0,double mu_x, double mu_y, double sigma_x, double sigma_y) {
    int ix, iy, i, n0;


    for (ix = 0; ix < Lx; ix++) {
        for (iy = 0; iy < Ly; iy++) {
            double gauss_x = exp(-0.5 * pow((ix - mu_x) / sigma_x, 2)) / (sigma_x * sqrt(2 * M_PI));
            double gauss_y = exp(-0.5 * pow((iy - mu_y) / sigma_y, 2)) / (sigma_y * sqrt(2 * M_PI));
            double rho = rho0 * gauss_x * gauss_y;
            for (i = 0; i < Q; i++) {
                n0 = n(ix, iy, i);
                f[n0] = feq(rho, Ux0, Uy0, i);
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

void LatticeBoltzmann::ImposeFields(int t){
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
void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0,Ux0,Uy0; int ix,iy;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true); Ux0=Jx(ix,iy,true)/rho0; Uy0=Jy(ix,iy,true)/rho0;
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}
//--------------- FUNCIONES GLOBALES ------------


int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=35;
  double rho0=1000.0, Ux0=0.1,Uy0=0.1;
  double mu_x = Lx / 2, mu_y = Ly / 2, sigma_x = Lx / 32, sigma_y = Ly / 32; // Parámetros de la distribución gaussiana

  //Start
  Aire.Start(rho0, Ux0, Uy0, mu_x, mu_y, sigma_x, sigma_y);
  //Run
  for(t=0;t<tmax;t++){
    Aire.Collision();
    Aire.ImposeFields(t);
    Aire.Advection();
  }
  //Show
  Aire.Print("intentando.dat");

return 0;
}
