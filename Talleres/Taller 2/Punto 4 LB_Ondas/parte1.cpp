#include <iostream>
#include <fstream>
#include <cmath>

//---------------Constantes globales-------------------
const int Lx=400;
const int Ly=200;
const int Q=5; // Numero de direcciones 
const double W0=1.0/3; // Peso de la direccion cero

const double tau=0.5; //Tiempo de relajacion
const double Utau=1.0/tau;
const double UmUtau=1-Utau;
//---------------Clase Lattice Boltzmann-------------------

class LatticeBoltzmann{
  private:
    double w[Q];      //Peso de cada direccion
    int Vx[Q],Vy[Q];  //Componentes de velocidad
    double *f, *fnew; //Funciones de distribucion
  public:
    LatticeBoltzmann(void); //Constructor
    ~LatticeBoltzmann(void); //Destructor
    int n(int ix,int iy,int i){return (ix*Ly+iy)*Q+i;};
    //Campos macroscopicos
    double rho(int ix,int iy,bool UseNew); //Calcula la densidad en la celda
    double Jx(int ix,int iy,bool UseNew); //Calcula la corriente en x
    double Jy(int ix,int iy,bool UseNew); //Calcula la corriente en y
    double Ccell(int ix, int iy); //Calcula la velocidad del sonido en la celda
    //Funciones de equilibrio
    double feq(double rho0,double Jx0,double Jy0,int i,int ix, int iy);
    //Evolucion temporal
    void Start(double rho0,double Jx0,double Jy0);
    //Reglas de evolución
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    //Impresión de resultados
    void Print(const char * NameFile);  
};
LatticeBoltzmann::LatticeBoltzmann(void){
  //Inicializando los pesos
  w[0]=W0; w[1]=w[2]=w[3]=w[4]=(1.0-W0)/4;
  //Inicializando los vectores de velocidad
  Vx[0]=0;  Vx[1]=1;  Vx[2]=0;  Vx[3]=-1; Vx[4]=0;
  Vy[0]=0;  Vy[1]=0;  Vy[2]=1;  Vy[3]=0;  Vy[4]=-1;
  //Creando las matrices f y fnew, se guardan en forma unidimensional en una lista
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
double  LatticeBoltzmann::feq(double rho0,double Jx0,double Jy0,int i,int ix, int iy){
  double C2 = pow(Ccell(ix,iy),2);
  if(i>0)
    return 3*w[i]*(C2*rho0+Vx[i]*Jx0+Vy[i]*Jy0);
  else
    return rho0*(1-3*C2*(1-W0));
}

void LatticeBoltzmann::Start(double rho0,double Jx0,double Jy0){
  int ix,iy,i,n0;
  for(ix=0;ix<Lx;ix++) //En cada celda
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){ //En cada direccion
	      n0=n(ix,iy,i);
	      f[n0]=feq(rho0,Jx0,Jy0,i,ix,iy);
      }
}

void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0;
  double rho0,Jx0,Jy0;
  for(ix=0;ix<Lx;ix++) //En cada celda
    for(iy=0;iy<Ly;iy++){
      //Calculo de los campos macroscopicos
      rho0=rho(ix,iy,false); Jx0=Jx(ix,iy,false); Jy0=Jy(ix,iy,false);
      for(i=0;i<Q;i++){ //En cada direccion
        n0=n(ix,iy,i);
        fnew[n0]=UmUtau*f[n0]+Utau*feq(rho0,Jx0,Jy0,i,ix,iy);
      }
    }
 }
void LatticeBoltzmann::ImposeFields(int t){
  //Ponemos una fuente plana en ix=0
  int ix=0,i,iy,n0;
  double lambda=10.0,omega, rho0;
  for(iy=0;iy<Ly;iy++){
    omega=(2*M_PI/lambda)*Ccell(ix,iy);
    rho0=10*sin(omega*t);
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(rho0,Jx(ix,iy,false),Jy(ix,iy,false),i,ix,iy);
      }
    }
  }
  
void LatticeBoltzmann::Advection(void){
  int ix,iy,ixnext, iynet,n0,n0next;
  for(ix=0;ix<Lx;ix++) //En cada celda
    for(iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){ //En cada direccion
        n0=n(ix,iy,i);
        ixnext=(ix+Vx[i]+Lx)%Lx;
        iynet=(iy+Vy[i]+Ly)%Ly;
        // ixnext=(ix+Vx[i]);
        // iynet=(iy+Vy[i]);
        // if (ixnext==Lx || ixnext<0) ixnext=ix;       
        // if (iynet==Ly || iynet<0) iynet=iy;
        n0next=n(ixnext,iynet,i);
        f[n0next]=fnew[n0];
      }
}
double LatticeBoltzmann::Ccell(int ix, int iy){
  double C = 0.5; // Velocidad del sonido [cells/click]
  //double ix0 = 100; // Posicion de la interfaz punto 4
  double theta = 20; //Angulo de incidencia en grados
  double ix0 = 100 + std::tan(theta*M_PI/180.0)*(iy-Ly/2);
  double n_i=0.5*std::tanh(ix-ix0)+1.5; //tanh para suavizar el paso entre los medios 
  double C_i = C/n_i; //Velocidad del sonido en el medio i
  return C_i; //retornar C para el punto 2
}


void LatticeBoltzmann::Print(const char * NameFile){
  std::ofstream MyFile(NameFile); double rho0; int ix,iy;
  for(ix=0;ix<Lx/2;ix++){ //Para esta configuración vamos a graficar hasta 200
    for(iy=0;iy<Ly;iy++){
      rho0=rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<std::endl;
    }
    MyFile<<std::endl;
  }
  MyFile.close();
}

int main(void){
  LatticeBoltzmann Waves;
  double t,tmax=400;
  double rho0=0,Jx0=0,Jy0=0;
  Waves.Start(rho0,Jx0,Jy0);

  for(t=0;t<tmax;t++){
    Waves.Collision();
    Waves.ImposeFields(t);
    Waves.Advection();
    if(int(t)%int(tmax/10)==0)std::clog<<"\nPorcentaje de avance: "<<(t/tmax)*100<<"%"<<std::endl;
  }
  Waves.Print("parte1.dat");
  //gnuplot 
  //splot "parte1.dat" w l
  return 0;
}