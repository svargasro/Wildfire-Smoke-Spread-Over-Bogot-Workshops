#include <iostream>
#include <cmath>
#include <tuple>

const int Lx=256, Ly=256;
const double p0=0.25,p=0.25;
const double pr = 1-2*p-p0;
const int Q=4;
const int N_automatas=1;

class LatticeGas{
private:
    int V[Q]; //V[i] i = 0 (izquierda) i = 1 (arriba) i = 2 (derecha) i = 3 (abajo)
    double f[Lx][Ly][Q], f_new[Lx][Ly][Q]; //n[ix][i] 
public:
    LatticeGas(void);
    void Borrese(void);
    void Inicie(int N, double mu_x, double mu_y, double sigma_x, double sigma_y);
    double rho(int ix, int iy, bool UseNew);
    void Colisione(void);
    void Adveccione(void);
    void Muestrese(void);
};

LatticeGas::LatticeGas(void){
    /*Convencion de direcciones:
    D2Q5
        1   
        | 
    0 - . - 2
        | 
        3   
    */
  //Definir los vectores velocidad
  //Derecha e izquierda
  V[2]=1;  V[0]=-1;
  //Arriba y abajo
  V[1]=1;  V[3]=-1;
}

void LatticeGas::Borrese(void){
    for(int ix=0;ix<Lx;ix++){
      for(int iy=0;iy<Ly;iy++){
          for(int i=0;i<Q;i++){
            f[ix][iy][i]=0;
        }
      }    
    }
}

void LatticeGas::Inicie(int N, double mu_x, double mu_y, double sigma_x, double sigma_y){
    
    double gauss_x, gauss_y;
    for(int ix=0;ix<Lx;ix++){
      gauss_x=exp(-0.5*pow((ix-mu_x)/sigma_x,2))/(sigma_x*sqrt(2*M_PI));
      for(int iy=0;iy<Ly;iy++){
          gauss_y=exp(-0.5*pow((iy-mu_y)/sigma_y,2))/(sigma_y*sqrt(2*M_PI));
          for(int i=0;i<Q;i++){
            f[ix][iy][i]=f_new[ix][iy][i]=(N/Q)*gauss_x*gauss_y;// La mitad de las partículas se iran a la derecha y la otra mitad a la izquierda
        }
      }
    }
}

double LatticeGas::rho(int ix, int iy, bool UseNew){
double sum=0.0;

  if(UseNew) for(int i=0;i<Q;i++) sum += f_new[ix][iy][i]; //Revisar scope
  else  for(int i=0;i<Q;i++) sum += f[ix][iy][i];
return sum;
}


void LatticeGas::Colisione(void){//Quede aquí  ----------------------------------------------
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++){ //Recorriendo cada celda
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){ //Recorriendo cada dirección
        f_new[ix][iy][i]= p0*f[ix][iy][i] + p*f[ix][iy][(i+1)%Q] 
        + p*f[ix][iy][(i+2)%Q] + pr*f[ix][iy][(i+3)%Q];
      }
    }
  }
}

void LatticeGas::Adveccione(void){
  int ix,iy,i;
    for(ix=0;ix<Lx;ix++){//Para cada celda 
      for(iy=0;iy<Ly;iy++){
        for(i=0;i<Q;i++){ //Y en cada dirección 
        if (i == 2 || i == 0)//Derecha e izquierda
            f[(ix+V[i]+Lx)%Lx][iy][i]=f_new[ix][iy][i];   
        else //Arriba y abajo
            f[ix][(iy+V[i]+Lx)%Lx][i]=f_new[ix][iy][i];
        }
      } 
    }
}


//------------------- FUNCIONES GLOBALES -------
std::tuple<double, double, double> Sigma2(LatticeGas * Difusion, bool last_i){
  int ix,iy,iautomata;
  bool UseNew=false;
  // bool UseNew = true;
  //Calcular cuántas bolitas hay
  double N=0.0, xprom=0.0, yprom=0.0;
  
  for(iautomata=0;iautomata<N_automatas;iautomata++){
    for(ix=0;ix<Lx;ix++){
      for(iy=0;iy<Ly;iy++){
        N+=Difusion[iautomata].rho(ix,iy,UseNew);
        xprom+=ix*Difusion[iautomata].rho(ix,iy,UseNew);
        yprom+=iy*Difusion[iautomata].rho(ix,iy,UseNew);
       } } }
  //Calcular la posición promedio
  xprom/=N;  yprom/=N;

  //Calcular la varianza promedio
  double sigma2_x=0.0, sigma2_y=0.0;
  double sigma2;
  for(iautomata=0;iautomata<N_automatas;iautomata++)
    for(ix=0;ix<Lx;ix++)
      for(iy=0;iy<Ly;iy++){
        sigma2_x+=pow(ix-xprom,2.0)*Difusion[iautomata].rho(ix,iy,UseNew);
        sigma2_y+=pow(iy-yprom,2.0)*Difusion[iautomata].rho(ix,iy,UseNew);
      }
  //clog<<"Sigma2: "<<Sigma2<<endl;
  sigma2_x/=(N-1);  sigma2_y/=(N-1);
  sigma2 = sigma2_x + sigma2_y;
   //Exportar datos distribucion final
  if(last_i) std::clog<<xprom<<"\t"<<yprom<<"\t"<<sqrt(sigma2_x)<<"\t"<<sqrt(sigma2_y)<<std::endl;
  return std::make_tuple(sigma2_x, sigma2_y, sigma2);
}

int main(void){
  LatticeGas Difusion[N_automatas];
  int N_par=2400; //N particulas
  double mu0_xi=Lx/2, sigma0_xi=16; //Centro y desviación estándar iniciales
  double mu0_yi=Ly/2, sigma0_yi=16;
  int tmax=350;

  for (int i = 0; i < N_automatas; i++) Difusion[i].Borrese();
  for (int i = 0; i < N_automatas; i++) Difusion[i].Inicie(N_par, mu0_xi, mu0_yi, sigma0_xi, sigma0_yi);

  //Exportar datos distribucion incial
  std::clog<<mu0_xi<<"\t"<<mu0_yi<<"\t"<<sqrt(sigma0_xi)<<"\t"<<sqrt(sigma0_yi)<<std::endl;

  for(int t=0;t<tmax;t++){
    double sigma_x, sigma_y, sigma2;
    std::tie(sigma_x, sigma_y,sigma2) = Sigma2(Difusion,false);
    std::cout<<t<<"\t"<<sigma_x<<"\t"<<sigma_y<<"\t"<<sigma2<<std::endl;
    for (int i = 0; i < N_automatas; i++) Difusion[i].Colisione();
    for (int i = 0; i < N_automatas; i++) Difusion[i].Adveccione();
    if(t == tmax-1) Sigma2(Difusion,true);
  }

  return 0;
}