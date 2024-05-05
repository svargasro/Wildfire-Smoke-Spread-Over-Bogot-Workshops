#include <iostream>
#include <cmath>
#include <tuple>
#include "../../Random64.h"

const int Lx=256, Ly=256;
const double p0=0.25,p=0.25;
const int Q=4;
const int N_automatas=1;

class LatticeGas{
private:
    int V[Q]; //V[i] i = 0 (izquierda) i = 1 (arriba) i = 2 (derecha) i = 3 (abajo)
    int n[Lx][Ly][Q], n_new[Lx][Ly][Q]; //n[ix][i] 
public:
    LatticeGas(void);
    void Borrese(void);
    void Inicie(int N, double mu_x, double mu_y, double sigma_x, double sigma_y, Crandom&ran64);
    double rho(int ix, int iy, bool UseNew);
    void Colisione(Crandom& ran64);
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
            n[ix][iy][i]=0;
        }
      }    
    }
}

void LatticeGas::Inicie(int N, double mu_x, double mu_y, double sigma_x, double sigma_y,Crandom&ran64){
    int ix,iy,i;
    while(N>0){
        //Escojo un sitio al azar usando una distribución gaussiana
        ix=(int) ran64.gauss(mu_x,sigma_x);// Tomar la parte entera
        iy=(int) ran64.gauss(mu_y,sigma_y);// Tomar la parte entera
        if(ix<0) ix=0; 
        if(ix>=Lx) ix=Lx-1;

        if(iy<0) iy=0;
        if(iy>=Ly) iy=Ly-1;
        //Escojo al azar entre las dos direcciones 
        i =(int) Q*ran64.r();
        //Si ese sitio esta vacio, lo lleno y decremento n en 1
        if(n[ix][iy][i]==0){n[ix][iy][i]=1; N--;}
    }
}

double LatticeGas::rho(int ix, int iy, bool UseNew){
double sum=0.0;

  if(UseNew) for(int i=0;i<Q;i++) sum += n_new[ix][iy][i]; //Revisar scope
  else  for(int i=0;i<Q;i++) sum += n[ix][iy][i];
return sum;
}


void LatticeGas::Colisione(Crandom & ran64){//Quede aquí  ----------------------------------------------
  int ix,iy,i;
  double n_aleatorio;
  for(ix=0;ix<Lx;ix++){ //Recorriendo cada celda
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){       //Recorriendo cada dirección 
      n_aleatorio = ran64.r();
      //Se rota la celda 0 grados (quieta)
      if(n_aleatorio<p0)n_new[ix][iy][i]=n[ix][iy][i];
      //Se rota la celda 90 grados (derecha)
      else if(n_aleatorio<p0+p)n_new[ix][iy][i]=n[ix][iy][(i+3)%Q];
      //Se rota la celda 270 grados (izquierda)
      else if(n_aleatorio<p0+2*p)n_new[ix][iy][i]=n[ix][iy][(i+1)%Q];
      //Se rota la celda 180 grados (abajo)
      else n_new[ix][iy][i]=n[ix][iy][(i+2)%Q];
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
            n[(ix+V[i]+Lx)%Lx][iy][i]=n_new[ix][iy][i];   
        else //Arriba y abajo
            n[ix][(iy+V[i]+Lx)%Lx][i]=n_new[ix][iy][i];
        }
      } 
    }
}


//------------------- FUNCIONES GLOBALES -------
tuple<double, double, double> Sigma2(LatticeGas * Difusion){
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
  clog<<"N: "<<N<<endl;
  // N = 2400;
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
  
  return make_tuple(sigma2_x, sigma2_y, sigma2);
}

int main(void){
  LatticeGas Difusion[N_automatas];
  Crandom ran64(1);
  int N_par=2400; //N particulas
  double mu_xi=Lx/2, sigma_xi=16; //Centro y desviación estándar iniciales
  double mu_yi=Ly/2, sigma_yi=16;
  int tmax=350;

  for (int i = 0; i < N_automatas; i++) Difusion[i].Borrese();
  for (int i = 0; i < N_automatas; i++) Difusion[i].Inicie(N_par, mu_xi, mu_yi, sigma_xi, sigma_yi, ran64);

  for(int t=0;t<tmax;t++){
    double sigma_x, sigma_y, sigma2;
    std::tie(sigma_x, sigma_y,sigma2) = Sigma2(Difusion);
    cout<<t<<"\t"<<sigma_x<<"\t"<<sigma_y<<"\t"<<sigma2<<endl;
    // cout<<t<<"\t"<<sigma2<<endl;   
    for (int i = 0; i < N_automatas; i++) Difusion[i].Colisione(ran64);
    for (int i = 0; i < N_automatas; i++) Difusion[i].Adveccione();
  }
  return 0;
}