#include <iostream>
#include <cmath>
#include "vector.h"
#include "Random64.h"
#include <algorithm>
#include <fstream>
#include <set>
using namespace std;

//---- declarar constantes ---
const double K=1.0e4;
const double Lx=60, Ly=60;
const int Nx=5, Ny=5, N=Nx*Ny;
const double m = 1.0, Ep = 1.0, r_0 = 10.0, R0=2.5;

const double epsilon=0.1786178958448091e00;
const double lambda=-0.2123418310626054e00;
const double chi=-0.6626458266981849e-1;
const double lambda2=(1.0-2.0*lambda)/2.0;
const double chiepsilon=1.0-2.0*(chi+epsilon);

//-------------------------------- Clases -----
class Cuerpo;
class Colisionador;
//-------------------------------- Clase cuerpo ---
class Cuerpo{
private:
  vector3D r,V,F; double m,R; double theta,omega,tau,I;
	int flag=0;
public:
  void Inicie(double x0,double y0,double Vx0,double Vy0,
	      double theta0,double omega0,double m0,double R0);
  void BorreFuerza(){F.load(0,0,0);};
  void AdicioneFuerza(vector3D F0){F+=F0;};
  void Mueva_r(double dt, double Coeficiente);
  void Mueva_V(double dt, double Coeficiente);
  void Dibujese(void);
	void Wall(int n){flag = n; }; 
  double Getx(void){return r.x();}; //inline
  double Gety(void){return r.y();}; //inline
	double GetVx(void){return V.x();}; //inline
	double GetVy(void){return V.y();}; //inline
	int GetFlag(void){return flag;}; //inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double x0,double y0,double Vx0,double Vy0,
		    double theta0,double omega0,double m0,double R0){
  r.load(x0,y0,0); V.load(Vx0,Vy0,0);  m=m0;  R=R0; flag = 0;
  theta=theta0; omega=omega0; I=2.0/5*m*R*R;
} 
void Cuerpo::Mueva_r(double dt, double Coeficiente){
  r+=V*(Coeficiente*dt); theta+=omega*(Coeficiente*dt);
}
void Cuerpo::Mueva_V(double dt, double Coeficiente){
  V+=F*(Coeficiente*dt/m); omega+=tau*(Coeficiente*dt/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)"; // ","
    //<<r.x()<<"+"<<R*cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*sin(theta)/7.0<<"*t";
}

//--- clase Colisionador ----
class Colisionador{
private:
public:
  void CalculeFuerzas(Cuerpo * Molecula);
  void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2);
};

void Colisionador::CalculeFuerzas(Cuerpo * Molecula){
   int i,j;
	vector3D Faux;
  //Borrar fuerzas
  for(i=0;i<N;i++){Molecula[i].BorreFuerza();}
  //Calcular las fuerzas entre todas las parejas de Moleculas
  for(i=0;i<N;i++){
		double h;
		// Fuerza entre paredes
		if (Molecula[i].Getx() > Lx -R0 ){ // Pared derecha
		 h = Molecula[i].Getx()-(Lx -R0);
		Faux.load((-1)*K*pow(abs(h),1.5),0,0);
		Molecula[i].AdicioneFuerza(Faux); 
		Molecula[i].Wall(1);
		}
		if (Molecula[i].Getx() < R0 ){ // Pared izquierda
		h = R0- Molecula[i].Getx() ;
		Faux.load(K*pow(abs(h),1.5),0,0);
		Molecula[i].AdicioneFuerza(Faux); 
		Molecula[i].Wall(1);
		}
		if (Molecula[i].Gety()  < R0 ){ //Pared de abajo
		h = R0 - Molecula[i].Gety();
		Faux.load(0,K*pow(abs(h),1.5),0);
		Molecula[i].AdicioneFuerza(Faux); 
		Molecula[i].Wall(2);
		}
		if (Molecula[i].Gety()  > (Ly - R0) ){ // Pared de arriba
		h = Molecula[i].Gety() -(Ly - R0);
		Faux.load(0,(-1)*K*pow(abs(h),1.5),0);
		Molecula[i].AdicioneFuerza(Faux); 
		Molecula[i].Wall(2);
		}
		
    for(j=i+1;j<N;j++){CalculeFuerzaEntre(Molecula[i],Molecula[j]);}
		}
	}
	
void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2){
   vector3D r21,n,F1; double d21,F;
  r21=Molecula2.r-Molecula1.r; d21=r21.norm(); n=r21/d21;
  F=(12*Ep/(d21*d21))*((pow((r_0/(d21)),12))-(pow(r_0/(d21),6)));
	F1 = F*n;
  Molecula2.AdicioneFuerza(F1);   Molecula1.AdicioneFuerza(F1*(-1));
  }   
//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
   cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'k_01.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-10:"<<Lx+10<<"]"<<endl;
  cout<<"set yrange[-10:"<<Ly+10<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
    cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
}
void TermineCuadro(void){
    cout<<endl;
}

void Grafica(void){
  cout<<"set term pdf"<<endl; 
  cout<<"set out 'I_PvsT.pdf'"<<endl;
	cout<<"set title 'Presión vs Temperatura'"<<endl;
  cout<<"set ylabel 'Presion'"<<endl;
	cout<<"set xlabel 'Temperatura'"<<endl;
	cout<<"set autoscale"<<endl;
	cout<<"set key"<<endl;
	cout<<"set font ',7'"<<endl;
	cout<<"plot 'I.txt' u 1:2 w l"<<endl;
}

//-----------  Programa Principal --------------  
int main(void){
	ofstream I_data;
  I_data.open ("I.txt");

	double KBT[10] = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
10.0, 15.0, 20.0};

	//for(int ii = 0; ii< 10; ii++){
  Cuerpo Molecula[N]; 
  Colisionador Hertz;
  Crandom ran64(1);
  double m0=1.0, R0=2.5,  kT=20.0, V0=sqrt(2*kT/m0);
  int i,ix,iy;
  double t,tdibujo,dt=1e-4,tmax=200.0,tcuadro=tmax/280; // 280 500
  double dx=Lx/(Nx+1), dy=Lx/(Ny+1);
  double Theta, OmegaMax=1.0;
	double Vel1x[N] = {0.0};
	double Vel1y[N] = {0.0};
	double intensidad = 0.0;
	set<double> svel;
	int ti = 0;
	double Vel = 0.0;

  //Inicializar las moléculas
  for(ix=0;ix<Nx;ix++){
    for(iy=0;iy<Ny;iy++){
      Theta=2*M_PI*ran64.r();
      //------(x0,y0,Vx0,Vy0, theta0,omega0  ,m0,R0)
      Molecula[Nx*iy+ix].Inicie((ix+1)*dx,(iy+1)*dy, V0,0,0,OmegaMax,m0,R0);//OJO
    }}

	////---------------///
  for(t=0,tdibujo=0 ; t<tmax ; t+=dt,tdibujo+=dt){

		for(int ii=0;ii<N;ii++){
		Vel1x[ii] = Vel1y[ii] =0.0; // Borrar velocidades
		Molecula[ii].Wall(0); // Reset flag
		Vel1x[ii] = Molecula[ii].GetVx();
		Vel1y[ii] = Molecula[ii].GetVy();
		}
    //Dibujar
    //--- Muevase por PEFRL ---
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chiepsilon);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,chi);
    Hertz.CalculeFuerzas(Molecula);
    for(i=0;i<N;i++)Molecula[i].Mueva_V(dt,lambda2);
    for(i=0;i<N;i++)Molecula[i].Mueva_r(dt,epsilon);  

		for(int ii=0;ii<N;ii++){
			if (Molecula[ii].GetFlag()==1){
				intensidad += abs( Molecula[ii].GetVx() - Vel1x[ii]);
				}
			else if (Molecula[ii].GetFlag()==2){
				intensidad += abs(Molecula[ii].GetVy() - Vel1y[ii]); 
				}
		}

		//Se Guardan las Velocidades
			if(tdibujo>tcuadro){ 
				if(ti >=80){
      		for(i=0;i<N;i++){ 
					svel.insert(Molecula[i].GetVx());
				}
			}
      tdibujo=0; ti+=1;
    }
  }    

	I_data  << intensidad/(200*6) << "\t"; 

	// Calcular la DESVIACION ESTANDAR
  double vprom,v2prom,sigma_v;
  //Calculo vprom, v2prom
	for(auto it=svel.begin(); it != svel.end(); it++){
		vprom+=*it;
		v2prom+=pow(*it,2);
  }
	vprom/=svel.size();
	v2prom/=svel.size(); 
  //Calculo sigma_v
  sigma_v=sqrt(v2prom-vprom*vprom);
  //Imprimir sigma_v
  I_data<<sigma_v*sigma_v<<endl;
		//}
	I_data.close();
	Grafica();
  return 0;
}
