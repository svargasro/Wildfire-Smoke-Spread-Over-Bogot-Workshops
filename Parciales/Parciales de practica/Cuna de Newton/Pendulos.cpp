#include <iostream>//Libreria para imprimir en la consola
#include <cmath>//Libreria para funciones matemáticas
#include "../vector.h"
#include "../Random64.h"
using namespace std;

//Constantes del problema físico
const double g=9.8;

//Número de moleculass
const int Nx=3; 
const int Ny=1; 
const double Lx=60, Ly=60; 
const int N = Nx*Ny;

const double K=1.0e4; //Elasticidad de la colisión

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
  //Atribtuos de la clase Cuerpo
  //vector3D r, rp, V, F; double theta, m, R, L0;
  vector3D r, rp, F; 
  double m, R, L0;
  double theta, omega, alpha;
  vector3D Vaux, Fg;
public:
  vector3D V;
  void Inicie(double x0,double y0,double xp,double yp,
	      double Vx0,double Vy0,double m0,double R0);
  void BorreFuerza(void){F.load(0,0,0);};// Inline
  void SumeFuerza(vector3D dF){F+=dF;};// Inline
  void Mueva_r(double dt,double coeficiente);
  void Mueva_V(double dt,double coeficiente);
  void SetAlpha(double alpha0){alpha = alpha0;};
  void Dibujese(void);
  friend class Colisionador;
};
class Colisionador{
private:
public:
  void CalculeTodasLasFuerzas(Cuerpo * moleculas);
  void CalculeFuerzaEntre(Cuerpo & moleculas1,Cuerpo & moleculas2);
};

//-------Implementar las funciones de las clases------
//------- Funciones de la clase cuerpo --------
void Cuerpo::Inicie(double x0,double y0,double xp,double yp,
	      double Vx0,double Vy0,double m0,double R0){
  r.load(x0,y0,0);  V.load(Vx0,Vy0,0); m=m0; R=R0;
  rp.load(xp,yp,0); 
  L0 = (r - rp).norm();
}
void Cuerpo::Mueva_r(double dt,double coeficiente){
  // double L = (r - rp).norm();
  // double pre = 10;

  // // theta = atan2((r-rp).y(),(r-rp).x());

  // // clog<<"theta"<<theta<<endl;
  // clog<<"Vx="<<V.x()<<", Vy="<<V.y()<<endl;
  // // Vaux.load(V.x()*cos(theta),V.y()*sin(theta),0);
  // // clog<<"Vaux="<<Vaux<<endl;
  // if ((L >= L0 - pre) && (L <= L0 + pre) ) {
  //   r+=V*(coeficiente*dt);
  //   // r+=Vaux*(coeficiente*dt);
  // }
  
  theta += omega * dt * coeficiente;

  r.load(L0*cos(theta),L0*sin(theta),0);
 
  
}
void Cuerpo::Mueva_V(double dt,double coeficiente){

  omega += alpha * dt * coeficiente;
  V.load(-1*L0*omega*sin(theta),L0*omega*cos(theta),0);

  // V+=F*(coeficiente*dt/m);
  // theta = atan2((r-rp).y(),(r-rp).x());

  // vector3D Faux;
  // Faux.load(V.x()*cos(theta),V.y()*sin(theta),0);
 
  // V+=Faux*(coeficiente*dt);
}
// Dibuja un circulo de radio R en la posición del cuerpo en base al parametro t
// t va de 0 a 7 aproximadamente 2pi
void Cuerpo::Dibujese(void){
  theta = atan2((r-rp).y(),(r-rp).x());
  cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
      <<rp.x()<<"+"<<L0*cos(theta)/7.0<<"*t,"<<rp.y()<<"+"<<L0*sin(theta)/7.0<<"*t";
}
//------- Funciones de la clase Colisionador --------
void Colisionador::CalculeTodasLasFuerzas(Cuerpo * moleculas){
  
  double theta0,alpha0,L0;
  int i,j;
  //Borro las fuerzas de todos los moleculass
  for(i=0;i<N;i++)
    moleculas[i].BorreFuerza();

  //Recorro por parejas, calculo la fuerza de cada pareja y se la sumo a los dos
  for(i=0;i<N;i++){

    theta0 = moleculas[i].theta;
    L0 = moleculas[i].L0;
    moleculas[i].alpha = -g / L0 * sin(theta0);

    for(j=0;j<i;j++)
      CalculeFuerzaEntre(moleculas[i],moleculas[j]);

  }
    
}

//Calcula la fuerza entre dos moleculas
void Colisionador::CalculeFuerzaEntre(Cuerpo & moleculas1, Cuerpo & moleculas2){
  
  //Calcular el vector normal
  vector3D r21=moleculas2.r-moleculas1.r; 
  double d=r21.norm();


  //Determinar si hay colisión
  double s = moleculas1.R + moleculas2.R - d;
  
  if(s>0){
    double aux= K*pow(s,1.5);

    //Calcular vector normal
    vector3D n = r21*(1.0/d);

    //Calcular la fuerza
    vector3D F=n*aux;

    //Sumar las fuerzas
    
    moleculas2.SumeFuerza(F);  moleculas1.SumeFuerza(F*(-1));

  }

 
}


//----------- Funciones Globales -----------
//---Funciones de Animacion---
void InicieAnimacion(void){
  cout<<"set terminal gif animate delay 50"<<endl; 
  cout<<"set output 'Pendulos.gif'"<<endl;
  cout<<"unset key"<<endl;//No dibuja la leyenda de los datos
  cout<<"set xrange[-20:"<<Lx + 20<<"]"<<endl;// Rango en x
  cout<<"set yrange[-20:"<<Ly + 20<<"]"<<endl;// Rango en y
  cout<<"set size ratio -1"<<endl; //Para que la escala en x y y sea la misma
  cout<<"set parametric"<<endl;// Coordenadas parametricas para x y y
  cout<<"set trange [0:7]"<<endl;//Rango del parametro t de 0 a 7 aprox 2pi
  cout<<"set isosamples 12"<<endl;//Establece el número de muestras de t 
  //tanto en X como en Y a 12. Esto determina cuántos puntos se calculan 
  //y se muestran en las gráficas paramétricas. 
}

//---Funciones de cuadros---
//Inicia el cuadro de la animación
void InicieCuadro(double tiempo, double velocidad){
    cout << "set label \'t=" << tiempo << "s    v_{pc}="<< velocidad <<"(m/s) \' at graph 0.5,0.9 center" << endl;
    cout<<"plot 0,0 ";
    // cout<<" , "<<Lx/7<<"*t,0";        //pared de abajo
    // cout<<" , "<<Lx/7<<"*t,"<<Ly;     //pared de arriba
    // cout<<" , 0,"<<Ly/7<<"*t";        //pared de la izquierda
    // cout<<" , "<<Lx<<","<<Ly/7<<"*t"; //pared de la derecha
    
}

//Termina el cuadro
void TermineCuadro(void){
    cout<<endl;
    cout<<"unset label" << endl;
}

int main(){
  
  Cuerpo moleculas[N+4];
  Colisionador Newton;
  Crandom ran64(1);//Para generar números aleatorios
  int i,ix,iy;//Variables auxiliares para los ciclos de inicializacion
  //Parametros de la simulación
  double m0=1; double R0=4;
  double kT=100; 
  //Variables auxiliares para la condición inicial
  double dx=Lx/(Nx+1),dy=Ly/(Ny+1);//Separación entre moleculas en x y en y
  double theta; double V0=10;//sqrt(kT/m0);//Velocidad inicial por el teorema de equipartición 
  //de la energía
  double x0,y0,xp,yp,Vx0,Vy0;
  //Variables auxiliares para correr la simulacion
  int Ncuadros=10; double t,tdibujo,dt=1e0,tmax=5,tcuadro=tmax/(Ncuadros+1); 
  // cout<<"tmax="<<tmax<<", tcuadro="<<tcuadro<<", dt="<<dt<<endl;
  // clog<<"Ncuadros="<<Ncuadros<<", tmax="<<tmax<<", tcuadro="<<tcuadro<<", dt="<<dt<<endl;

  //Variables auxiliares para las paredes
  double Rpared=10*Lx, Mpared=100*m0;

  InicieAnimacion();
  
  //INICIO
  //Inicializar las paredes
  //----------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  // moleculas[N].Inicie(Lx/2,Ly+Rpared, 0,  0, 0,  0,Mpared,Rpared,0); //Pared arriba
  // moleculas[N + 1].Inicie(Lx/2,-Rpared, 0,  0, 0,  0,Mpared,Rpared,0); //Pared abajo
  // moleculas[N + 2].Inicie(Lx+Rpared,Ly/2, 0,  0, 0,  0,Mpared,Rpared,0); //Pared derecha
  // moleculas[N + 3].Inicie(-Rpared,Ly/2, 0,  0, 0,  0,Mpared,Rpared,0); //Pared izquierda
  

  //Inicializar las moleculas
  //Solo inicializar 
  // for(iy=0;iy<Ny;iy++)
  //   for(ix=0;ix<Nx;ix++){
  //     theta=2*M_PI*ran64.r();
  //     x0=(ix+1)*dx; y0=(iy+1)*dy; //Vx0=V0*cos(theta); Vy0=V0*sin(theta);
  //     //----------------(x0,y0,z0,Vx0,Vy0,Vz0,m0,R0)
  //     moleculas[iy*Nx+ix].Inicie(x0,y0, 0,Vx0,Vy0,  0,m0,R0);	
  //   }
  // Vx0=V0; Vy0=0;
  // theta=M_PI/12;
  // x0=L0*cos(theta); y0=L0*sin(theta);

  for(ix=0;ix<Nx;ix++){
      // theta=2*M_PI*ran64.r();
      xp=(ix+1)*dx; yp=Ly;
      x0=(ix+1)*dx; y0=0;
      Vx0=V0*cos(theta); Vy0=V0*sin(theta);
      //----------------(x0,y0,xp,yp,Vx0,Vy0,m0,R0)
      moleculas[iy*Nx+ix].Inicie(x0,y0,xp,yp,Vx0,Vy0,m0,R0);
  } 

  //---------------(x0,y0,z0,Vx0,   Vy0,Vz0,m0,R0)
  // moleculas[0].Inicie(x0, 0, 0,  0, 0.5*V0,  0,m0,1.0);
  // moleculas[1].Inicie(x1, 0, 0,  0, 0.5*V1,  0,m1,0.5);
  //CORRO
  double v = V0;
  for(t=tdibujo=0;t<=tmax;t+=dt,tdibujo+=dt){

    if(tdibujo>=tcuadro){
      
      InicieCuadro(t,v);
      for(i=0;i<N;i++) {
        moleculas[i].Dibujese();
        if(i == N/2) v = moleculas[i].V.norm();
        }
      TermineCuadro();
      
      tdibujo=0;
    }
    // cout<<moleculas[1].Getx()<<" "<<moleculas[1].Gety()<<endl;
    // clog<<"\nTiempo: "<<t<<"\n"<<endl;

    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,xi);    
    Newton.CalculeTodasLasFuerzas(moleculas); 
    for(i=0;i<N;i++) moleculas[i].Mueva_V(dt,Um2lambdau2);

    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(moleculas); 
    for(i=0;i<N;i++) moleculas[i].Mueva_V(dt,lambda);
    
    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,Um2chiplusxi);
    Newton.CalculeTodasLasFuerzas(moleculas); 
    for(i=0;i<N;i++)moleculas[i].Mueva_V(dt,lambda);

    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,chi);
    Newton.CalculeTodasLasFuerzas(moleculas); 
    for(i=0;i<N;i++)moleculas[i].Mueva_V(dt,Um2lambdau2);

    for(i=0;i<N;i++) moleculas[i].Mueva_r(dt,xi);
    
  }
  return 0;

}
