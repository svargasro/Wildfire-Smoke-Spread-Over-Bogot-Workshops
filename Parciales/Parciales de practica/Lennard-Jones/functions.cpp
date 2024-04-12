#include "functions.hpp"


//--------------------Funciones de la clase Cuerpo--------------------
void Cuerpo::Inicie(double x0, double y0, double Vx0, double Vy0,
                    double m0, double R0){
                    r.load(x0,y0,0); v.load(Vx0,Vy0,0); m = m0; R = R0;
                    }
void Cuerpo::Mueva_r(double dt, double coeficiente){
    r += v*(coeficiente*dt);
}
void Cuerpo::Mueva_v(double dt, double coeficiente){
    v += F*(coeficiente*dt/m);
}

void Cuerpo::Dibujese(void){
    std::cout<<" , "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//--------------------Funciones de la clase Colisionador--------------------

void Colisionador::CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2, double epsilon, double r0){
    r21 = Molecula2.r - Molecula1.r;
    double r = r21.norm();
    double F = 12*epsilon*(std::pow(r0/r,12) - std::pow(r0/r,6))/r;
    F1 = r21*F/r;
    F2 = -1*r21*F/r;

    Molecula1.F += F1; Molecula2.F += F2;
}

void Colisionador::CalculeTodasLasFuerzas(Cuerpo * moleculas, double epsilon, double r0, int N){
    for(int i = 0; i < N; i++){
        moleculas[i].F.load(0,0,0);
    }

    // for(int i = 0; i < N; i++){
    //     for(int j = i + 1; j < N; j++){
    //         CalculeFuerzaEntre(moleculas[i], moleculas[j], epsilon, r0);
    //     }       

    for (int i = 0; i< N; i++){
        for( int j = 0; j < i; j++){
            CalculeFuerzaEntre(moleculas[i], moleculas[j], epsilon, r0);
        }
    }
}



void InicieAnimacion(void){
    std::cout<<"set terminal gif animate\n";
    std::cout<<"set output 'animacion.gif'\n";
    std::cout<<"unset key\n";
    std::cout<<"set xrange [-10:10]\n";
    std::cout<<"set yrange [-10:10]\n";
    std::cout<<"set size square\n";
    std::cout<<"set parametric\n";
    std::cout<<"set trange [0:7]\n";
}
void InicieCuadro(void){
    std::cout<<"plot 0,0";
}
void TermineCuadro(void){
    std::cout<<std::endl;
}