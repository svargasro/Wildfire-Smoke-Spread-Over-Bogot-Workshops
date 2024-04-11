#pragma once // Directiva del preprocesador para evitar la doble inclusión de archivos de cabecera

#include <iostream>
#include "../../vector.h"
#include "../../random64.h"

//--------------------Declaracion de la clase Cuerpo-------------------
class Cuerpo{
private:
    vector3D r,v,F;
    double m,R;
public:
    // vector3D r,v,F;
    // double m,R;
    void Inicie(double x0,double y0, double Vx0, double Vy0
                ,double m0,double R0);
    void Mueva_r(double dt, double coeficiente);
    void Mueva_v(double dt, double coeficiente);
    void Dibujese(void);
    friend class Colisionador; //Para que Colisionador pueda acceder a r, V y m
};

//--------------------Declaracion de la clase Colisionador-------------------
class Colisionador{
    private:
    vector3D r21,F1,F2;
    public:
    void CalculeFuerzaEntre(Cuerpo & Molecula1, Cuerpo & Molecula2, double epsilon, double r0);
    void CalculeTodasLasFuerzas(Cuerpo * moleculas, double epsilon, double r0, int N);
};


//--------------------Funciones para la animacion -------------------

void InicieAnimacion(void);
void InicieCuadro(void);
void TermineCuadro(void);