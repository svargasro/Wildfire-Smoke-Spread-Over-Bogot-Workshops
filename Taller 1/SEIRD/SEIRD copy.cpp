#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

//Constantes 
double n=9000000.0;// Poblacion total
double b=0.166;// tasa de contagios 
double rb=0.5; //asintomaticos son menos infecciosos que la constante B
double e = 1.0/5.0;// tasa de exposición
double pa = 0.8;// fracción de las personas que van a ser asintomáticos
double m = 1./21.0;//tasa de recuperados 
double pd =0.03;//fracción de casos asintomáticos que mueren


//Funciones del Modelo Seird

double fSt( double t, double St, double Ia, double Is, double Et, double Rt, double Dt){
	
	return (double) -1*(St/n)*(b*rb*Ia+Is);
	
}


double fEt( double t, double St, double Ia, double Is, double Et, double Rt, double Dt){
	
	return (double) (St/n)*(b*rb*Ia+b*Is)-e*(pa*Et+(1-pa)*Et);
	
}


double fIa( double t, double St, double Ia, double Is, double Et, double Rt, double Dt){
	
	return (double) e*pa*Et-m*Ia;
	
}


double fIs( double t, double St, double Ia, double Is, double Et, double Rt, double Dt){
	
	return (double) e*(1-pa)*Et-m*((1-pd)*Is-pd*Is);
	
}


double fRt( double t, double St, double Ia, double Is, double Et, double Rt, double Dt){
	
	return (double) m*(1-pd)*Is + m*Ia;
	
}

double fDt( double t, double St, double Ia, double Is, double Et, double Rt, double Dt){
	
	return (double) m*pd*Is;
	
}

void RK4S(double ci[6], double a, double b , double h)  {
	
	int nf=6,n=0,d=0 ;
	double (*(fns[6]))( double t, double St, double Ia, double Is, double Et, double Rt, double Dt) = {fSt,fIa,fIs,fEt,fRt,fDt};
   
    
	n= (int)(b-a)/h;

	double t[n];//Arrglo la variabiale t
	double x[nf][n];//Arreglo de variables
	double k[nf][4]; //Arreglo con las funciones auxiliares
	
	
    for(int i = 0; i<nf; i++){ //Incialización de variables
		
		x[i][0] = ci[i];
		
	}
   
    cout<<" -----------------------------------------------------------------------------------------------\n"; 
    printf("| %3s |  %10s  |  %10s  |  %10s  |  %10s  |  %10s  |  %10s  |\n","t","S","Ia","Is","E","R","D");
    cout<<" -----------------------------------------------------------------------------------------------\n"; 
    
    for(int i = 0; i<n ; i++) {
    	
    	t[i]=a+i*h;
    
	    for(int j = 0; j<nf ; j++) {//Ciclo de funciones auxiliares 1
	    	
			k[j][0] = h*fns[j](t[i], x[0][i], x[1][i], x[2][i], x[3][i], x[4][i], x[5][i]);
			
		}
		
	
		for(int j = 0; j<nf ; j++) {//Ciclo de funciones auxiliares 2
	    	
			k[j][1] = h*fns[j](t[i] + h/2, x[0][i] + k[0][0]/2, x[1][i] + k[1][0]/2, x[2][i] + k[2][0]/2,
                                x[3][i] + k[3][0]/2 , x[4][i] + k[4][0]/2 , x[5][i] + k[5][0]/2);
		
		}
		
		for(int j = 0; j<nf ; j++) {//Ciclo de funciones auxiliares 3
	    	
			k[j][2] = h*fns[j](t[i] + h/2, x[0][i] + k[0][1]/2, x[1][i] + k[1][1]/2, x[2][i] + k[2][1]/2,
                                x[3][i] + k[3][1]/2, x[4][i] + k[4][1]/2, x[5][i] + k[5][1]/2);
	    		
		}
		
		for(int j = 0; j<nf ; j++) {//Ciclo de funciones auxiliares 4
	    	
			k[j][3] = h*fns[j](t[i] + h, x[0][i] + k[0][2], x[1][i] + k[1][2], x[2][i] + k[2][2],
                                 x[3][i] + k[3][2], x[4][i] + k[4][2], x[5][i] + k[5][2]);
	    		
		}
		
		for(int j = 0; j<nf ; j++) {//Ciclo funciónes aproximadas
	    	
			x[j][i + 1] = x[j][i] + (1.0/6.0)*(k[j][0] + 2*k[j][1] + 2*k[j][2] + k[j][3]);
	    	
		}
	      
		if(i%50==0){
			d++;
	    	printf("| %3d |  %10.2f  |  %10.2f  |  %10.2f  |  %10.2f  |  %10.2f  |  %10.2f  |\n",d, x[0][i], x[1][i], x[2][i], x[3][i], x[4][i],x[5][i]);
		
		}
			//printf("| %6.2f |  %10.2f  |  %10.2f  |  %10.2f  |  %10.2f  |  %10.2f  |  %10.2f  |\n",t[i], x[0][i], x[1][i], x[2][i], x[3][i], x[4][i],x[5][i]);
	}
	
cout<<" -----------------------------------------------------------------------------------------------\n"; 
   
}





int main(int argc, char * argv[]) {

	double ci[6]={n,1400.0,140.0,5368.44,3.486,0};
	
	RK4S(ci,0,366,0.02);
	
	
}