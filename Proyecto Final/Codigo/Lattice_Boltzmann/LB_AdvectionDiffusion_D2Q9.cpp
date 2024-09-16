#include <iostream>
#include <cmath>
#include <fstream>
#include <string>  // Add this line to include the <string> header
#include <iomanip> // iomap sirve para setw y setfill
#include <algorithm> // std::fill

//-------------------------------CONSTANTES GLOBASLES------------------------
const int Lx = 10;
const int Ly = 14;

const int Q = 9;// Número de direcciones en el espacio de velocidades

const double C = 1;//
const double Cs = C/sqrt(3);// Velocidad de la onda sonora
const double Cs2 = Cs * Cs;// Velocidad de la onda sonora al cuadrado

const double theta = 1; // Coeficiente de relajación para el término de fuente

//-------------------------------VARIBLES GLOBALES------------------------

// Término de fuente para cada celda
double S_glob; //Termino de fuente para toda la cuadricula
double *S = new double[Lx*Ly]();  // Inicializa todos los valores a 0.0

double D; // Coeficiente de difusión
double tau; // Tiempo de relajación 
double Utau;  // Inverso del tiempo de relajación
double UmUtau; // 1 - 1/tau

//-------------------------------CLASES--------------------------------------
//--------------------- class LatticeBoltzman ------------
class LatticeBoltzman{
private:
    double w[Q];           // Pesos
    int Vx[Q], Vy[Q];      // Vectores de velocidad
    double *f, *fnew;      // Funciones de distribución
public:
    LatticeBoltzman(void);
    ~LatticeBoltzman(void);
    int n(int ix, int iy, int i) { return (ix * Ly + iy) * Q + i;};
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew);
    double Jy(int ix, int iy, bool UseNew);
    double feq(double rho0, double Ux0, double Uy0, int i);
    //double fsource(double rho0, double Ux0, double Uy0, int i, int ind);
    //double dif_fsource(double rho0, double Ux0, double Uy0, int i, int ind, double delta_t);
    void Collision(double delta_t);
    void ImposeFields(int t);
    void Advection(void);
    void Start(double rho0, double Ux0, double Uy0, double mu_x,
               double mu_y, double sigma_x, double sigma_y);
    void Print(std::string NameFile, double t);
    void Printframe(double t);
};

LatticeBoltzman::LatticeBoltzman(void){ 
    // Set the weights
    w[0] = 4.0 / 9;
    w[1] = w[2] = w[3] = w[4] = 1.0 / 9;
    w[5] = w[6] = w[7] = w[8] = 1.0 / 36;
    // Set the velocity vectors
    Vx[8] = 1;  Vx[1] = 1;  Vx[5] = 1;
    Vx[4] = 0;  Vx[0] = 0;  Vx[2] = 0;
    Vx[7] = -1; Vx[3] = -1; Vx[6] = -1;

    Vy[8] = -1; Vy[1] = 0;  Vy[5] = 1;
    Vy[4] = -1; Vy[0] = 0;  Vy[2] = 1;
    Vy[7] = -1; Vy[3] = 0;  Vy[6] = 1;

    // Create the dynamic arrays
    int ArraySize = Lx * Ly * Q;
    f = new double[ArraySize]();
    fnew = new double[ArraySize]();
    //std::cout<<"f: "<<fnew[100]<<std::endl;
}

LatticeBoltzman::~LatticeBoltzman(void){
    delete[] f;
    delete[] fnew;
}

//---------------------Campos macroscópicos---------------------
double LatticeBoltzman::rho(int ix, int iy, bool UseNew){
    double sum = 0;
    for (int i = 0; i < Q; i++){
        int n0 = n(ix, iy, i);
        if (UseNew) sum += fnew[n0];
        else sum += f[n0];
    }  
    //std::cout<<"rho: "<<sum<<std::endl;
    return sum;
}

double LatticeBoltzman::Jx(int ix, int iy, bool UseNew){
    double sum = 0;
    for (int i = 0; i < Q; i++){
        int n0 = n(ix, iy, i);

        if (UseNew) sum += Vx[i] * fnew[n0];
        else sum += Vx[i] * f[n0];
    }
    return sum;
}

double LatticeBoltzman::Jy(int ix, int iy, bool UseNew){
    double sum = 0;
    for (int i = 0; i < Q; i++){
        int n0 = n(ix, iy, i);

        if (UseNew) sum += Vy[i] * fnew[n0];
        else sum += Vy[i] * f[n0];
    }
    return sum;
}

//---------------------Función de equilibrio---------------------
double LatticeBoltzman::feq(double rho0, double Ux0, double Uy0, int i){
    double UdotVi = Ux0 * Vx[i] + Uy0 * Vy[i];
    double U2 = Ux0 * Ux0 + Uy0 * Uy0;
    double result = rho0 * w[i] * (1 + UdotVi / Cs2 + (UdotVi * UdotVi) / (2 * Cs2 * Cs2) - U2 / (2 * Cs2));

    return result;
}
/*
//-------------------Función de fuente-------------------
double LatticeBoltzman::fsource(double rho0, double Ux0, double Uy0, int i, int ind){
    double UdotVi = Ux0 * Vx[i] + Uy0 * Vy[i];
    double result = w[i]*S[ind]*(1 + ((tau - 0.5)*UdotVi)/((tau - theta*0.5)*Cs2));
    //std::cout<<"S: "<<S[ind]<<std::endl;

    return result;
}

double LatticeBoltzman::dif_fsource(double rho0, double Ux0, double Uy0, int i, int ind, double delta_t){
    double UdotVi = Ux0 * Vx[i] + Uy0 * Vy[i];
    double DtUdotVi = Ux0*(Vx[i] - Vx[i]*delta_t) + Uy0*(Vy[i] - Vy[i]*delta_t);

    double result = S[ind]*w[i]*(1 + (tau - 0.5)/((tau - theta*0.5)*Cs2)*(UdotVi-delta_t*DtUdotVi))/delta_t;
    return result;
}
*/

// ---------------------Evolución temporal---------------------
void LatticeBoltzman::Start(double rho0, double Ux0, double Uy0,
                            double mu_x, double mu_y, double sigma_x, double sigma_y){
    int ix, iy, i, n0;

    for (ix = 0; ix < Lx; ix++){
        for (iy = 0; iy < Ly; iy++){

            double gauss_x = exp(-0.5 * pow((ix - mu_x) / sigma_x, 2)) / (sigma_x * sqrt(2 * M_PI));
            double gauss_y = exp(-0.5 * pow((iy - mu_y) / sigma_y, 2)) / (sigma_y * sqrt(2 * M_PI));
            double rho = rho0; //* gauss_x * gauss_y;

            for (i = 0; i < Q; i++){
                n0 = n(ix, iy, i);
                f[n0] = feq(rho, Ux0, Uy0, i);
                //std::cout<<"f: "<<f[n0]<<std::endl;
            }
        }
    }
}

void LatticeBoltzman::Collision(double delta_t){
    int ix, iy, i, n0;
    double rho0, Ux0, Uy0, source, dif_source;
    for (ix = 0; ix < Lx; ix++){
        for (iy = 0; iy < Ly; iy++){
            rho0 = rho(ix, iy, false);
            Ux0 = Jx(ix, iy, false) / rho0;
            Uy0 = Jy(ix, iy, false) / rho0;
            int ind = (ix * Ly + iy); //Linealizacion de la cuadrícula
            for (i = 0; i < Q; i++){
                n0 = n(ix, iy, i);
                //source = delta_t*fsource(rho0, Ux0, Uy0, i, ind);
                //dif_source = delta_t*delta_t*dif_fsource(rho0, Ux0, Uy0, i, ind, delta_t)*0.5; 
                fnew[n0] = UmUtau * f[n0] + Utau * feq(rho0, Ux0, Uy0, i);
            }
        }
    }
}

void LatticeBoltzman::ImposeFields(int t){
    // Implementación para imponer un campo de velocidad constante
    double rho0, Ux0, Uy0;
    int n0;
    double rho_incendio;
    for (int ix = 0; ix < Lx; ix++){
        for (int iy = 0; iy < Ly; iy++){
            rho0 = rho(ix, iy, true); // Usar fnew para obtener la densidad
            Ux0 = 0.3;                // Velocidad en x
            Uy0 = 0.3;                // Velocidad en y
            for (int i = 0; i < Q; i++){
                n0 = n(ix, iy, i);
                fnew[n0] = feq(rho0, Ux0, Uy0, i);
            }
            if(t<=9){
                if(ix==0 && iy==4){ //Usme(0,4)
                    for (int i = 0; i < Q; i++){
                        //std::cout<<"Incendio en: "<<ix<<", "<<iy<<std::endl;
                        rho_incendio = 1;
                        n0 = n(ix, iy, i);
                        fnew[n0] = feq(rho_incendio, 0, 0, i);
                    }
                }
                if(ix==6 && iy==7){ //Quebrada la vieja(6,7)
                    for (int i = 0; i < Q; i++){
                        //std::cout<<"Incendio en: "<<ix<<", "<<iy<<std::endl;
                        rho_incendio = 1;
                        n0 = n(ix, iy, i);
                        fnew[n0] = feq(rho_incendio, 0, 0, i);
                    }
                }
            }
            
                /*
            if(ix == 0 || ix == Lx-1 || iy == 0 || iy == Ly-1){
                for (int i = 0; i < Q; i++){
                    int n0 = n(ix, iy, i);
                    fnew[n0] = feq(1e7, 0, 0, i);
                }
            }*/
        }   
    }
}

void LatticeBoltzman::Advection(void){
    int ix, iy, i, ixnext, iynext, n0, n0next;

    for (ix = 0; ix < Lx; ix++){ // for each cell
        for (iy = 0; iy < Ly; iy++){
            for (i = 0; i < Q; i++){ // on each direction
                ixnext = (ix + Vx[i] + Lx) % Lx;
                iynext = (iy + Vy[i] + Ly) % Ly;
                n0 = n(ix, iy, i);
                n0next = n(ixnext, iynext, i);
                f[n0next] = fnew[n0]; // periodic boundaries
            }
        }    
    }    
}

//---------------------Impresión de resultados---------------------
void LatticeBoltzman::Print(std::string NameFile, double t){
    std::ofstream MyFile(NameFile);
    double rho0, Ux0, Uy0;

    for (int ix = 0; ix < Lx; ix++){
        for (int iy = 0; iy < Ly; iy++){
            rho0 = rho(ix, iy, false);
            MyFile << ix << " " << iy << " " << rho0 << std::endl;
        }
        MyFile << std::endl;
    }
    MyFile.close();
}

void LatticeBoltzman::Printframe(double t){
    // Check if the "frames" directory exists, if not, create it
    int ret;
    if (system("test -d frames") != 0){
        ret = system("mkdir frames");
        if (ret != 0) std::cerr << "Error: No se pudo crear el directorio frames." << std::endl;
    }
    std::ofstream GnuplotScript("frame_script.gp");
    GnuplotScript << "set terminal pngcairo size 800,800 enhanced font 'Verdana,10'" << std::endl;
    GnuplotScript << "set output 'frames/density_" << std::setw(3) << std::setfill('0') << t << ".png'" << std::endl;
    GnuplotScript << "set pm3d map" << std::endl;
    GnuplotScript << "set size ratio -1" << std::endl;
    GnuplotScript << "set xrange [0:" << Lx << "]" << std::endl;
    GnuplotScript << "set yrange [0:" << Ly << "]" << std::endl;
    GnuplotScript << "set cbrange [0:*]" << std::endl;
    GnuplotScript << "set palette defined (0 'black', 1 'red', 2 'yellow', 3 'white', 4 'red')" << std::endl;
    GnuplotScript << "set title 'Densidad en t = " << t << "'" << std::endl;
    GnuplotScript << "plot 'data/density_" << std::setw(3) << std::setfill('0') << t << ".dat' u 1:2:3 w image" << std::endl;
    GnuplotScript.close();

    // Ejecutar el script de Gnuplot
    ret = system("gnuplot frame_script.gp");
    if (ret != 0) std::cerr << "Error: No se pudo ejecutar el script de gnuplot" << std::endl;
}

int main(int argc, char* argv[]){
    // Parámetros de la simulación
    int tframe = 1, tmax = 20, delta_t = 1, ret;
    double rho0 = 0.001, Ux0 = 0.3, Uy0 = 0.3;                               // Densidad inicial y velocidad
    double mu_x = Lx / 2, mu_y = Ly / 2, sigma_x = Lx/4, sigma_y = Ly/4; // Parámetros de la distribución gaussiana
    LatticeBoltzman Air;

    // Parámetros a ajustar
    S_glob = std::stod(argv[1]); //Termino de fuente global
    D = std::stod(argv[2]); //Coeficiente de difusión 0.016 -> tau = 0.548

    std::fill(S, S + Lx*Ly, S_glob); // Inicializa todos los valores de S 

    tau = (D / delta_t*Cs2) + 0.5; // Calculo tiempo de relajación en función de D
    Utau = 1.0 / tau;    UmUtau = 1 - Utau; //Otros valores útiles 

    // Iniciar la simulación
    Air.Start(rho0, Ux0, Uy0, mu_x, mu_y, sigma_x, sigma_y);

    // Ejecutar la simulación
    for (int t = 0; t <= tmax; t++)
    {
        
        Air.Collision(delta_t);
        Air.ImposeFields(t); // Ux0, Uy0);
        Air.Advection();
        if (t % tframe == 0){

            if (system("test -d data") != 0){
                ret = system("mkdir data");
                if (ret != 0) std::cerr << "Error: No se pudo crear el directorio data." << std::endl;
            }
            // Imprimir los datos en archivo
            std::stringstream ss;
            ss << "data/density_" << std::setw(3) << std::setfill('0') << t << ".dat";
            Air.Print(ss.str(), t);
            // Generar y guardar el frame con Gnuplot
            Air.Printframe(t);
            std::cout << "Porcentaje de avance: " << (t * 100) / tmax << "%" << std::endl;
        }
    }
    delete[] S;
    return 0;
}
