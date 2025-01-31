#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//------------------- CONSTANTES GLOBALES -------------
const int Lx = 400;
const int Ly = 200;

const int Q = 5;
const double W0 = 1.0 / 3;

const double tau = 0.5;
const double Utau = 1.0 / tau;
const double UmUtau = 1 - Utau;
const double theta = 20 * M_PI / 180.0;

//------------------- FUNCIONES GLOBALES -------------
class LatticeBoltzmann {
private:
    double w[Q];      // Weights
    int Vx[Q], Vy[Q]; // Velocity vectors
    double *f, *fnew; // Distribution Functions
public:
    //----Constructor y destructor----
    LatticeBoltzmann(void);
    ~LatticeBoltzmann(void);
    int n(int ix, int iy, int i) { return (ix * Ly + iy) * Q + i; };
    //----Campos macroscopicos----
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew);
    double Jy(int ix, int iy, bool UseNew);
    //----Funciones de equilibrio----
    double feq(double rho0, double Jx0, double Jy0, int i, double C);
    //----Evolucion temporal----
    void Start(double rho0, double Jx0, double Jy0);
    void Collision(void);
    void ImposeFields(int t);
    void Advection(void);
    void Print(const char *NameFile);
    double Ccelda(int ix, int iy); // Nueva función para calcular Ccelda
};

LatticeBoltzmann::LatticeBoltzmann(void) {
    // Set the weights
    w[0] = W0; w[1] = w[2] = w[3] = w[4] = (1.0 - W0) / 4;
    // Set the velocity vectors
    Vx[0] = 0;  Vx[1] = 1;  Vx[2] = 0;  Vx[3] = -1; Vx[4] = 0;
    Vy[0] = 0;  Vy[1] = 0;  Vy[2] = 1;  Vy[3] = 0;  Vy[4] = -1;
    // Create the dynamic arrays
    int ArraySize = Lx * Ly * Q;
    f = new double[ArraySize];  fnew = new double[ArraySize];
}

LatticeBoltzmann::~LatticeBoltzmann(void) {
    delete[] f;  delete[] fnew;
}

double LatticeBoltzmann::rho(int ix, int iy, bool UseNew) {
    double sum; int i, n0;
    for (sum = 0, i = 0; i < Q; i++) {
        n0 = n(ix, iy, i);
        if (UseNew) sum += fnew[n0]; else sum += f[n0];
    }
    return sum;
}

double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew) {
    double sum; int i, n0;
    for (sum = 0, i = 0; i < Q; i++) {
        n0 = n(ix, iy, i);
        if (UseNew) sum += Vx[i] * fnew[n0]; else sum += Vx[i] * f[n0];
    }
    return sum;
}

double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew) {
    double sum; int i, n0;
    for (sum = 0, i = 0; i < Q; i++) {
        n0 = n(ix, iy, i);
        if (UseNew) sum += Vy[i] * fnew[n0]; else sum += Vy[i] * f[n0];
    }
    return sum;
}

double LatticeBoltzmann::feq(double rho0, double Jx0, double Jy0, int i, double C) {
    double C2 = C * C;
    double AUX0 = 1 - 3 * C2 * (1 - W0);
    if (i > 0)
        return 3 * w[i] * (C2 * rho0 + Vx[i] * Jx0 + Vy[i] * Jy0);
    else
        return rho0 * AUX0;
}


double LatticeBoltzmann::Ccelda(int ix, int iy) {
    double n1 = 1.0;
    double n2 = 2.0;
    //double ix0 = 100.0;
    double ix0 = 100.0 - (iy*tan(theta)) ;
    double ancho = 5.0; // Ancho de transición para suavizar la función

    double n = 1.5 + tanh((ix - (ix0)) / ancho)/ 2.0;
    double C = 0.707 / n;
    return C;
}


void LatticeBoltzmann::Start(double rho0, double Jx0, double Jy0) {
    int ix, iy, i, n0;
    for (ix = 0; ix < Lx; ix++) // for each cell
        for (iy = 0; iy < Ly; iy++)
            for (i = 0; i < Q; i++) { // on each direction
                n0 = n(ix, iy, i);
                f[n0] = feq(rho0, Jx0, Jy0, i, Ccelda(ix, iy));
            }
}

void LatticeBoltzmann::Collision(void) {
    int ix, iy, i, n0; double rho0, Jx0, Jy0;
    for (ix = 0; ix < Lx; ix++) // for each cell
        for (iy = 0; iy < Ly; iy++) {
            // compute the macroscopic fields on the cell
            rho0 = rho(ix, iy, false); Jx0 = Jx(ix, iy, false); Jy0 = Jy(ix, iy, false);
            for (i = 0; i < Q; i++) { // for each velocity vector
                n0 = n(ix, iy, i);
                fnew[n0] = UmUtau * f[n0] + Utau * feq(rho0, Jx0, Jy0, i, Ccelda(ix, iy));
            }
        }
}

void LatticeBoltzmann::ImposeFields(int t) {
    int i, ix, iy, n0;
    double lambda, omega, rho0, Jx0, Jy0; lambda = 10;
    // an oscillating source in the middle
    for (iy = 0; iy < Ly; iy++) {
        ix = 0; // iy = Ly / 2;
        omega = 2 * M_PI / lambda * Ccelda(ix, iy);
        rho0 = 10 * sin(omega * t); Jx0 = Jx(ix, iy, false); Jy0 = Jy(ix, iy, false);
        for (i = 0; i < Q; i++) {
            n0 = n(ix, iy, i);
            fnew[n0] = feq(rho0, Jx0, Jy0, i, Ccelda(ix, iy));
        }
    }
}

void LatticeBoltzmann::Advection(void) {
    int ix, iy, i, ixnext, iynext, n0, n0next;
    for (ix = 0; ix < Lx; ix++) // for each cell
        for (iy = 0; iy < Ly; iy++)
            for (i = 0; i < Q; i++) { // on each direction
                ixnext = (ix + Vx[i] + Lx) % Lx; iynext = (iy + Vy[i] + Ly) % Ly;
                n0 = n(ix, iy, i); n0next = n(ixnext, iynext, i);
                f[n0next] = fnew[n0]; // periodic boundaries
            }
}

void LatticeBoltzmann::Print(const char *NameFile) {
    ofstream MyFile(NameFile); double rho0; int ix, iy,Jx0, Jy0,i,n0;
    for (ix = 0; ix < Lx / 2; ix++) {
        for (iy = 0; iy < Ly; iy++) {
          rho0 = rho(ix, iy, true);
          Jx0=Jx(ix, iy, false);
          Jy0=Jy(ix, iy, false);
          for (i = 0; i < Q; i++) {
            n0 = n(ix, iy, i);

            MyFile << ix << " " << iy << " " << rho0 <<endl;//" " << Ccelda(ix,iy)  << " "  <<Jx0<<" "<< Jy0<<" "<<feq(rho0, Jx0, Jy0, i, Ccelda(ix, iy))<< endl;
          }
        }
        MyFile << endl;
    }
    MyFile.close();
}

//--------------- FUNCIONES GLOBALES ------------


int main(void) {
    LatticeBoltzmann Ondas;
    int t, tmax = 400;
    double rho0 = 0, Jx0 = 0, Jy0 = 0;

    // Start
    Ondas.Start(rho0, Jx0, Jy0);
    // Run
    for (t = 0; t < tmax; t++) {
        Ondas.Collision();
        Ondas.ImposeFields(t);
        Ondas.Advection();
    }
    // Show
    Ondas.Print("Ondas.dat");

    return 0;
}
