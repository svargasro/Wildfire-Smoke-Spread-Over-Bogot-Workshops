#include <iostream>
#include <cmath>
#include <fstream>
#include <string>    // Para el manejo de strings
#include <iomanip>   // Para setw y setfill
#include <algorithm> // Para std::fill

//-------------------------------CONSTANTES GLOBALES------------------------
// Dimensiones de la cuadrícula
// const int Lx = 500; // 100
const int Lx = 100; // 100
const int Ly = Lx*(1.4); // 140


int t_hour = 240;
double *Ux = new double[Lx * Ly * t_hour];
double *Uy = new double[Lx * Ly * t_hour]; // Velocidades en x y y

// Número de direcciones en el espacio de velocidades
const int Q = 9;

// Velocidad del lattice y velocidad del sonido
const double C = 1;            // Velocidad característica del lattice
const double Cs = C / sqrt(3); // Velocidad de la onda sonora
const double Cs2 = Cs * Cs;    // Velocidad de la onda sonora al cuadrado

//-------------------------------VARIABLES GLOBALES------------------------

// Parámetros de difusión y relajación
double D = 0.016; // Coeficiente de difusión
double tau;       // Tiempo de relajación
double Utau;      // Inverso del tiempo de relajación
double UmUtau;    // 1 - 1/tau

//-------------------------------CLASES--------------------------------------
// Clase para el método de Lattice Boltzmann
class LatticeBoltzman
{
private:
    double w[Q];      // Pesos de las direcciones
    int Vx[Q], Vy[Q]; // Vectores de velocidad en las direcciones x e y
    double *f, *fnew; // Funciones de distribución para los diferentes estados (previo, actual y nuevo)

public:
    LatticeBoltzman(void);                                                                                     // Constructor
    ~LatticeBoltzman(void);                                                                                    // Destructor
    int n(int ix, int iy, int i) { return (ix * Ly + iy) * Q + i; };                                           // Conversión de índices 2D a 1D
    double rho(int ix, int iy, bool UseNew);                                                                   // Cálculo de la densidad
    double Jx(int ix, int iy, bool UseNew);                                                                    // Cálculo del flujo en x
    double Jy(int ix, int iy, bool UseNew);                                                                    // Cálculo del flujo en y
    double feq(double rho0, double Ux0, double Uy0, int i);                                                    // Función de equilibrio
    double fsource(double rho0, double Ux0, double Uy0, int i, int ind);                                       // Término de fuente
    void Collision();                                                                                          // Fase de colisión
    void ImposeFields(double Ux0, double Uy0, int t_hour);                                                     // Imponer campos de velocidad
    void Advection(void);                                                                                      // Fase de advección
    void Start(double rho0, double Ux0, double Uy0, double mu_x, double mu_y, double sigma_x, double sigma_y); // Inicialización
    void PrintData(std::string NameFile, double t);                                                            // Imprimir resultados a archivo
    void PrintFrame(double t);                                                                                 // Generar y guardar el frame de la simulación
};

//-------------------------------IMPLEMENTACIÓN DE LA CLASE LatticeBoltzman------------------------

// Constructor: Inicializa pesos, velocidades y distribuciones
LatticeBoltzman::LatticeBoltzman(void)
{
    // Asignación de los pesos para cada dirección
    w[0] = 4.0 / 9;
    w[1] = w[2] = w[3] = w[4] = 1.0 / 9;
    w[5] = w[6] = w[7] = w[8] = 1.0 / 36;

    // Asignación de los vectores de velocidad
    Vx[8] = 1;
    Vx[1] = 1;
    Vx[5] = 1;
    Vx[4] = 0;
    Vx[0] = 0;
    Vx[2] = 0;
    Vx[7] = -1;
    Vx[3] = -1;
    Vx[6] = -1;

    Vy[8] = -1;
    Vy[1] = 0;
    Vy[5] = 1;
    Vy[4] = -1;
    Vy[0] = 0;
    Vy[2] = 1;
    Vy[7] = -1;
    Vy[3] = 0;
    Vy[6] = 1;

    // Creación de los arreglos dinámicos para las funciones de distribución
    int ArraySize = Lx * Ly * Q;

    f = new double[ArraySize]();
    fnew = new double[ArraySize]();
}

// Destructor: libera la memoria dinámica
LatticeBoltzman::~LatticeBoltzman(void)
{
    delete[] f;
    delete[] fnew;
}

// Función para calcular la densidad macroscópica (rho)
//---------------------Campos macroscópicos---------------------
double LatticeBoltzman::rho(int ix, int iy, bool UseNew)
{
    double sum = 0;
    for (int i = 0; i < Q; i++)
    {
        int n0 = n(ix, iy, i);

        if (UseNew)
            sum += fnew[n0];
        else
            sum += f[n0];
    }
    // std::cout<<"rho: "<<sum<<std::endl;
    return sum;
}
// Función para calcular el flujo en la dirección x
double LatticeBoltzman::Jx(int ix, int iy, bool UseNew)
{
    double sum = 0;
    for (int i = 0; i < Q; i++)
    {
        int n0 = n(ix, iy, i);

        if (UseNew)
            sum += Vx[i] * fnew[n0];
        else
            sum += Vx[i] * f[n0];
    }
    return sum;
}
// Función para calcular el flujo en la dirección y
double LatticeBoltzman::Jy(int ix, int iy, bool UseNew)
{
    double sum = 0;
    for (int i = 0; i < Q; i++)
    {
        int n0 = n(ix, iy, i);

        if (UseNew)
            sum += Vy[i] * fnew[n0];
        else
            sum += Vy[i] * f[n0];
    }
    return sum;
}

// Función de equilibrio para la colisión
double LatticeBoltzman::feq(double rho0, double Ux0, double Uy0, int i)
{
    double UdotVi = Ux0 * Vx[i] + Uy0 * Vy[i];
    double U2 = Ux0 * Ux0 + Uy0 * Uy0;
    double result = rho0 * w[i] * (1 + UdotVi / Cs2 + (UdotVi * UdotVi) / (2 * Cs2 * Cs2) - U2 / (2 * Cs2));
    return result;
}

// ---------------------Evolución temporal---------------------

// Inicialización de la función de distribución y condiciones iniciales
void LatticeBoltzman::Start(double rho0, double Ux0, double Uy0,
                            double mu_x, double mu_y, double sigma_x, double sigma_y)
{
    int ix, iy, i, n0;

    // Recorre toda la cuadrícula de celdas
    for (ix = 0; ix < Lx; ix++) // Para cada celda en el eje x
    {
        for (iy = 0; iy < Ly; iy++) // Para cada celda en el eje y
        {
            // Calcular el valor de la función gaussiana en las direcciones x e y
            double gauss_x = exp(-0.5 * pow((ix - mu_x) / sigma_x, 2)) / (sigma_x * sqrt(2 * M_PI));
            double gauss_y = exp(-0.5 * pow((iy - mu_y) / sigma_y, 2)) / (sigma_y * sqrt(2 * M_PI));
            double rho = rho0; //* gauss_x * gauss_y;

            // Inicializar las funciones de distribución en todas las direcciones de velocidad
            for (i = 0; i < Q; i++) // Para cada dirección de la función de distribución
            {
                n0 = n(ix, iy, i); // Convertir los índices 2D a 1D
                // Establecer la función de distribución de equilibrio como valor inicial
                f[n0] = feq(rho, Ux0, Uy0, i); // Inicializa f
            }
        }
    }
}

// Fase de colisión: actualiza las funciones de distribución
void LatticeBoltzman::Collision()
{
    int ix, iy, i, n0;
    double rho0, Ux0, Uy0;

    for (ix = 0; ix < Lx; ix++)
    {
        for (iy = 0; iy < Ly; iy++)
        {
            rho0 = rho(ix, iy, false);      // Calcula la densidad actual
            Ux0 = Jx(ix, iy, false) / rho0; // Flujo en x
            Uy0 = Jy(ix, iy, false) / rho0; // Flujo en y

            for (i = 0; i < Q; i++)
            {
                n0 = n(ix, iy, i);
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

// Fase de advección: actualiza las funciones de distribución con condiciones periódicas
void LatticeBoltzman::Advection(void)
{
    int ix, iy, i, ixnext, iynext, n0, n0next;

    for (ix = 0; ix < Lx; ix++) // Itera sobre la cuadrícula
    {
        for (iy = 0; iy < Ly; iy++)
        {
            for (i = 0; i < Q; i++) // Sobre cada dirección de velocidad
            {
                ixnext = (ix + Vx[i] + Lx) % Lx; // Condiciones periódicas en x
                iynext = (iy + Vy[i] + Ly) % Ly; // Condiciones periódicas en y
                n0 = n(ix, iy, i);               // Índice actual
                n0next = n(ixnext, iynext, i);   // Índice siguiente

                f[n0next] = fnew[n0]; // Actualiza el estado actual
            }
        }
    }
}

//----------------------Carga de resultados----------------------
// Funcion para imprimir los resultados de la simulación en un archivo
void LoadData(std::string NameFile)
{
    // Abre un archivo de entrada para cargar los datos de la simulación
    std::ifstream MyFile(NameFile);
    std::ofstream MyFile2("Debug.txt");

    for (int ix = 0; ix < Lx; ix++)
    {
        for (int iy = 0; iy < Ly; iy++)
        {
            for (int t = 0; t < t_hour; t++)
            {
                int index = (iy * Lx + ix) * t_hour + t;
                // int index = ix * Ly * t_hour + iy * t_hour + t;
                MyFile >> Ux[index] >> Uy[index];
                MyFile2 << index + 1 << " " << Ux[index] << " " << Uy[index] << std::endl;
            }
        }
    }
    std::cout << "Datos cargados exitosamente" << std::endl;
    // Cierra el archivo de entrada
    MyFile.close();
    MyFile2.close();
}
// Funciones de impresión de resultados y visualización
//---------------------Impresión de resultados---------------------
// Funcion para imprimir los resultados de la simulación en un archivo
void LatticeBoltzman::PrintData(std::string NameFile, double t)
{
    // Abre un archivo de salida para guardar los datos de la simulación
    std::ofstream MyFile(NameFile);
    double rho0, Ux0, Uy0;

    // Recorre la cuadrícula en pasos de 4 para ahorrar espacio en los archivos de salida
    for (int ix = 0; ix < Lx; ix += 4)
    {
        for (int iy = 0; iy < Ly; iy += 4)
        {
            // Calcula la densidad (rho0) y las velocidades (Ux0, Uy0) para cada celda
            rho0 = rho(ix, iy, false);      // Densidad en la celda (ix, iy)
            Ux0 = Jx(ix, iy, false) / rho0; // Velocidad en x
            Uy0 = Jy(ix, iy, false) / rho0; // Velocidad en y

            // Escribe los datos en el archivo de salida
            MyFile << ix << " " << iy << " "std::endl; //<< rho0 << " " << Ux0 << " " << Uy0 << std::endl;
        }
        MyFile << std::endl; // Inserta un salto de línea después de cada fila de la cuadrícula
    }

    // Cierra el archivo de salida
    MyFile.close();
}

// Funcion para imprimir un frame de la simulación
void LatticeBoltzman::PrintFrame(double t)
{
    // Verifica si el directorio "frames" existe, de lo contrario, lo crea
    int ret;
    if (system("test -d frames") != 0) // Comando del sistema para verificar si el directorio existe
    {
        ret = system("mkdir frames"); // Comando para crear el directorio
        if (ret != 0)
            std::cerr << "Error: No se pudo crear el directorio frames." << std::endl;
    }

    // Crear el archivo de script para Gnuplot
    std::ofstream GnuplotScript("frame_script.gp");

    GnuplotScript << "set terminal pngcairo size 800,800 enhanced font 'Verdana,10'" << std::endl;                    // Configuración del terminal de salida para Gnuplot, generando imágenes PNG
    GnuplotScript << "set output 'frames/density_" << std::setw(3) << std::setfill('0') << t << ".png'" << std::endl; // Nombre del archivo de salida PNG, numerado según el tiempo de simulación
    GnuplotScript << "set pm3d map" << std::endl;                                                                     // Usar un mapa de calor (pm3d) para visualizar los datos
    GnuplotScript << "set size ratio -1" << std::endl;                                                                // Mantener una relación de aspecto cuadrada
    GnuplotScript << "set xrange [0:" << Lx << "]" << std::endl;                                                      // Definir los límites de los ejes x e y, basados en las dimensiones de la cuadrícula
    GnuplotScript << "set yrange [0:" << Ly << "]" << std::endl;
    GnuplotScript << "set cbrange [0:*]" << std::endl;                                                                        // Configurar la escala de colores (color bar) para los valores de densidad
    GnuplotScript << "set palette defined (0 'black', 1 'red', 2 'orange', 3 'yellow', 4 'white')" << std::endl;              // Definir la paleta de colores, de negro a blanco pasando por rojo, naranja y amarillo
    GnuplotScript << "set title 'Densidad en t = " << t << "'" << std::endl;                                                  // Definir el título del gráfico, basado en el tiempo de simulación actual
    GnuplotScript << "plot 'data/density_" << std::setw(3) << std::setfill('0') << t << ".dat' u 1:2:3 w image" << std::endl; // Instrucción para graficar los datos de densidad desde el archivo correspondiente
    GnuplotScript.close();                                                                                                    // Cerrar el archivo de script

    // Ejecutar el script de Gnuplot para generar la imagen PNG
    ret = system("gnuplot frame_script.gp");
    if (ret != 0)
        std::cerr << "Error: No se pudo ejecutar el script de gnuplot" << std::endl;
}

// Funcion main
int main(int argc, char *argv[])
{

    // Parámetros generales de la simulación
    int t, tframe = 1000, tmax = t_hour * 100, delta_t = 1, ret; // tframe: intervalo entre frames, tmax: tiempo máximo de simulación
    // int t, tframe = 100, tmax = 1000, delta_t = 1, ret;              // tframe: intervalo entre frames, tmax: tiempo máximo de simulación
    // double rho0 = 100000000000, Ux0 = Lx / 100, Uy0 = Ly / 100;             // Densidad inicial y velocidades iniciales en x e y
    double rho0 = 0.001, Ux0 = 0.01, Uy0 = 0.01;  

    // Parámetros para la distribución gaussiana que inicializa la densidad
    double mu_x = Lx / 2.0, mu_y = Ly / 2.0, sigma_x = Lx / 32.0, sigma_y = Ly / 32.0; // Parámetros de la gaussiana: centro (mu_x, mu_y), sigma: ancho

    // std::cout << mu_x << " " << mu_y << " " << sigma_x << " " << sigma_y << std::endl;

    // Crear una instancia de la clase LatticeBoltzman
    LatticeBoltzman Air;

    // LoadData("velocity.txt");

    // Leer parámetros desde la línea de comandos: término de fuente y coeficiente de difusión
    D = std::stod(argv[1]); // Coeficiente de difusión (ejemplo: 0.016 -> tau = 0.548)

    // Calcular el tiempo de relajación tau basado en el coeficiente de difusión y el paso de tiempo
    tau = (D / delta_t * Cs2) + 0.5;

    // Calcular otros valores útiles basados en tau
    Utau = 1.0 / tau;
    UmUtau = 1 - Utau; // 1 - 1/tau

    // Inicializar la simulación con las condiciones iniciales (densidad y velocidades)
    Air.Start(rho0, Ux0, Uy0, mu_x, mu_y, sigma_x, sigma_y);

    // Bucle principal de la simulación
    for (t = 0; t <= tmax; t++)
    {
        // Etapas de la simulación: colisión, imposición de campos y advección
        Air.Collision();
        Air.ImposeFields(Ux0, Uy0, t);
        Air.Advection();

        // Guardar resultados cada tframe pasos
        if (t % tframe == 0)
        {
            // Verificar si el directorio "data" existe, de lo contrario, crearlo
            if (system("test -d data") != 0)
            {
                ret = system("mkdir data");
                if (ret != 0)
                    std::cerr << "Error: No se pudo crear el directorio data." << std::endl;
            }

            // Crear un archivo para guardar los resultados de la simulación
            std::stringstream ss;
            ss << "data/density_" << std::setw(3) << std::setfill('0') << t << ".dat"; // Nombre del archivo basado en el tiempo de simulación
            Air.PrintData(ss.str(), t);                                                // Llamada a la función que guarda los datos

            // Generar un frame de la simulación usando Gnuplot
            Air.PrintFrame(t);

            // Mostrar el porcentaje de avance de la simulación en la consola
            std::cout << "Porcentaje de avance: " << (t * 100) / tmax << "%" << std::endl;
        }
    }

    return 0;
}
