#include <iostream>
#include <cmath>
#include <fstream>
#include <string>    // Para el manejo de strings
#include <iomanip>   // Para setw y setfill
#include <algorithm> // Para std::fill
#include <vector>

//-------------------------------CONSTANTES GLOBALES------------------------
// Dimensiones de la cuadrícula
const int Lx = 100; // 100
const int Ly = Lx*(1.4); // 140
int iter_per_hour = 9604;
//int iter_per_hour = 9604;

const int Q = 9;// Número de direcciones en el espacio de velocidades

int t_hour = 241;
//int t_hour = 240;
double *Ux = new double[Lx * Ly * t_hour];
double *Uy = new double[Lx * Ly * t_hour]; // Velocidades en x y y

//Variables para añadir fuentes del incendio
std::vector<int> id;
std::vector<int> rho_f;
int index_f = 0;
bool IsData = true;

const double C = 1;//
const double Cs = C/sqrt(3);// Velocidad de la onda sonora
const double Cs2 = Cs * Cs;// Velocidad de la onda sonora al cuadrado

//-------------------------------VARIBLES GLOBALES------------------------

double D = 0.016; // Coeficiente de difusión
double tau; // Tiempo de relajación 
double Utau;  // Inverso del tiempo de relajación
double UmUtau; // 1 - 1/tau

//-------------------------------CLASES--------------------------------------
//--------------------- class LatticeBoltzman ------------
class LatticeBoltzman{
    private:
        double w[Q];      // Pesos de las direcciones
        int Vx[Q], Vy[Q]; // Vectores de velocidad en las direcciones x e y
        double *f, *fnew; // Funciones de distribución para los diferentes estados (previo, actual y nuevo)
    public:
    LatticeBoltzman(void);
    ~LatticeBoltzman(void);
    int n(int ix, int iy, int i) { return (ix * Ly + iy) * Q + i;};
    double rho(int ix, int iy, bool UseNew);
    double Jx(int ix, int iy, bool UseNew);
    double Jy(int ix, int iy, bool UseNew);
    double feq(double rho0, double Ux0, double Uy0, int i);
    void Collision();
    void ImposeFields(int t);
    void Advection(void);
    void Start(double rho0, double Ux0, double Uy0);
    void PrintData(std::string NameFile);                                                            // Imprimir resultados a archivo
    void PrintFrame(double t, int t_a);                                                                                 // Generar y guardar el frame de la simulación
};

//-------------------------------FUNCIONES GLOBALES------------------------
//Declaración de funciones
void LoadData(std::string NameFile);                 // Cargar los datos de velocidad desde un archivo
void LoadDataBinary(const std::string &filename, double *Ux, double *Uy, int size); // Cargar los datos desde un archivo binario
void SaveDataBinary(const std::string &filename, double *Ux, double *Uy, int size); // Guardar los datos en un archivo binario
void LoadDataFires(std::string NameFile);


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
    double result = rho0 * w[i] * (1 + UdotVi / Cs2 + (UdotVi * UdotVi) / (2.0 * Cs2 * Cs2) - U2 / (2.0 * Cs2));

    return result;
}

// ---------------------Evolución temporal---------------------

// Inicialización de la función de distribución y condiciones iniciales
void LatticeBoltzman::Start(double rho0, double Ux0, double Uy0)
{
    int ix, iy, i, n0;

    // Recorre toda la cuadrícula de celdas
    for (ix = 0; ix < Lx; ix++) // Para cada celda en el eje x
    {
        for (iy = 0; iy < Ly; iy++) // Para cada celda en el eje y
        {
            double rho = rho0;

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

void LatticeBoltzman::ImposeFields(int t)
{
    double rho0;
    int auxT=(t+1)/iter_per_hour;
    double Ux0,Uy0;
    int index_tmp = 0;
    for (int ix = 0; ix < Lx; ix++)
    {
        for (int iy = 0; iy < Ly; iy++)
                    // int index = ix*Ly+iy +Lx*Ly*auxT;
            // if((t+1)%iter_per_hour == 0)  index_f = index_f + fire_per_hour;
            // Ux0 = Ux[index];
            // Uy0 = Uy[index];
            // rho0 = rho(ix, iy, true); // Usar fnew para obtener la densidad
            // if(ix*Ly+iy == id[index_f]){
            //     rho0 = rho_f[index_f];
            //     index_f++;
            // }
        {
            int index = ix*Ly+iy +Lx*Ly*auxT;
            int index_tmp =0;
            Ux0 = Ux[index];
            Uy0 = Uy[index];
            rho0 = rho(ix, iy, true); // Usar fnew para obtener la densidad
            for (int i = 0; i < Q; i++)
            {
                int n0 = n(ix, iy, i);
                if(((ix*Ly)+iy) == id[index_f] && IsData){
                    //std::cout<<"Indice i_f "<<index_f<<std::endl;
                    //std::cout<<"Id_ "<<id[index_f]<<std::endl;
                    //std::cout<<"t "<<t<<std::endl;
                    fnew[n0] = feq(rho_f[index_f], Ux0, Uy0, i);
                    if(i==Q-1){
                        index_f++;
                        index_tmp++;}
                    if((t+1)%(iter_per_hour*24)!=0)index_f=index_f-index_tmp;
                    if(index_f >= id.size()) IsData = false;

                }
                //else if(ix*Ly+iy == id[index_f-index_tmp]){

                //}
                else{
                    fnew[n0] = feq(rho0, Ux0, Uy0, i);
                }
            }
        }   
    }
}

void LatticeBoltzman::Advection(void){
    int ix, iy, i, ixnext, iynext, n0, n0next;

    for (ix = 0; ix < Lx; ix++){ // for each cell
        for (iy = 0; iy < Ly; iy++){
            for (i = 0; i < Q; i++){ // on each direction
                ixnext = ix + Vx[i];
                iynext = iy + Vy[i];
                //Free boundaries
                if (ixnext >= 0 && ixnext < Lx && iynext >= 0 && iynext < Ly) {
                    n0 = n(ix, iy, i);                
                    n0next = n(ixnext, iynext, i);
                    f[n0next] = fnew[n0]; 
                }
                else{
                    continue;
                }
            }
        }
    }
}

// Funciones de impresión de resultados y visualización
//---------------------Impresión de resultados---------------------
// Funcion para imprimir los resultados de la simulación en un archivo
void LatticeBoltzman::PrintData(std::string NameFile)
{
    // Abre un archivo de salida para guardar los datos de la simulación
    std::ofstream MyFile(NameFile);
    double rho0, Ux0, Uy0;
    int step = 1;

    // Recorre la cuadrícula en pasos de a step para ahorrar espacio en los archivos de salida
    for (int ix = 0; ix < Lx; ix += step)
    {
        for (int iy = 0; iy < Ly; iy += step)
        {
            // Calcula la densidad (rho0) y las velocidades (Ux0, Uy0) para cada celda
            rho0 = rho(ix, iy, false);      // Densidad en la celda (ix, iy)

            // Escribe los datos en el archivo de salida
            MyFile << ix << " " << iy << " "<< rho0 << std::endl;
        }
        MyFile << std::endl; // Inserta un salto de línea después de cada fila de la cuadrícula
    }

    // Cierra el archivo de salida
    MyFile.close();
}

// Funcion para imprimir un frame de la simulación
void LatticeBoltzman::PrintFrame(double t, int t_a)
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
    GnuplotScript << "set output 'frames/density_" << std::setw(3) << std::setfill('0') << t_a << ".png'" << std::endl; // Nombre del archivo de salida PNG, numerado según el tiempo de simulación
    GnuplotScript << "set pm3d map" << std::endl;                                                                     // Usar un mapa de calor (pm3d) para visualizar los datos
    GnuplotScript << "set size ratio -1" << std::endl;                                                                // Mantener una relación de aspecto cuadrada
    GnuplotScript << "set xrange [0:" << Lx << "]" << std::endl;                                                      // Definir los límites de los ejes x e y, basados en las dimensiones de la cuadrícula
    GnuplotScript << "set yrange [0:" << Ly << "]" << std::endl;
    GnuplotScript << "set cbrange [0:5000000]" << std::endl;                                                                        // Configurar la escala de colores (color bar) para los valores de densidad
    GnuplotScript << "set palette defined (0 'black', 1 'red', 2 'orange', 3 'yellow', 4 'white')" << std::endl;              // Definir la paleta de colores, de negro a blanco pasando por rojo, naranja y amarillo
    GnuplotScript << "set title 'Densidad en t = " << t << "'" << std::endl;                                                  // Definir el título del gráfico, basado en el tiempo de simulación actual
    GnuplotScript << "plot 'data/density_" << std::setw(3) << std::setfill('0') << t_a << ".dat' u 1:2:3 w image" << std::endl; // Instrucción para graficar los datos de densidad desde el archivo correspondiente
    GnuplotScript.close();                                                                                                    // Cerrar el archivo de script

    // Ejecutar el script de Gnuplot para generar la imagen PNG
    ret = system("gnuplot frame_script.gp");
    if (ret != 0)
        std::cerr << "Error: No se pudo ejecutar el script de gnuplot" << std::endl;
}

//----------------------Carga de resultados----------------------
// Función para cargar los datos necesarios para la simulación.
// La función intenta cargar los datos desde un archivo binario previamente guardado para acelerar la carga.
// Si el archivo binario no existe, carga los datos desde un archivo de texto y luego los guarda en formato binario para futuras ejecuciones.
void LoadData(std::string NameFile)
{
    // Nombre del archivo binario para almacenar los datos
    const std::string binaryFilename = "velocity.bin";

    // Calcula el tamaño total de los datos a cargar (número de elementos)
    int totalSize = Lx * Ly * t_hour;

    // Verifica si el archivo binario ya existe
    std::ifstream binFile(binaryFilename, std::ios::binary);
    if (binFile.good())
    {
        // Si el archivo binario existe, carga los datos desde él
        std::cout << "Cargando datos desde archivo binario..." << std::endl;
        LoadDataBinary(binaryFilename, Ux, Uy, totalSize); // Llama a la función para cargar los datos binarios
    }
    else
    {
        // Si el archivo binario no existe, carga los datos desde los archivos de texto originales
        std::cout << "Cargando datos desde archivos de texto..." << std::endl;
        std::ifstream Velocity(NameFile); // Abre el archivo de texto con los datos de velocidad
        std::ofstream debug("debug.txt"); // Abre un archivo de depuración para verificar la carga

        // Itera sobre la cuadrícula y el tiempo para cargar los datos en los arrays Ux y Uy
        for (int ix = 0; ix < Lx; ix++)
        {
            for (int iy = 0; iy < Ly; iy++)
            {
                for (int t = 0; t < t_hour; t++)
                {
                    int index = (ix * Ly + iy) * t_hour + t;                                 // Calcula el índice lineal para los arreglos
                    Velocity >> Ux[index] >> Uy[index];                                      // Lee los datos de velocidad desde el archivo de texto
                    debug << index + 1 << " " << Ux[index] << " " << Uy[index] << std::endl; // Escribe en el archivo de depuración
                }
            }
        }

        // Cierra los archivos abiertos
        Velocity.close();
        debug.close();

        // Guarda los datos en un archivo binario para acelerar futuras ejecuciones
        SaveDataBinary(binaryFilename, Ux, Uy, totalSize);
        std::cout << "Datos guardados en archivo binario para futuras ejecuciones." << std::endl;
    }

    std::cout << "Datos cargados exitosamente" << std::endl;
}

void LoadDataFires(std::string NameFile){
        // Carga los datos desde el archivo de texto
        std::ifstream coordFile(NameFile);

        if (!coordFile.is_open()) {
            std::cerr << "No se pudo abrir el archivo de coordenadas." << std::endl;
            return;
        }

        int temp_id, temp_rho_f;
        while (coordFile >> temp_id >> temp_rho_f) {
            id.push_back(temp_id);
            rho_f.push_back(temp_rho_f);
        }
        coordFile.close();
}

// Función para cargar los datos desde un archivo binario.
// Esta función se llama cuando el archivo binario ya existe, lo que permite cargar los datos más rápidamente.
void LoadDataBinary(const std::string &filename, double *Ux, double *Uy, int size)
{
    // Abre el archivo binario en modo lectura
    std::ifstream file(filename, std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Error: No se pudo abrir el archivo binario para cargar los datos." << std::endl;
        return;
    }

    // Lee los datos desde el archivo binario y los almacena en los arreglos Ux y Uy
    file.read(reinterpret_cast<char *>(Ux), size * sizeof(double));
    file.read(reinterpret_cast<char *>(Uy), size * sizeof(double));

    // Cierra el archivo binario
    file.close();
}

// Función para guardar los datos en un archivo binario.
// Esta función se llama después de cargar los datos desde el archivo de texto por primera vez, para guardar los datos en formato binario.
void SaveDataBinary(const std::string &filename, double *Ux, double *Uy, int size)
{
    // Abre el archivo binario en modo escritura
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open())
    {
        std::cerr << "Error: No se pudo abrir el archivo binario para guardar los datos." << std::endl;
        return;
    }

    // Escribe los datos de los arreglos Ux y Uy en el archivo binario
    file.write(reinterpret_cast<const char *>(Ux), size * sizeof(double));
    file.write(reinterpret_cast<const char *>(Uy), size * sizeof(double));

    // Cierra el archivo binario
    file.close();
}

// Funcion main
int main(int argc, char* argv[]){

    int t_a = 0;
    // Parámetros de la simulación
    int tframe = iter_per_hour/2, tmax = iter_per_hour*t_hour, delta_t = 1, ret; // tframe: intervalo entre frames, tmax: tiempo máximo de simulación

    // Densidad inicial y velocidad 
    double rho0 = 0.00001, Ux0 = 0.0, Uy0 = 0.0; 
    // Crear una instancia de la clase LatticeBoltzman
    LatticeBoltzman Air;

    LoadData("velocity.txt");
    LoadDataFires("coordenadasfuentesgrilla1410.txt");

    // Leer parámetros desde la línea de comandos: término de fuente y coeficiente de difusión
    D = std::stod(argv[1]); // Coeficiente de difusión (ejemplo: 0.016 -> tau = 0.548)

    // Calcular el tiempo de relajación tau basado en el coeficiente de difusión y el paso de tiempo
    tau = (D / delta_t * Cs2) + 0.5;

    // Calcular otros valores útiles basados en tau
    Utau = 1.0 / tau;
    UmUtau = 1 - Utau; // 1 - 1/tau

    // Inicializar la simulación con las condiciones iniciales (densidad y velocidades)
    Air.Start(rho0, Ux0, Uy0);

    // Bucle principal de la simulación
    for (int t = 0; t <= tmax; t++)
    {
        
        Air.Collision();
        Air.ImposeFields(t); // Ux0, Uy0);
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
            ss << "data/density_" << std::setw(3) << std::setfill('0') << t_a << ".dat"; // Nombre del archivo basado en el tiempo de simulación
            Air.PrintData(ss.str());                                                // Llamada a la función que guarda los datos

            // Generar un frame de la simulación usando Gnuplot
            Air.PrintFrame(t,t_a);
            t_a++;

            // Mostrar el porcentaje de avance de la simulación en la consola
            std::cout << "Porcentaje de avance: " << (t * 100) / tmax << "%" << std::endl;
        }
    }
    return 0;
}