#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip> // Include the <iomanip> header for std::setprecision

struct DataPoint {
    double K;
    double tMax;
    double TauMax;
};

int main() {
    
    std::vector<double> Ks = {1e8, 2e8, 5e8, 1e9, 2e9, 5e9, 1e10};
    // std::vector<double> Ks = {1e8, 1e9, 1e10};

    std::vector<DataPoint> dataPoints;

    double scale = 1e-10; // Escala para convertir a microsegundos

    // Abrir y leer los archivos
    for (const auto& K : Ks) {
       

        std::ostringstream oss;
        oss << std::scientific << std::setprecision(1) << K;
        std::string fileName = "Punto(tMax,TauMax)_K=" + oss.str() + ".txt";
        std::ifstream inFile(fileName);

        if (inFile.is_open()) {
            std::cout << "Leyendo archivo el archivo "<<K<< fileName << std::endl;
            // Leer los valores y guardarlos en dataPoints
            double tMax, TauMax;
            inFile >> tMax >> TauMax;// Leer los valores del archivo
            dataPoints.push_back({K*scale, tMax, TauMax*scale});
            inFile.close();

        } else {
            std::cout << "No se pudo abrir el archivo " << fileName << std::endl;
        }
    }
    std::cout << "Se leyeron " << dataPoints.size() << " puntos" << std::endl;
    // std::cout << "K tMax TauMax" << std::endl;
    // std::cout << dataPoints[0].K << " " << dataPoints[0].tMax << " " << dataPoints[0].TauMax << std::endl;
    // std::cout << dataPoints[1].K << " " << dataPoints[1].tMax << " " << dataPoints[1].TauMax << std::endl;
    // std::cout << dataPoints[2].K << " " << dataPoints[2].tMax << " " << dataPoints[2].TauMax << std::endl;
    // std::cout << dataPoints[3] << std::endl; 
    // std::vector<DataPoint> dataPoints = {
    //     {1e8, 0.17678, 1.15038e+07},
    //     {2e8, 0.17629, 1.51807e+07},
    //     {5e8, 0.1757, 2.19021e+07},
    //     {1e9, 0.17541, 2.88999e+07},
    //     {2e9, 0.17519, 3.81315e+07},
    //     {5e9, 0.17498, 5.50088e+07},
    //     {1e10, 0.17487, 7.26017e+07}
    //     // Agrega aquí más puntos si los tienes
    // };

    std::ofstream outFile("Tabla.txt");
    for (const auto& point : dataPoints) {
        outFile << point.K << " " << point.tMax << " " << point.TauMax << "\n";
    }

    return 0;
}
