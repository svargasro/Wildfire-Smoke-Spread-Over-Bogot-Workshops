#include <iostream>
#include <vector>

using namespace std;

// Definición de la función SIR
vector<double> SIR(double t, vector<double>& statevector, vector<double>& params) {
    double gamma = params[0];
    double beta = params[1];

    double s = statevector[0];
    double i = statevector[1];
    double r = statevector[2];

    double ds = -beta * s * i;
    double di = beta * s * i - gamma * i;
    double dr = gamma * i;

    return {ds, di, dr};
}

// Implementación de Runge-Kutta de cuarto orden
vector<vector<double>> RK4(vector<double>(*func)(double, vector<double>&, vector<double>&),
                                              vector<double>& initial_state, double t0, 
											  double tf, double h, vector<double>& params) {
    
	int num_steps = static_cast<int>((tf - t0) / h) + 1;
    vector<vector<double>> states(num_steps, vector<double>(initial_state.size()));
    vector<double> t(num_steps);

    states[0] = initial_state;
    t[0] = t0;

    for (int i = 1; i < num_steps; ++i) {
        vector<double> k1 = h*func(t[i - 1], states[i - 1], params);
        vector<double> k2 = h*func(t[i - 1] + 0.5 * h, states[i - 1], params);
        vector<double> k3 = h*func(t[i - 1] + 0.5 * h, states[i - 1], params);
        vector<double> k4 = h*func(t[i - 1] + h, states[i - 1], params);

        for (size_t j = 0; j < initial_state.size(); ++j) {
            states[i][j] = states[i - 1][j] + (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) / 6.0;
        }

        t[i] = t[i - 1] + h;
    }

    return states;
}

int main() {
    // Parámetros y condiciones iniciales
    vector<double> params = {0.1, 0.2}; // gamma, beta
    vector<double> initial_state = {0.9, 0.1, 0.0}; // s, i, r
    double t0 = 0.0;
    double tf = 500.0;
    double h = 5.0e1;

    // Aplicar RK4 para resolver el sistema de ecuaciones
    vector<vector<double>> states = RK4(SIR, initial_state, params, t0, tf, h);

    // Imprimir los resultados
    for (size_t i = 0; i < states.size(); ++i) {
        cout << "t = " << states[i][0] << ", s = " << states[i][0] << ", i = " << states[i][1] << ", r = " << states[i][2] << endl;
    }

    return 0;
}

//Plotearlo en gnuplot