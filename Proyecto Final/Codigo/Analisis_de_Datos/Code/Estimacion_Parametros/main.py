import ajuste
import obs_data
import numpy as np


if __name__=="__main__":
    # Ejemplo de uso
    Lx = 100
    Ly = (1.4)*Lx
    sources_1D = []  # Lectura ubicación de todas las fuentes para todas las horas.
    initial_rho_sources = [] #Valores iniciales de rho de las fuentes para todas las horas.
    velocityData = [] #Valores de velocidades para todas las horas.

    stations = obs_data.coord_est('../../Data/estaciones_coord.txt') #Ubicación de las estaciones.
    PM_estaciones = obs_data.PM_base('../../Data/PM_base.csv') #Datos de contaminación en estaciones.

    iter_per_hour = 9604
    timesAdjusment = np.arange(0,9604*24,9604/2) #Lectura de todos los tiempos para ajuste.

    #timesAdjusment[0,9604/4,(9604/4)*2]

    for i, t in enumerate(timesAdjusment):
        # Optimización de los parámetros de las fuentes
        tInicial = timesAdjusment[i]
        if(i != len(timesAdjusment)):
            tFinal = timesAdjusment[i+1]

        if(i==0): #El Lattice Boltzmann debe usar start solo cuando t=0, en otro caso, usa Load para cargar fnew.
            LatticeBoltzmann(tInicial, 1, ID_sources, rho_sources, velocidades) #Genera un archivo de fnew para todas las celdas linealizadas.

        rho_opt = optimize_parameters_for_cells(stations, PM_estaciones, sources_1D,initial_rho_sources)

        density_opt= LB_density_with_sources(Lx, Ly, sources_1D, rho_opt, sigma)

        LatticeBoltzmann(int tinicial,int tfinal, ID_Sources,  rho_sources, velocidades)

        print("Valores óptimos de rho:", rho_opt, rho_sources_random)




    rho_opt = ajuste.optimize_parameters_for_cells(stations, rho_obs, sources_1D, initial_rho_sources, Lx, Ly, sigma)
    density_opt = ajuste.LB_density_with_sources(Lx, Ly, sources_1D, rho_opt, sigma)

    print("Valores óptimos de rho:", rho_opt, rho_sources_random)
