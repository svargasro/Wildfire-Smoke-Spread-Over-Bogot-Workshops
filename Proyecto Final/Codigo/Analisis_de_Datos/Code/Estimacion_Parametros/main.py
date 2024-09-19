import ajuste 
import obs_data
import numpy as np

if __name__=="__main__":
    # Ejemplo de uso
    Lx = 15
    Ly = 15
    sigma = 5.0
    sources_1D = [33, 141, 175,23]  # Índices de celdas con incendio.
    initial_rho_sources = [1.0, 2.0, 1.5, 3.4] #Valores iniciales de rho (Investigación Alejandra)
    stations = obs_data.coord_est()   #Ubicación de las estaciones.
    #stations = [21, 70, 81, 122, 133, 3] #Ubicación de las estaciones.
    rho_sources_random = np.random.uniform(0.5, 4.0, size=4) #Valores de rho reales (Random para probar varias veces)

    # Simular la densidad y obtener una grilla para plot (Densidades en el espacio)
    density_grid_simulated = ajuste.LB_density_with_sources(Lx, Ly, sources_1D, initial_rho_sources, sigma)
    density_grid_obs = ajuste.LB_density_with_sources(Lx, Ly, sources_1D, rho_sources_random, sigma)

    #Plots de las densidades.
    #plot_density_grid(density_grid_simulated,Lx,Ly)
    #plot_density_grid(density_grid_obs,Lx,Ly)

    rho_simulated = ajuste.linearize_grid(density_grid_simulated)
    rho_obs = ajuste.linearize_grid(density_grid_obs)

    # Optimización de los parámetros de las fuentes, pasando Lx, Ly y sigma
    rho_opt = ajuste.optimize_parameters_for_cells(stations, rho_obs, sources_1D, initial_rho_sources, Lx, Ly, sigma)
    density_opt = ajuste.LB_density_with_sources(Lx, Ly, sources_1D, rho_opt, sigma)

    print("Valores óptimos de rho:", rho_opt, rho_sources_random)