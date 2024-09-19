import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def index_to_coordinates(index, Ly):
    """Convierte un índice unidimensional a coordenadas bidimensionales dado el número de filas (Ly)."""
    x = index // Ly
    y = index % Ly
    return (x, y)

def coordinates_to_index(x, y, Ly):
    """Convierte coordenadas bidimensionales a un índice unidimensional dado el número de columnas (Lx)."""
    return x * Ly + y

def gaussian_2d(x, y, x0, y0, rho, sigma):
    """Calcula la gaussiana en el punto (x, y) dada una fuente en (x0, y0) con valor rho y ancho sigma."""
    return rho * np.exp(-((x - x0)**2 + (y - y0)**2) / (2 * sigma**2))

def LB_density_with_sources(Lx, Ly, sources_1D, rho_sources, sigma):
    """Simula la densidad en una cuadrícula con fuentes gaussianas y devuelve el arreglo bidimensional."""
    # Crear una cuadrícula de densidad inicializada en ceros
    density_grid = np.zeros((Lx, Ly)) #CAMBIO.  # Nota el orden (Ly, Lx) para coincidir con (Y, X)
    num_sources = len(sources_1D)
    sources = [(0, 0, 0)] * num_sources #Tripleta que guardará (x,y,rho)

    # Asignar valores a cada tripleta
    for i in range(num_sources):
        x, y = index_to_coordinates(sources_1D[i], Ly) #CAMBIO.
        sources[i] = (x, y, rho_sources[i])  #CAMBIO. # Nota: (y, x) en lugar de (x, y) porque se usa (x, y) en gaussian_2d
    # Recorrer todas las celdas de la cuadrícula
    for ix in range(Lx):
        for iy in range(Ly):
            #Para cada fuente, sumar la contribución gaussiana a la cuadrícula
            for (x0, y0, rho) in sources:
                density_grid[ix, iy] += gaussian_2d(ix, iy, x0, y0, rho, sigma)

    #Código de prueba funcionamiento.
    # for ix in range(Lx):
    #     for iy in range(Ly):
    #         density_grid[ix,iy]=0
    # for (x0,y0,rho) in sources:
    #     density_grid[x0,y0] = rho
    return density_grid

def linearize_grid(density_grid):
    """Convierte un arreglo bidimensional en un arreglo unidimensional recorriendo por filas."""
    # Obtener las dimensiones de la cuadrícula
    Lx, Ly = density_grid.shape #CAMBIO

    # Crear un arreglo unidimensional para almacenar los datos linealizados
    linearized_array = np.zeros(Lx * Ly)

    # Rellenar el arreglo unidimensional recorriendo la cuadrícula fila por fila. #CAMBIO
    for x in range(Lx):
        for y in range(Ly):
            linearized_array[x * Ly + y] = density_grid[x, y]

    return linearized_array

def unlinearize_grid(linearized_array, Lx, Ly):
    """Convierte un arreglo unidimensional en un arreglo bidimensional dado el número de columnas (Lx) y filas (Ly)."""
    # Crear un arreglo bidimensional para almacenar los datos
    density_grid = np.zeros((Lx, Ly))

    # Rellenar el arreglo bidimensional recorriendo el arreglo unidimensional fila por fila
    for x in range(Lx):
        for y in range(Ly):
            density_grid[x, y] = linearized_array[x * Ly + y]

    return density_grid

def objective(rho_sources_opt, stations, density_observed, Lx, Ly, sources_1D, sigma):
    """
    Calcula el error entre la densidad observada y la densidad simulada solo en las celdas que tienen estaciones.

    Parameters:
    - rho_sources_opt: valores optimizados de rho para las celdas fuente.
    - stations: lista de índices de celdas con estaciones.
    - density_observed: arreglo unidimensional de densidad observada.
    - Lx, Ly: dimensiones de la grilla.
    - sources_1D: lista de índices de celdas con fuentes.
    - sigma: parámetro que afecta la distribución de densidad.
    
    Returns:
    - error: suma de los cuadrados de las diferencias en las celdas de las estaciones.
    """
    # Recalcular la densidad simulada usando los valores optimizados de rho
    density_grid_simulated_opt = LB_density_with_sources(Lx, Ly, sources_1D, rho_sources_opt, sigma)
    density_simulated_opt = linearize_grid(density_grid_simulated_opt)

    # Extraer las densidades en las estaciones de la densidad observada y simulada optimizada
    observed_stations = density_observed[stations]
    simulated_stations_opt = density_simulated_opt[stations]

    # Calcular el error cuadrático medio entre las densidades observadas y simuladas en las estaciones
    error = np.sum((observed_stations - simulated_stations_opt) ** 2)
    return error

def optimize_parameters_for_cells(stations, density_observed, sources_1D, rho_sources, Lx, Ly, sigma):
    """
    Optimiza los valores de rho para minimizar el error en las estaciones.
    
    Parameters:
    - stations: lista de índices de celdas con estaciones.
    - density_observed: densidad observada en forma de arreglo unidimensional.
    - sources_1D: índices de las celdas fuente.
    - rho_initial: valores iniciales de rho para las celdas fuente.
    - Lx, Ly: dimensiones de la grilla.
    - sigma: parámetro que afecta la simulación de la densidad.
    
    Returns:
    - rho_opt: valores optimizados de rho.
    """
    rho_initial = rho_sources  # Valores iniciales para las celdas fuente
    bounds = [(0.5, 4.0)] * len(sources_1D)  # Restricciones sobre los valores de las celdas fuente

    # Ejecutar la optimización
    result = minimize(objective, rho_initial, args=(stations, density_observed, Lx, Ly, sources_1D, sigma), 
                      method='L-BFGS-B', bounds=bounds)

    # Resultados óptimos
    rho_opt = result.x
    return rho_opt

def plot_density_grid(density_grid,Lx,Ly):
    """Grafica el arreglo bidimensional de densidad."""
    plt.imshow(density_grid, cmap='hot', origin='lower')
    plt.colorbar(label='Densidad')
    plt.title('Distribución de Densidad en la Cuadrícula')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('TemporalImage.png')
    #plt.show()
