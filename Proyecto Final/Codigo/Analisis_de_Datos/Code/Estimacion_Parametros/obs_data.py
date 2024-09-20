import numpy as np
import pandas as pd

def coord_est(nombre_archivo):
    estaciones = pd.read_csv(nombre_archivo)
    return estaciones['linealization']


def PM_base(nombre_archivo):
    # Leer el archivo CSV y extraer los valores numéricos
    df = pd.read_csv(nombre_archivo, header=None)
    
    # Extraer los valores numéricos en formato científico
    datos = df.iloc[:, 0].astype(float)
    
    # Crear el array unidimensional
    array_unidimensional = np.array(datos.tolist())
    
    return array_unidimensional

