import pandas as pd

def coord_est():
    estaciones = pd.read_csv('../../Data/estaciones_coord.txt')
    return estaciones['linealization']

