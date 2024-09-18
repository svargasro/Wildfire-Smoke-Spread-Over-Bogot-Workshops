#Este archivo se encarga de crear una grilla en la imagen del mapa, 
#donde se encuentran las estaciones de monitoreo de la calidad del aire. 
#Ademas, se guardan las coordenadas de las estaciones en un archivo .txt

#Observaciones: Este archivo se ha creado para poder ejecutarlo desede
#una terminal de comandos, ya que no se puede ejecutar desde un notebook.

#Librerias
import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches


#path de la imagen
image_path = '../../Data/Map.png'

#Carga de la imagen
image = cv2.imread(image_path)

#Conversion de la imagen a RGB
image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)

#Conversion de la imagen a escala de grises
image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)

# Aplicar suavizado para mejorar la detección de círculos
image_blur = cv2.medianBlur(image_gray, 7)

# Definir el rango de valores de gris oscuro
lower_gray = 125
upper_gray = 255

# Crear una máscara para los píxeles oscuros
mask = cv2.inRange(image_blur, lower_gray, upper_gray)

# Aplicar la máscara a la imagen borrosa
masked_image = cv2.bitwise_and(image_blur, image_blur, mask=mask)


# Parámetros de la grilla
Lx = 300
Ly = int(Lx * 1.4)

# Detectar círculos usando la transformada de Hough
circles = cv2.HoughCircles(masked_image, cv2.HOUGH_GRADIENT, dp=1.5, minDist=20,
                           param1=50, param2=30, minRadius=13, maxRadius=21)

# Verificar si se han detectado círculos
if circles is not None:
    print(f'Se detectaron {circles.shape[1]} círculos.')
    
    # Redondear las coordenadas y el radio de los círculos
    circles = np.round(circles[0, :]).astype("int")
    
    # Crear una matriz para almacenar la grilla de posiciones de las estaciones
    height, width, _ = image.shape
    
    # Crear la matriz de posiciones de las estaciones
    stations_positions = np.zeros((Ly, Lx))

    # Crear una figura para visualizar la imagen y la grilla
    fig, ax = plt.subplots(figsize=(8, 12))
    ax.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))

    cell_height = height / Ly  # Número de filas
    cell_width = width / Lx  # Número de columnas
    
    # Configurar los ticks en los ejes según la grilla
    interval_ticks = Lx // 10  # Intervalo adaptativo en función de Lx
    x_ticks = [i * cell_width for i in range(0, Lx, interval_ticks)]
    y_ticks = [j * cell_height for j in range(0, Ly, interval_ticks)]

    # Establecer los ticks y las etiquetas en función de la grilla
    ax.set_xticks(x_ticks)
    ax.set_yticks(y_ticks)
    ax.set_xticklabels([str(i) for i in range(0, Lx, interval_ticks)])
    ax.set_yticklabels([str(j + 1) for j in reversed(range(0, Ly, interval_ticks))])

    # Procesar cada círculo detectado
    for (x, y, r) in circles:
        # Celda en la que se encuentra el círculo
        cell_x = int(x / cell_width)  # Convertir a entero
        cell_y = int((height - y) / cell_height)  # Invertir la coordenada y

        # Marcar la celda en la grilla, asegurándose de no exceder los límites
        if 0 <= cell_y < Ly and 0 <= cell_x < Lx:
            stations_positions[cell_y, cell_x] = 1
            # Anotar la posición del círculo con los índices corregidos

    # Colorear la grilla sobre la imagen
    for j in range(Lx):
        for i in range(Ly):
            # Dibujar un rectángulo azul en la celda si hay una estación
            if stations_positions[i, j] == 1:
                rect = patches.Rectangle((j * cell_width, (Ly -(i+1)) * cell_height), cell_width, cell_height, linewidth=1, edgecolor='b', facecolor='none', alpha=0.8)
                ax.add_patch(rect)
                # plt.text(j * cell_width + cell_width / 2, (Ly -(i+1)) * cell_height + cell_height / 2, f'({j}, {i})', color='black', fontsize=8, ha='center', va='center', alpha=1)
            else:
                # Dibujar un rectángulo rojo en la celda si no hay una estación
                rect = patches.Rectangle((j * cell_width, (Ly -(i+1)) * cell_height), cell_width, cell_height, linewidth=1, edgecolor='r', facecolor='none', alpha=0.1)
                ax.add_patch(rect)


#path del archivo donde se guardara las estaciones
path_data = '../../Data/estaciones.txt'

#Cargar las estaciones
data_estaciones = pd.read_csv(path_data)                

# Filtrar el dataset por la columna "VelViento" igual a True
data_VelViento = data_estaciones[data_estaciones['VelViento'] == True]

# Crear una matriz para almacenar la grilla de posiciones de las estaciones
coords_X, coords_Y = [], []
coords_X, coords_Y = np.where(stations_positions == 1)                

# Crear una matriz para almacenar las coordenadas de las estaciones
data_estaciones_coords = data_VelViento[['Estación', 'x_cell', 'y_cell']].copy()
data_estaciones_coords.reset_index(drop=True, inplace=True)

# Asignar las coordenadas de las estaciones a la matriz
data_VelViento.to_csv('../../Data/estaciones_coord.txt', sep=',', index=False)