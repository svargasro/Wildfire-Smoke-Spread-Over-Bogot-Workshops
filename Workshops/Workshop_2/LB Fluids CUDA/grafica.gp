# Establece el terminal de salida para una imagen PNG
set terminal pngcairo enhanced font 'Verdana,10'
set output 'Ondas.png'

# Define los rangos para las variables x e y
set xrange [-0:128]
set yrange [-0:128]
set zrange [-20:20]

# Configuración de la vista y la malla de la superficie
set view 60, 30, 1, 1
set isosample 50, 50


# Configuración de etiquetas y título
set xlabel 'X axis'
set ylabel 'Y axis'
set zlabel 'Z axis'


# Grafica la función
splot "Ondas.dat" w l
