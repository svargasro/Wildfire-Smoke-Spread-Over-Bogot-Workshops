Ly = 256
Lx = 256
rho = 100.0

# Establece el terminal de salida para una imagen PNG
set terminal pngcairo enhanced font 'Verdana,10'
set output 'AdvectionDifusion.png'

# Define los rangos para las variables x e y
set xrange [-0:Lx]
set yrange [-0:Ly]
set zrange [0:rho]

# Configuración de la vista y la malla de la superficie
set view 60, 30, 1, 1
set isosample 50, 50


# Configuración de etiquetas y título
set title 'Advección y Difusión'
set xlabel 'X axis'
set ylabel 'Y axis'
set zlabel 'Z axis'


# Grafica la función
splot "AdvectionDifusion.dat" using 2:3:4 wit lines palette
