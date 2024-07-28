Ly = 256
Lx = 256
tmax = 100
rho = 100.0

# Establece el terminal de salida para una imagen PNG
set terminal pngcairo size 800,800 enhanced font 'Verdana,10'

# Define los rangos para las variables x e y
set xrange [0:Lx]
set yrange [0:Ly]
set zrange [0:rho]  # Ajusta esto según el rango esperado de rho0

# Configuración de la vista y la malla de la superficie
set view 60, 30, 1, 1

# Configuración de etiquetas y título
set xlabel 'X axis'
set ylabel 'Y axis'
set zlabel 'Density'

# Define el archivo de datos
file = 'Density.dat'

# Bucle para generar cada frame
do for [t=0:tmax] {
    set title sprintf("Time step: %d", t)
    set output sprintf('frames/frame_%04d.png', t)
    splot file using (column(1)==t ? column(2) : 1/0):(column(1)==t ? column(3) : 1/0):(column(1)==t ? column(4) : 1/0) with points pt 7 ps 1 palette
}
