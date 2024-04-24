# Configuración del terminal
set terminal png size 800,600
set output 'TauMax&tMaxVS_K.png'

# Títulos y etiquetas
set title 'Comparación de Datos en Escala Logarítmica'
set xlabel 'Columna 1'
set ylabel 'Columna 2 y 3'
set grid

# Escala logarítmica
set logscale xy

# Estilo de datos y línea de ajuste
set style data points
set style line 1 lc rgb 'red' lt 1 lw 2 pt 7 ps 1.5  # Estilo para los puntos de la columna 2
set style line 2 lc rgb 'blue' lt 1 lw 2 pt 5 ps 1.5 # Estilo para los puntos de la columna 3

# Funciones de ajuste
f1(x) = b * x + ct
f2(x) = a * x + ck

# Ajuste de regresión
fit f1(x) 'Tabla.txt' using (log10($1)):(log10($2)) via b,ct
fit f2(x) 'Tabla.txt' using (log10($1)):(log10($3)) via a,ck

# Plot
plot 'Tabla.txt' using 1:2 title sprintf('t') with points ls 1, \
     f1(x) title sprintf('Fit t: y = %.2fx + %.2f', b, ct) with lines ls 1, \
     'Tabla.txt' using 1:3 title sprintf('k') with points ls 2, \
     f2(x) title sprintf('Fit k: y = %.2fx + %.2f', a, ck) with lines ls 2
