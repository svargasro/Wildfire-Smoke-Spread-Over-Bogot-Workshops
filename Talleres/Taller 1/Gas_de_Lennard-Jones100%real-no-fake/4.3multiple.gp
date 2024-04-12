set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'RadiosDeGiro.png'
set xlabel 't' font 'Arial,14'
set ylabel 'x(t)' font 'Arial,14'
set style line 1 lt rgb '#1BA65A' lw 2 pt 7 ps 0.5

set grid linecolor rgb '#252525' linewidth 1.5
set key left top spacing 1.2 font 'Arial,12'
set border linewidth 1.5

set xtics font 'Arial,12'
set ytics font 'Arial,12'

plot 'GasKbT=10.txt' using 1:2 title 'GasKbT=10.0' with lines linecolor rgb '#5C3F81', \
     'LiquidoKbT=0.5.txt' using 1:2 title 'LiquidoKbT=0.5' with lines linecolor rgb '#338A4C', \
     'SolidoKbT=0.05.txt' using 1:2 title 'SolidoKbT=0.05' with lines linecolor rgb '#CF6E22' 

#Como vemos en general el radio de giro tiene en orden descendente los valores de los gases, liquidos y solidos, 
#lo cual concuerda con lo observado anteriormente. Cabe resaltar que los valores aleatorios dentro de la simulación
#pueden variar ligeramente los resultados, sin embargo, la mayoria de las veces se mantiene este orden.git add 