set terminal pngcairo size 800,600 enhanced font 'Arial,12'
set output 'grafico4.1.png'
set xlabel 't' font 'Arial,14'
set ylabel 'x(t)' font 'Arial,14'
set style line 1 lt rgb '#1BA65A' lw 2 pt 7 ps 0.5

set grid linecolor rgb '#252525' linewidth 1.5
set key left top spacing 1.2 font 'Arial,12'
set border linewidth 1.5

set xtics font 'Arial,12'
set ytics font 'Arial,12'

set title "Partícula bidimensional bajo potencial de Lennard-Jones" font 'Arial,18'
plot '4.1.txt' using 1:2 with linespoints linestyle 1 title 'Datos' 