set terminal pngcairo size 800,600 
set output 'grafico4.3.png'
set xlabel 't' font 'Verdana,14'
set ylabel 'R_{giro}' font 'Verdana,14'
set style line 1 lt rgb '#7F4A9B' lw 2 pt 7 ps 0.35
set key left top spacing 1.2 font 'Verdana,12'

# Obtener el último valor de la columna 1 (asumiendo que el título está en la última línea)
last_value = system("tail -n 1 4.3.txt | awk '{print $1}'")
penultimate_value = system("tail -n 2 4.3.txt | head -n 1 | awk '{print $1}'")

set title "Gráfica de R_{giro} en función del tiempo, para KbT = "  . penultimate_value font 'Verdana,18'
set label "R_{giro} promedio = " . last_value at screen 0.93,0.83 right textcolor rgb '#222756'

set arrow 1 from graph 0, first last_value to graph 1, first last_value nohead lc rgb '#222756' lw 2
plot '4.3.txt' using 1:2 with linespoints linestyle 1 title 'Datos'
