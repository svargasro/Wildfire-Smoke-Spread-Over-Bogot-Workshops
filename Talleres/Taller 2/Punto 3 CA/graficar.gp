set terminal png size 800,600
set output 'CA_Difusion1D.png'
set title "Varianza vs time"
set xlabel "time"
set ylabel "Varianza"

plot "punto_a.txt" w l
