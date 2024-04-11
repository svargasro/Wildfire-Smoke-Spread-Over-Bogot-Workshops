set terminal pdfcairo enhanced color font 'Verdana,10'
set output 'Modelo_SIR.pdf'

set xlabel 'Tiempo'
set ylabel 'Susceptibles, Infectados, Recuperados'
set title 'Modelo SIR'

#using 1:2 with lines title 'Susceptibles', \ #columna 1 y 2
#using 1:3 with lines title 'Infectados', \ #columna 1 y 3
#using 1:4 with lines title 'Recuperados' #columna 1 y 4
plot 'data.txt' using 1:2 with lines title 'Susceptibles', \
     'data.txt' using 1:3 with lines title 'Infectados', \
     'data.txt' using 1:4 with lines title 'Recuperados'

unset output
