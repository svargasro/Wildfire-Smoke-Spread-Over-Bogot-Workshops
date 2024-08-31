set terminal pngcairo size 800,800 enhanced font 'Verdana,10'
set output 'frames/density_990.png'
set pm3d map
set size ratio -1
set xrange [0:100]
set yrange [0:100]
set cbrange [0:*]
set palette defined (0 'black', 1 'red', 2 'yellow', 3 'white', 4 'red')
set title 'Densidad en t = 990'
plot 'data/density_990.dat' u 1:2:3 w image
