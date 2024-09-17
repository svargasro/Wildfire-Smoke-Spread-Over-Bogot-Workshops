set terminal pngcairo size 800,800 enhanced font 'Verdana,10'
set output 'frames/density_23000.png'
set pm3d map
set size ratio -1
set xrange [0:256]
set yrange [0:256]
set cbrange [0:*]
set palette defined (0 'black', 1 'red', 2 'orange', 3 'yellow', 4 'white')
set title 'Densidad en t = 23000'
plot 'data/density_23000.dat' u 1:2:3 w image
