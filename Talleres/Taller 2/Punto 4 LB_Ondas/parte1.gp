set pm3d map
set size ratio 1
set terminal jpeg enhanced font "Arial,12"
set output "parte1.jpg"

# Paleta 'viridis' emulada
set palette defined (0 "#440154", 1 "#482777", 2 "#3E4A89", 3 "#31688E", 4 "#26838F", 5 "#1F9D8A", 6 "#6CCE59", 7 "#B6DE2B", 8 "#FDE725")
set title "Refracción de ondas planas" font "Arial,16"
set xlabel "Eje X" font "Arial,14"
set ylabel "Eje Y" font "Arial,14"
splot "parte1.dat"
