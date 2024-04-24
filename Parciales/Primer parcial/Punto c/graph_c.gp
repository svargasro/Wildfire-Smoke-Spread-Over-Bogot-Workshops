set terminal png size 800,600 enhanced font 'Verdana,10'
set output 'CunaNewtonTau_AllK.png'
set title 'Time vs Tau para Varios K'
set xlabel 'Time'
set ylabel 'Tau '
set grid
set style data points # It uses points to plot the data


plot \
'CunaNewtonTau_K=1.0e+08.txt' using 1:2 title sprintf('K = %.2e', 1.0e+08) with linespoints pointtype 7 pointsize 0.8, \
'CunaNewtonTau_K=2.0e+08.txt' using 1:2 title sprintf('K = %.2e', 2.0e+08) with linespoints pointtype 7 pointsize 0.8, \
'CunaNewtonTau_K=5.0e+08.txt' using 1:2 title sprintf('K = %.2e', 5.0e+08) with linespoints pointtype 7 pointsize 0.8, \
'CunaNewtonTau_K=1.0e+09.txt' using 1:2 title sprintf('K = %.2e', 1.0e+09) with linespoints pointtype 7 pointsize 0.8, \
'CunaNewtonTau_K=2.0e+09.txt' using 1:2 title sprintf('K = %.2e', 2.0e+09) with linespoints pointtype 7 pointsize 0.8, \
'CunaNewtonTau_K=5.0e+09.txt' using 1:2 title sprintf('K = %.2e', 5.0e+09) with linespoints pointtype 7 pointsize 0.8, \
'CunaNewtonTau_K=1.0e+10.txt' using 1:2 title sprintf('K = %.2e', 1.0e+10) with linespoints pointtype 7 pointsize 0.8, 


