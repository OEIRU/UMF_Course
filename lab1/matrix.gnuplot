set title 'Sparsity pattern'
set xlabel 'Column'
set ylabel 'Row'
plot 'matrix.dat' using 2:1 with points pt 7 ps 0.5 title ''
pause -1
