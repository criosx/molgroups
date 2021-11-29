
set timestamp
set title 'Current best profile for each model'
set xlabel 'Depth (Ang)'
set ylabel 'rho (number density)'



plot 'profile0.dat' u 1:2 t 'rho 0' w l

