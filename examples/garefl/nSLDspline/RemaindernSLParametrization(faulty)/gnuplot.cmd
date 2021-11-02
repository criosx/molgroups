
set timestamp
set title 'Current best profile for each model'
set xlabel 'Depth (Ang)'
set ylabel 'rho (number density)'



plot 'profile0.dat' u 1:2 t 'rho 0' w l, 'profile1.dat' u 1:2 t 'rho 1' w l, 'profile2.dat' u 1:2 t 'rho 2' w l, 'profile3.dat' u 1:2 t 'rho 3' w l

