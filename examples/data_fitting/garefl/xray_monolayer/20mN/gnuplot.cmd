
set timestamp
set bars
set title 'Current best fit'
set xlabel 'Q (inverse Angstroms)'
set ylabel 'RQ^4'
set logscale y




plot 'fit0.dat' u 1:($3*$1**4):2:($4*$1**4) t 'Model 0' w xyerrorbars, 'fit0.dat' u 1:($5*$1**4) not w lines

