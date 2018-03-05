#!/usr/bin/gnuplot -persist
set size ratio -1
set key off
set title ARG1
unset border
#unset tics
splot './dat/'.ARG1.'_0.dat' using 1:2:3 w l lw 3 lt rgb 'black', './dat/'.ARG1.'_1.dat' using 1:2:3 w l lw 3 lt rgb 'red'
set view equal xyz
pause -1
