# for run with gnuplot
set view map
set contour
set cntrparam levels 50
unset key
unset surface
set xtics 0,60,360
set xlabel "critical angle"
set ylabel "a_{pla}(au)"
splot [0:360][] "hamiltplares.dat" u 2:1:3 w l
pause -1

