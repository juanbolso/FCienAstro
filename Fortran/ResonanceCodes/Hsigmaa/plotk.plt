set view map
set contour
set cntrparam levels 200
unset key
unset surface
set xtics 0,60,360
set xlabel "sigma"
set ylabel "a (au)"
splot [0:360][] "hamilto.dat" u 2:1:3 w l
pause -1

