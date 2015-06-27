dim=`cat dim-superficie.dat`  #dimension de la superficie
set parametric
set xrange [0:*]
set yrange [0:*]
plot "Mfrente.dat" u 1:($3==0.5 ? $2:1/0) w l t "",\\
     _i,t w l t ""

#
#set multiplot layout 2,1
#unset key
#_i=30
#set pm3d corners2color c1 map
#set palette defined (-1 "white", 0 "black", 0.5 "red" ,1 "green")
#set size ratio -1
#splot "estado.dat" i _i u 1:2:3 w pm3d
#unset pm3d
#set arrow from first _i, graph 0 to first _i, graph 1 nohead
#set size ratio 0
#set yrange [0:dim+0.5]







