#!/usr/bin/gnuplot

#set term gif animate transparent opt delay 100 size 600,600 x000000
#set out "avance.gif"
#reset

cd "../HAMF_0809"
pwd

print a
i=0

#anima
imax="`grep '#nt' estado.dat|wc -l`"
imax=100
print imax
di=1
load 'anima'

i=i+1; replot

#grafico del avance en un plano
reset
set pm3d corners2color c1 map
set size ratio -1
i=0
splot "estado.dat" i i u 1:2:3 w pm3d


#grafico del avance en una superficie
reset
set pm3d corners2color c1
set size ratio -1
i=0
splot "estado.dat" i i u 1:2:4:3 w pm3d


#grafico de la dimension del frente
reset
unset multiplot
sep=30                          #separacion entre graficos
dim="`cat dim-superficie.dat`"  #dimension de la superficie
set yrange[*:*]
set key b r
plot "frente.dat" u 1:($2+sep*$3) w l t "",\
    "" u 1:(sep*$3) w l   t "d=0",\
    "" u 1:(sep*$3+1) w l t "d=1",\
    "" u 1:(sep*$3+2) w l t "d=2",\
    "" u 1:(sep*$3+dim) w l t 'd=dim'


#grafico de la dimension de los sitios quemados
reset
unset multiplot
sep=30
dim="`cat dim-superficie.dat`"  #dimension de la superficie
set yrange[*:*]
set key b r
plot "quemados.dat" u 1:($2+sep*$3) w l t "",\
    "" u 1:(sep*$3) w l   t "d=0",\
    "" u 1:(sep*$3+1) w l t "d=1",\
    "" u 1:(sep*$3+2) w l t "d=2",\
    "" u 1:(sep*$3+dim) w l t 'd=dim'


#grafico de la dimension de los sitios activos
reset
unset multiplot
sep=30
dim="`cat dim-superficie.dat`"  #dimension de la superficie
set yrange[*:*]
set key b r
plot "activos.dat" u 1:($2+sep*$3) w l t "",\
    "" u 1:(sep*$3) w l   t "d=0",\
    "" u 1:(sep*$3+1) w l t "d=1",\
    "" u 1:(sep*$3+2) w l t "d=2",\
    "" u 1:(sep*$3+dim) w l t 'd=dim'

#grafico de la "densidad" de percolacion
reset
set yrange[*:*]
plot "percolacion.dat" u 1:2 w lp t ""

#grafico de la "derivada de la densidad" de percolacion
reset
set yrange[*:*]
!echo "0 0 0" > .gnuplot.tmp
!cat percolacion.dat >> .gnuplot.tmp
plot "<paste percolacion.dat .gnuplot.tmp" u 1:($5-$2)/($4-$1) w lp t ""

#
#grafico de la dimension y el avance del fuego
#

reset
set multiplot layout 2,1
unset key
dim="`cat dim-superficie.dat`"  #dimension de la superficie
i=30
set pm3d corners2color c1 map
set palette defined (-1 "white", 0 "black", 0.5 "red" ,1 "green")
set size ratio -1
splot "estado.dat" i i u 1:2:3 w pm3d
unset pm3d
set arrow from first i, graph 0 to first i, graph 1 nohead
set size ratio 0
set yrange [0:dim+0.5]
plot "frente.dat" u 1:($3==0.6 ? $2:1/0) w l t "",\
     1 w l t "d=1",\
     2 w l t "d=2",\
     2.443 w l t "d=dim"

i=i+1; replot










