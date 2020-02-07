set xlabel '$\delta$'
set ylabel '$E$ (MeV)'
set xrange [-0.15:0.15]
set yrange [-0.5:1.0]
set mxtics 
set mytics
set xtics -0.3,0.1,0.3
set key autotitle columnheader
set key right bottom
set title 'Zr isotopes'
pl for [col=2:7] 'pes.out' u 1:col w l lw 3 
pause -1
set term epslatex size 10.5cm,16cm color colortext standalone
set output "gp.tex"
repl
