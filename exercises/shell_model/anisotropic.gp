set xlabel '$\delta$'
set ylabel '$E(\MeV\)$'
set xrange [-0.3:0.3]
!set yrange [10:100]
set mxtics 
set mytics
set xtics -0.3,0.1,0.3
set nokey
pl for [col=2:20] 'anisotropic_ho2.out' u 1:col w l lw 1
pause -1
!set term epslatex size 10cm,16.2cm color colortext standalone
!set output "gp.tex"
repl

