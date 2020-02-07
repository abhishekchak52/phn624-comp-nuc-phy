set xlabel '$\delta$'
set ylabel '$E(\hbar\omega_0)$'
set xrange [-0.3:0.3]
set yrange [2:5]
set mxtics 
set mytics
set xtics -0.3,0.1,0.3
set nokey
pl for [col=2:35] 'nilsson.out' u 1:col w l lw 3 
pause -1
set term epslatex size 10cm,16.2cm color colortext standalone
set output "gp.tex"
repl
