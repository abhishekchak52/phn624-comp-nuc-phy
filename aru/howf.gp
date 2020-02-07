set xlabel '$x$'
set ylabel '$\psi_n(x)$'
set xrange [-5:5]
set yrange [-0.8:0.8]
set mxtics 
set mytics
set xtics -5,1,5
set ytics -0.8,0.2,0.8
  pl 'howf.out' u 1:2 w l lw 3 t '$n=0$
repl 'howf.out' u 1:3 w l lw 3 t '$n=1$
repl 'howf.out' u 1:4 w l lw 3 t '$n=2$
repl 'howf.out' u 1:5 w l lw 3 t '$n=3$
repl 'howf.out' u 1:6 w l lw 3 t '$n=4$
repl 'howf.out' u 1:22 w l lw 3 t '$n=20$
pause -1
set term epslatex size 11.3cm,7cm color colortext standalone
set output "gp.tex"
repl
