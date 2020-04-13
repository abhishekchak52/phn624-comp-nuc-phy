set term png;
set output 'fes.png'
set xlabel "Deformation"
set ylabel "Free Energy"
set yrange [-0.5:1];
set xrange [-0.13:0.13];
set key autotitle columnheader;
pl for [col=2:7] 'pes_data.txt' u 1:col w l
set term x11
repl
