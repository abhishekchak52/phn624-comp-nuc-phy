set nokey
set title 'Estimate of Largest Positive Charge Radius from alpha scattering by Au-197'
set xlabel '$E_{\alpha}(MeV)$' 
set ylabel '$R_{max}(fm)$'
pl 'data.txt' u 1:2 w l lw 2
pause -1
set term epslatex size 15cm,10cm color colortext standalone
set output "rfsc_mte.tex"
repl
