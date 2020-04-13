set title "$\Delta vs G$"
set nokey
set xlabel "G multiplier"
set ylabel "$\Delta$"
set yrange [0.0:4]
set xrange [0:1.5]
set term  png
set output "delta_vs_G.png"
pl 'delta_vs_g.txt' u 1:3 w l
clear;
set key
set output 'Epc_def.png'
set xrange [-0.3: 0.3]
set yrange [0:4]
set xlabel "$\delta$"
set ylabel "$E_{pc}$"
set title "Epc vs deformation"
pl 'bcs_delta.txt' u 1:3 w l title "Neutrons", '' u 1:4 w l title "Protons"
clear;
set output "pes_correction.png"
set xrange [-0.29:0.29]
set yrange [-1:5]
set xlabel "$\delta$"
set ylabel "$E (MeV)$"
set title "Correction to PES due to pairing correlation"
pl 'bcs_delta.txt' u 1:2 w l title "PES", '' u 1:($2+$3+$4) w l title "PES w/ correlation"
#repl
