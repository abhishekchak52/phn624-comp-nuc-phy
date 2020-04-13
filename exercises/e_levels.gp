set xrange [-1:10]
set yrange [-40:-24]
set ylabel "Energy [MeV]" font "arial, 24" rotate by 90
unset xtics
unset key
plot 'e_levels.txt' using 1:2 with lines lw 4
