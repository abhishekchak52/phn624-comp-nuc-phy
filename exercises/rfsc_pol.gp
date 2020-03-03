unset key
set xrange [-100:100]
set yrange [-200:200]

set polar
set size square

set object 1 circle at 0,0 size 7 fc rgb "black" lw 1

pl for [i=1:*] 'rfsc_pol.txt' u i:(column(i+1)) w l lw 3
