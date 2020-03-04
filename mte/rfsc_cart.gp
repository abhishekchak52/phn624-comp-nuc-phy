unset key
set xrange [-100:100]
set yrange [-200:200]

set object 1 circle at 0,0 size 7 fc rgb "black" lw 1
pl for [i=1:40] 'rfsc_cart.txt' u (column(2*i)):(column(2*i+1)) w l lw 3
