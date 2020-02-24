swe:
	ifort -o test_swe numerical_integration.f90 swe_inf_pot_well.f90

swe_plot:	swe
	./test_swe 
	gnuplot -e "set key autotitle columnhead;pl for [col=2:6] 'wf_data.txt' u 1:col w l lw 3; pause -1"

rfsc:
	ifort -o test_rfsc numerical_integration.f90 rfscat_cart.f90

rfsc_plot: rfsc
	./rfsc_cart.sh
	gnuplot -e "unset key;set xrange [-100:100];set yrange [-200:200];pl for [i=3:*:2] 'rfsc_cart.txt' u (column(i-1)):i w l; pause -1"
	
	
