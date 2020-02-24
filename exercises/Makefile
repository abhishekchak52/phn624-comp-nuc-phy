swe:
	ifort -o test_swe numerical_integration.f90 swe_inf_pot_well.f90

swe_plot:	swe
	./test_swe 
	gnuplot -e "set key autotitle columnhead;pl for [col=2:6] 'wf_data.txt' u 1:col w l lw 3; pause -1"

rfsc:
	ifort -o test_rfsc numerical_integration.f90 rfscat_cart.f90

rfsc_plot: rfsc
	./rfsc_cart.sh
	gnuplot rfsc_cart.gp
	
	
