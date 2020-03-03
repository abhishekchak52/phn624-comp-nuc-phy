fort_comp = ifort -no-wrap-margin
swe:
	${fort_comp} -o test_swe numerical_integration.f90 swe_inf_pot_well.f90

swe_plot:	swe
	./test_swe 
	@gnuplot -p -e "set key autotitle columnhead;pl for [col=2:6] 'wf_data.txt' u 1:col w l lw 3"

rfsc:
	${fort_comp} -o test_rfsc numerical_integration.f90 rfscat_cart.f90

rfsc_plot: rfsc
	./rfsc_cart.sh
	gnuplot -p rfsc_cart.gp

rfsc_pol:
	${fort_comp} -o rfsc_pol numerical_integration.f90 rfsc_pol.f90

rfsc_pol_plot:	rfsc_pol
	./rfsc_pol.sh
	gnuplot -p rfsc_pol.gp
howf:
	${fort_comp} -o howf_test numerical_integration.f90 ho_wf.f90

howf_plot: howf
	./howf_test
	@gnuplot -p -e "pl for [col=2:*] 'howf_data.txt' u 1:col w l lw 3"
	@gnuplot -p -e "pl 'ham_data.txt' matrix w image"
