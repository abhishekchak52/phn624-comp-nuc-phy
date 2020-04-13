fort_comp = ifort -no-wrap-margin
lib = numerical_constants.f90 numerical_integration.f90 nuclear_models.f90
extra_libs = cgc.f diag.f factlog.f sort.f
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
	@gnuplot -p -e "pl for [col=2:10] 'howf_data.txt' u 1:col w l lw 3"
	@gnuplot -p -e "pl 'ham_data.txt' matrix w image"

alpa:
	${fort_comp} -o alpt special_functions.f90 assoc_leg_pol_test.f90

test_alpa: alpa
	./alpt > alpt_data.txt
	gnuplot -p -e "pl for [col=2:6] 'alpt_data.txt' u 1:col w l lw 3 "

spha:
	${fort_comp} -o spharmonics special_functions.f90 spherical_harmonics.f90

plot_spha: spha
	./spharmonics > spha_data.txt
	# gnuplot -p -e "set mapping spherical; spl 'spha_data.txt' u 2:1:3 w l"
	# gnuplot -p -e "set mapping spherical; spl 'spha_data.txt' u 2:1:4 w l"
	gnuplot -e "set mapping spherical; spl 'spha_data.txt' u 2:1:5; pause -1"

multipole_mom:
	${fort_comp} -o mul_mom numerical_constants.f90 special_functions.f90 numerical_integration.f90 multipole_moment.f90

tov_src:
	${fort_comp} -o tov numerical_constants.f90 numerical_integration.f90 tov.f90

tov_plot: tov_src
	./tov

nilsson_src:
	${fort_comp} -o nilsson ${extra_libs} ${lib} nilsson.f90 -mkl

nilsson_plot: nilsson_src
	./nilsson

bcs_src:
	${fort_comp} -o nil_bcs bcs.f nil_bcs.f ${extra_libs} -mkl

bcs_plots: bcs_src
	./nil_bcs
	@gnuplot -p bcs.gp
hot_src:
	${fort_comp} -o nil_therm secant.f nil_therm.f ${extra_libs} -mkl
hot_plots: hot_src
	./nil_therm
#	@gnuplot -p -e "pl 'hot_nuc.txt' u 0:2 w l"
#	 @gnuplot -p -e "pl 'ef_delta.txt' u 0:1 w l,'' u 0:2 w l"
# @gnuplot -p -e "pl 'chempot_temp.txt' u 1:2 w l"
#	 @gnuplot -p -e "set term png;set output 'f_dist.png';pl 's_data.txt' u 2:3 w l,'' u 2:4 w l"
#	 @gnuplot -p -e "set term png;set output 's_dist.png';pl 's_data.txt' u 2:5 w l,'' u 2:6 w l"
	@gnuplot -p hot.gp

hot_aru_src:
	${fort_comp} -o nil_hot secant.f nil_hot.f ${extra_libs} -mkl
hot_aru: hot_aru_src
	./nil_hot
