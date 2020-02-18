default:
	ifort -o test_swe numerical_integration.f90 swe_inf_pot_well.f90

plot:	default
	./test_swe > test_swe.txt
	gnuplot -e "set key autotitle columnhead;pl 'test_swe.txt' u 1:2; pause -1"
