fort_comp = ifort
simp_test_src = test_simpson.f90 simpson.f90

test_simpson: $(simp_test_src)
	$(fort_comp) -o test_simpson $(simp_test_src)
