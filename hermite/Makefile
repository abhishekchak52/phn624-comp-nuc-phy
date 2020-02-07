fort_comp = ifort
simp_test_src = test_simpson.f90 simpson.f90
simp_wf_test_src = test_simpson_wf.f90 simpson.f90
ortho_test_src = simpson.f90 ho_wf.f90 ho_wf_ortho_test.f90

ho_wf.mod: ho_wf.f90
	$(fort_comp) -c ho_wf.f90 -o ../obj/ho_wf.o

ortho_test: $(ortho_test_src)
	$(fort_comp) -o ortho_test $(ortho_test_src)

test_simpson: $(simp_test_src)
	$(fort_comp) -o test_simpson $(simp_test_src)

test_simpson_wf: $(simp_wf_test_src)
	$(fort_comp) -o test_simpson_wf $(simp_wf_test_src)
