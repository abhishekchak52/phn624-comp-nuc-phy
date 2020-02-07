#makefile sets option -mkl
Comp = ifort
sources = dsyevr_ex.f
FFLAGS = -mkl
test : $(sources)
	$(Comp) -o test $(sources) ${FFLAGS}
