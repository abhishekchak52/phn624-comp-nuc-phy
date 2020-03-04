rfsc:
	ifort -no-wrap-margin -o rfsc_mte rfsc_mte.f90
rfsc_mte: rfsc
	./rfsc_mte.sh
	@gnuplot -p rfsc_mte.gp
	latex rfsc_mte.tex
	dvipdf rfsc_mte.dvi
