plots = reactor_spectra.pdf
data = sinev_al.txt mueller_al.txt 

all: $(plots)

sinev_al: % : %.f90
# static pattern rule
	@gfortran -o $@ $^
	@echo $@ done.

mueller_al: % : %.c %.h
	@gcc -o $@ -lm $^
	@echo $@ done.

$(data) : %.txt : %
	@./$< > $@
	@echo $@ done.

reactor_spectra.pdf : %.pdf : %.tex $(data)
	@pdflatex $<

