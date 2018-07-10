all: nu_e_nu_e.pdf

nu_e_nu_e : % : %.c %.h
	@gcc -o $@ -lm $<
	@echo $@ done.

# example target to run the program. Other .txt files will be done automatically.
data := nu_e_nu_e_sin2_theta13_dm2m0.txt
#data := $(wildcard *.txt)

$(data): nu_e_nu_e
	@./$<
	@echo $@ done.

nu_e_nu_e.pdf : %.pdf : %.tex % $(data)
	@pdflatex $<
	@echo $@ done.

