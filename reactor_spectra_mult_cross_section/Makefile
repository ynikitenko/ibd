prog := reactor_mult_cs
data := reactor_mueller_mult_cs_cvc.txt reactor_mueller_mult_cs_no_cvc.txt reactor_mult_cs_cvc_div_no_cvc.txt
plots := $(prog).pdf reactor_mult_cs_cvc_div_no_cvc.pdf
all : $(plots)

$(plots) : %.pdf : %.tex $(data)
	@pdflatex $<
	@echo $@ done.

#$(prog).txt : %.txt : %
%.txt : $(prog)
	@./$< 
	@echo $@ done.

$(prog) : % : %.cpp 
	@g++ -O3 -o $@ $<
	@echo $@ done.

