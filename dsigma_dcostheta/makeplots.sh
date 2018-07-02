export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:..
python dsigma_dcostheta.py
echo "dsigma_dcostheta.tex" | pdflatex -jobname=dsigma_dcostheta main.tex
