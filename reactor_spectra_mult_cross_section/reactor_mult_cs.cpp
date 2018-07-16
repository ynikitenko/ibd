#include <iostream>
#include <fstream>
#include "../utils.h"
#include "../reactor_spectra/mueller_al.h"
#include "../CrossSection.h"
using namespace std;

spectrum_vector *reactor_mult_cs_spectrum(const spectrum_vector *reactor_spectr) {
    int nbins = reactor_spectr->size();
    spectrum_vector *sp = new spectrum_vector(nbins);
    for(int i = 0; i < nbins; ++i) {
        double sp_val = (*reactor_spectr)[i].second;
        double Enu = (*reactor_spectr)[i].first;
        (*sp)[i].first = Enu;
        (*sp)[i].second = sp_val * get_sigma(Enu);
    }
    return sp;
}

void print_reactor_mult_cs_spectrum() {
    E_nu_min  = 1.8; 
    // E_nu_max  = 4; 
    E_nu_max  = 10; 
    step = 0.1;
    // Mueller et al reactor spectrum
    spectrum_vector *init_spectrum = get_spectrum_vector();

    switchCVC = true;
    spectrum_vector *react_cvc = reactor_mult_cs_spectrum(init_spectrum);
    switchCVC = false;
    spectrum_vector *react_no_cvc = reactor_mult_cs_spectrum(init_spectrum);

    const char *outfilename_cvc = "reactor_mueller_mult_cs_cvc.txt";
    print_spectrum(outfilename_cvc, "# Reactor spectrum convolved with cross section at the given neutrino energy", react_cvc);
    const char *outfilename_no_cvc = "reactor_mueller_mult_cs_no_cvc.txt";
    print_spectrum(outfilename_no_cvc, "# Reactor spectrum convolved with cross section at the given neutrino energy", react_no_cvc);

    const char *outfilename_sp_div = "reactor_mult_cs_cvc_div_no_cvc.txt";
    spectrum_vector *div_sp = spectrum_div(react_cvc, react_no_cvc);
    print_spectrum(outfilename_sp_div, "# Reactor spectrum convolved with cross section: CVC over non-CVC at the given neutrino energy", div_sp);

    delete react_cvc, react_no_cvc;
    delete div_sp;
    delete init_spectrum;
}

int main() {
    print_reactor_mult_cs_spectrum();
    return 0;
}

