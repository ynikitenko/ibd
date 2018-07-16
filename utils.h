#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
using namespace std;

typedef std::vector<std::pair<double, double> > spectrum_vector; // spectrum<Energy, NEvents>, can be used for other functions.

void print_spectrum(const char *foutname, const char *header, spectrum_vector *spectr, const char *sep = " ") {
    ofstream fout(foutname);
    fout << header << endl;
    for(int i = 0; i < spectr->size(); ++i) {
        fout << (*spectr)[i].first << sep << (*spectr)[i].second << endl;
    }
}

spectrum_vector *spectrum_div(const spectrum_vector *num_sp, const spectrum_vector *denom_sp) {
    // check that two spectra scale and size correspond
    int nbins = num_sp->size();
    if (denom_sp->size() != nbins) 
        throw std::length_error("spectrum vectors should have same size.");
    if (nbins)
        if ((*num_sp)[0].first != (*denom_sp)[0].first
            || (*num_sp)[nbins - 1].first != (*denom_sp)[nbins - 1].first)
                throw std::logic_error("spectrum vectors should have same scale.");

    // skip interval where denominator spectrum is zero
    int nzeroes = 0;
    while ((*denom_sp)[nzeroes].second == 0) 
        ++nzeroes;
    nbins -= nzeroes;

    spectrum_vector *sp = new spectrum_vector(nbins);
    int nnonzeroes = 0;
    // divide spectra while denominator spectrum is nonzero
    for (;nnonzeroes < nbins; ++nnonzeroes) {
        int i = nnonzeroes + nzeroes;
        double denom_sp_val = (*denom_sp)[i].second;
        if (denom_sp_val == 0)
            break;
        (*sp)[nnonzeroes].first = (*num_sp)[i].first;
        (*sp)[nnonzeroes].second =(*num_sp)[i].second / denom_sp_val;
    }
    sp->resize(nnonzeroes);
    return sp;
}

#endif

