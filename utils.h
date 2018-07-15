#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

typedef std::vector<std::pair<double, double> > spectrum_vector; // spectrum<Energy, NEvents>, can be used for other functions.

void print_spectrum(const char *foutname, const char *header, spectrum_vector *spectr, const char *sep = " ") {
    ofstream fout(foutname);
    fout << header << endl;
    for(int i = 0; i < spectr->size(); ++i) {
        fout << (*spectr)[i].first << sep << (*spectr)[i].second << endl;
    }
}

#endif

