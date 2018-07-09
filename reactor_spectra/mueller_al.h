#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Mueller et al., Improved Predictions of Reactor Antineutrino Spectra
// Phys.Rev.C 83:054615, 2011. 
// http://lanl.arxiv.org/abs/1101.2663v3

enum ISO {
    U235 = 0,
    U238,
    Pu239,
    Pu241
};
//const short NISO = 4; // number of isotopes
#define NISO 4 // number of isotopes

// Options
double iso_fraction[NISO] = {0.07, 0.3, 0.03, 0.6}; // isotope fractions in the whole flux

double E_nu_min = 1.8; // minimum antineutrino energy in MeV
double E_nu_max = 9.25; // maximum antineutrino energy in MeV
double step = 0.01;
// int nbins = 1200;       // precision of the spectrum. 
// E_nu_max is evaluated. npoints = nbins + 1.

// internal array for fixed parameters
const double *sp_coeff[NISO];
const int NPT = 6;
char spectra_initialized = 0;

int init_spectra() {
    // polynomial coefficients from the article.
    static const double ar_U235[]  = { 3.217, -3.111, 1.395, -3.690E-1, 4.445E-2, -2.053E-3 };
    static const double ar_U238[]  = { 4.833E-1, 1.927E-1, -1.283E-1, -6.762E-3, 2.233E-3, -1.536E-4 };
    static const double ar_Pu239[] = { 6.413, -7.432, 3.535, -8.820E-1, 1.025E-1, -4.550E-3 };
    static const double ar_Pu241[] = { 3.251, -3.204, 1.428, -3.675E-1, 4.254E-2, -1.896E-3 };
    sp_coeff[U238]  = ar_U238; 
    sp_coeff[U235]  = ar_U235;
    sp_coeff[Pu239] = ar_Pu239;
    sp_coeff[Pu241] = ar_Pu241;
    // memcpy(sp_coeff[U238], ar_U238, NPT * sizeof(double));
    spectra_initialized = 1; 
    return 0;
}

double get_spectrum_iso(double E_nu, enum ISO iso) {
    if (!spectra_initialized)
        init_spectra();
    double res = 0;
    for (int i = 0; i < NPT; ++i) {
        res += sp_coeff[iso][i] * pow(E_nu, i);
        // printf("i = %d, sp %f, pow %f for iso = %d\n", i, sp_coeff[iso][i], pow(E_nu, i), iso);
    }
    // printf("%f, %f\n", E_nu, res);
    return exp(res);
}

double get_spectrum(double E_nu) {
    double res = 0;
    for (enum ISO iso = 0; iso < NISO; ++iso)
        res += iso_fraction[iso] * get_spectrum_iso(E_nu, iso);
    return res;
}

double *get_spectrum_arr() {
    printf("# Reactor spectrum for antineutrino energy\n");
    int nbins = (E_nu_max - E_nu_min) / step; 
    // double step = (E_nu_max - E_nu_min) / nbins;
    double *spectrum = malloc((nbins + 1) * sizeof(double));
    int i = 0;
    for (double E_nu = E_nu_min; E_nu <= E_nu_max; E_nu += step) {
    // for (int i = 0; i <= nbins; ++i) {
    //     double E_nu = E_nu_min + step * i;
        double sp = get_spectrum(E_nu);
        spectrum[i] = sp;
        printf("%f %f\n", E_nu, sp);
        ++i; 
    }
    return spectrum;
}

