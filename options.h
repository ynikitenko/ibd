#ifndef OPTIONS_H
#define OPTIONS_H
#include "GConstants.h"

///--------------options-------------------------------------------------------
extern "C" {
    bool switchAntineutrino( true );///true = antineutrino, false = neutrino
    NeutrinoFlav flav( flav_e );///set the neutrino flavor, flav_e = electron, flav_mu = muon, flav_tau = tau, flav_mless = massless
    bool switchMassDifference( true );///take into account the neutron-proton mass difference, true = yes, false = no
    bool switchCVC( true );///conservation of the vector current, true = Ankowski, false = Strumia&Vissani
    bool switchRC( true );///radiative corrections, only for electron neutrino, true = yes, false = no

    double neutrinoE( 2.0*MeV );///neutrino energy
    int cosThetaRes(200);///cosTheta resolution
}
// on extern keyword: https://en.cppreference.com/w/cpp/language/language_linkage

#endif

