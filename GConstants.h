#ifndef G_CONSTANTS_H
#define G_CONSTANTS_H

#include <cmath>


enum NucleonType
{
    proton,
	neutron
};

enum NeutrinoFlav
{
	flav_e,
	flav_mu,
	flav_tau,
	flav_mless
};


const double MeV( 1.0 );
const double MeV2( MeV*MeV );
const double GeV( 1000.0*MeV );
const double GeV2( GeV*GeV );

const double fm( 1.0/197.3269718*MeV );//PDG 2014, p. 109
const double cm( 1.0e+13*fm );
const double cm2( cm*cm );
const double sr( 1.0 );

const double pi(4*atan(1.0));
const double pi2(pi*pi);

const double GF( 1.1663787e-5/GeV2 );//Fermi constant, PDG 2014, p. 109
const double cosThetaC( 0.97425 );//Cabibbo angle, PDG 2014, p. 214
const double GFcosTheta( GF*cosThetaC );//Fermi constant including the Cabibbo angle

const double alpha( 7.2973525698e-3 );//fine structure constant, PDG 2014, p. 109

const double pMass( 938.272046*MeV );//proton mass, PDG 2014, p. 80
const double nMass( 939.565379*MeV );//neutron mass, PDG 2014, p. 80
const double M( 0.50*(pMass+nMass) );

const double piCMass( 139.57018*MeV );//charged pion mass, PDG 2014, p. 34
const double pi0Mass( 134.9766*MeV );//neutral pion mass, PDG 2014, p. 34

const double eMass(     0.510998928*MeV );//electron mass, PDG 2014, p. 30
const double muMass(  105.6583715*MeV );//muon mass, PDG 2014, p. 30
const double tauMass(1776.82*MeV );//tau mass, PDG 2014, p. 30

const double pMass2( pMass*pMass );
const double nMass2( nMass*nMass );
const double M2( M*M );
const double piCMass2( piCMass*piCMass );

const double mu_p( 2.792847356 );//PDG 2014, p. 80
const double mu_n(-1.9130427 );//PDG 2014, p. 80
const double gA( -1.2723 );//PDG 2014, p. 1384
const double MA( 1026.0*MeV );

const double MA2( MA*MA );
const double MV2( 0.71*GeV2 );

#endif
