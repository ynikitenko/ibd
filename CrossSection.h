#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#include "GConstants.h"
#include "options.h"
#include "FormFactors.h"
#include "SpecialFunctions.h"

extern "C" {
    double find_LH_free(const double& eK, const double& qSq);
    double find_DCS_free_cosTheta(const double& eK, const double& cosTheta);
    double find_DCS_free_cosTheta_RC(const double eK, const double cosTheta);
}

///----------------------------------------------------------------------------

void print_dsigma_dcostheta(const char *outfilename) {
    std::ofstream result(outfilename);///output-file name
    const double cosThetaStep( 2.0/cosThetaRes );
    double cosTheta( -1.0 + 0.5*cosThetaStep );
    double dcs(0.0);

    for (int lCnt(0); lCnt < cosThetaRes; ++lCnt)
    {
        dcs = switchRC ? find_DCS_free_cosTheta_RC(neutrinoE, cosTheta) : find_DCS_free_cosTheta(neutrinoE, cosTheta);

        result<<std::showpoint<<std::fixed<<std::setw(8)<<std::setprecision(4)<<cosTheta<<" "
                <<std::showpoint<<std::scientific<<std::setw(12)<<std::setprecision(4)<<dcs<<std::endl;

        cosTheta += cosThetaStep;
    }
}

// options for cross-section output
extern double neutrinoEMin( 1.8*MeV );///minimal neutrino energy
extern double neutrinoEMax( 20.0*MeV );///maximal neutrino energy
extern int neutrinoERes(180);///energy resolution
// std::ofstream result( "cs.txt" );///output-file name

double get_sigma(double neutrinoE) {
    // get total IBD cross section in 10^{-38} cm^2.
    // neutrinoE is in MeV.
    double cs(0.0);

    // integrate on \cos\theta
    const double cosThetaStep( 2.0/cosThetaRes );
    for (int lCnt(0); lCnt < cosThetaRes; ++lCnt)
    {
        const double cosTheta( -1.0 + (0.5 + lCnt)*cosThetaStep );
        const double dcs = switchRC ? find_DCS_free_cosTheta_RC(neutrinoE, cosTheta) : find_DCS_free_cosTheta(neutrinoE, cosTheta);

        cs += dcs;
    }

    cs *= cosThetaStep;
    return cs;
}

void print_sigma(const char *outfilename) {
    std::ofstream result(outfilename);///output-file name
    const double neutrinoEStep( (neutrinoEMax - neutrinoEMin)/neutrinoERes );

    for (int nCnt(0); nCnt <= neutrinoERes; ++nCnt )
    {
        double neutrinoE(neutrinoEMin + nCnt*neutrinoEStep);
        double cs = get_sigma(neutrinoE);
        result<<std::showpoint<<std::fixed<<std::setw(5)<<std::setprecision(2)<<neutrinoE<<" "<<std::showpoint<<std::scientific<<std::setw(14)<<std::setprecision(5)<<cs<<std::endl;
    }
}

const double leptMass( (flav == flav_e) ? eMass : (flav == flav_mu) ? muMass : (flav == flav_tau) ? tauMass : 0.0 );
const double leptMassSq( leptMass*leptMass );
const double initNuclMass(   not switchMassDifference ? M : switchAntineutrino ? pMass : nMass   );
const double finalNuclMass(  not switchMassDifference ? M : switchAntineutrino ? nMass : pMass   );

const double initNuclMassSq( initNuclMass*initNuclMass );
const double finalNuclMassSq( finalNuclMass*finalNuclMass );


double find_LH_free(const double& eK, const double& qSq)
// The function calculating the contraction of the lepton and hadron tensors, corresponding to the expresion in the square brackets of Eq. (12),
// that is M^4 A + M^2 B (s-u) + C (s- u)^2.
// eK is the value of (anti)neutrino energy in MeV and qSq is the value of 4-momentum transfer squared, q^2, given in MeV^2.
{
    static const double avNuclMass(   M   );
    static const double avNuclMassSq(   M2   );

    static const double deltaM( finalNuclMass - initNuclMass );
    static const double deltaDiv( deltaM/avNuclMass );
    static const double mu( 0.25*leptMassSq/avNuclMassSq );

    const double p_o_k( initNuclMass*eK );

    const double tau( -0.25*qSq/avNuclMassSq );

    const double geP( ProtonKelly_GE(qSq) );
    const double geN( NeutronKelly_GE(qSq) );
    const double gmP( ProtonKelly_GM(qSq) );
    const double gmN( NeutronKelly_GM(qSq) );

    const double fpDivFA( 2.0*avNuclMassSq/(piCMass2 - qSq) );///pseudoscalar FF divided by the axial one

    const double gAxial( FA(qSq) );///axial form factor

    const double ge( geP - geN );
    const double gm( gmP - gmN );

    const double fa( gAxial );

    const double f1 = (ge + tau*gm)/(1 + tau);
    const double f2 = (gm - ge)/(1 + tau);

    const double gmSq = gm*gm;
    const double geSq = ge*ge;
    const double ff = (geSq + tau*gmSq)/(1 + tau);///( f1*f1 + f2*f2*tau );
    const double faGM = fa*gm;
    const double faSq = fa*fa;

    const double s_u( 4*p_o_k + qSq - leptMassSq - 2*avNuclMass*deltaM );

    const double z( switchCVC ? 1.0 : 0.0 );

    const double invA(   4*(tau + mu)*( (tau - mu)*(gmSq + faSq) - geSq + faSq + 4*mu*faSq*fpDivFA*(tau*fpDivFA - 1.0)  )
                        -8*mu*deltaDiv*faGM + deltaDiv*deltaDiv*(tau + mu)*( gmSq - (1.0 + tau)*f2*f2 - faSq + 4*mu*faSq*fpDivFA*fpDivFA )
                        + deltaDiv*deltaDiv*( tau*gmSq - (1.0 + tau)*f1*f1 - faSq - 4*mu*faSq*fpDivFA + 2*z*mu*f1*f2*leptMassSq/qSq )
                        + z*deltaDiv*deltaDiv*( leptMassSq/qSq*(1 + leptMassSq/qSq) - mu*(1 - leptMassSq/qSq) )*f1*f1 );
    const double invB(   4*tau*faGM - mu*deltaDiv*(f2*f2 + (1.0 - z)*f1*f2 + z*f1*f1/tau + 2*faSq*fpDivFA)    );
    const double invC(   0.25*( faSq + ff )   );

    return switchAntineutrino ? avNuclMassSq*avNuclMassSq*invA + avNuclMassSq*invB*s_u + invC*s_u*s_u : avNuclMassSq*avNuclMassSq*invA - avNuclMassSq*invB*s_u + invC*s_u*s_u;
}



double find_DCS_free_cosTheta(const double& eK, const double& cosTheta)
// The differential cross section as a function of cos(theta), given in unit 10^-38 cm^2,
// calculated without radiative corrections.
// eK is the value of (anti)neutrino energy in MeV and cosTheta is the value of cos(theta), theta being the charged lepton production angle.
{
    static const double coeffCC(   GFcosTheta*GFcosTheta/(8*pi*initNuclMassSq)/(1.0e-38*cm2/MeV2)   );
    static const double mmm(   0.5*leptMassSq + 0.5*initNuclMassSq - 0.5*finalNuclMassSq   );

    const double kcos( eK*cosTheta );
    const double eKM( eK + initNuclMass );
    const double dummy_MxK_p_m(  initNuclMass*eK + mmm );
    const double dummy_KxK_m_m( eK*eK - mmm );
    const double denomin( eKM*eKM - kcos*kcos );
    const double root(   sqrt( dummy_MxK_p_m*dummy_MxK_p_m - leptMassSq*denomin )   );
    const double omega(   ( eKM*dummy_KxK_m_m - eK*kcos*kcos - kcos*root )/denomin   );

    if (omega > eK - leptMass or omega < 0.0)
        return 0.0;

    const double qSq( -2*initNuclMass*omega - initNuclMassSq + finalNuclMassSq );
    const double kPrimeNorm(   sqrt( std::pow(eK - omega, 2) - leptMassSq )   );
    const double jacobian = std::abs(   2*initNuclMass*eK*std::pow(kPrimeNorm, 2)/( (eK - omega)*kcos - kPrimeNorm*eKM )   );


    const double coefE( initNuclMass + 2*eK );
    const double coefD( initNuclMass*coefE - finalNuclMassSq - leptMassSq );
    const double delta(   ( coefD*coefD - 4*leptMassSq*finalNuclMassSq )*eK*eK   );

    if ( delta <= 0.0 )
        return 0.0;

    const double dummyA(   leptMassSq*initNuclMass - eK*coefD   );
    const double dummyB(   sqrt(delta)   );

    const double qSqMin(   ( dummyA - dummyB )/coefE   );
    const double qSqMax(   ( dummyA + dummyB )/coefE   );///close to zero

    return ( qSqMin <= qSq  and qSq <= qSqMax ) ? jacobian*coeffCC/std::pow(eK, 2)*find_LH_free(eK, qSq) : 0.0;
}


double find_DCS_free_cosTheta_RC(const double eK, const double abscissaValue)
// The differential cross section as a function of cos(theta), given in unit 10^-38 cm^2,
// calculated with radiative corrections.
// The cos(theta) -even and cos(theta)-odd parts of the cross section are treated separately,
// as they are subject to different outer radiative corrections.
// eK is the value of (anti)neutrino energy in MeV 
// abscissaValue is the value of cos(theta), theta being the charged lepton production angle.
{
    if ( not (flav == flav_e) )
        return find_DCS_free_cosTheta( eK, abscissaValue );

    const double cosPTheta(   find_DCS_free_cosTheta( eK, abscissaValue )   );
    const double cosMTheta(   find_DCS_free_cosTheta( eK, -abscissaValue )   );

    if ( cosPTheta <= 0.0 and cosMTheta <= 0.0  )
        return 0.0;

    static const double twoMDelta( finalNuclMassSq - initNuclMassSq );
    const double coef_ASq( 2*eK*initNuclMass - twoMDelta + leptMassSq );
    const double coef_B( 2*(eK + initNuclMass) );
    const double coef_C( 2*eK*abscissaValue );
    const double denomin( coef_B*coef_B - coef_C*coef_C );

    const double rootSq(   coef_ASq*coef_ASq - leptMassSq*denomin   );

    if ( rootSq < 0.0  )
        return 0.0;

    const double eKPrime(   ( coef_ASq*coef_B + coef_C*sqrt(rootSq) )/denomin   );

    const double betaSq( 1.0 - std::pow(eMass/eKPrime, 2) );
    const double beta( sqrt(betaSq) );
    const double sqrtP(   sqrt(1.0 + beta)   );
    const double sqrtM(   sqrt(1.0 - beta)   );
    const double dummy(   (1.0 + beta)/(1.0 - beta)   );
    const double lne(   log(dummy)   );
    static const double lnConst(   3*log(pMass/eMass) + 23.0/4   );

    const double deltaOut(  lnConst + 8.0/beta*spence( 2*beta/(1 + beta) ) - 2.0/beta*lne*lne
                       + 4*log( 4*betaSq/(1 - betaSq) )*( 0.5/beta*lne - 1.0 ) + (0.75*beta + 1.75/beta)*lne );

    const double deltaTilde(   3*log(pMass/eMass) + 3.0/4 + 4.0/betaSq*( 1.0 - sqrt(1.0 - betaSq) ) + 8.0/beta*spence( 1.0 - sqrtM/sqrtP ) + ( 0.5/beta - 0.375 - 0.125/betaSq )*lne*lne
                            + (0.5/beta - 2.0)*lne - 4*(0.5/beta*lne - 1.0)*log( 0.5/beta*(1.0 + beta)*( sqrtP + sqrtM )/( sqrtP - sqrtM ) )   );

    const double radCorrection1( 1.0 + 0.02250 + alpha/(2*pi)*deltaOut );
    const double radCorrection2( 1.0 + 0.02250 + alpha/(2*pi)*deltaTilde );

    return 0.5*( cosPTheta + cosMTheta )*radCorrection1 + 0.5*( cosPTheta - cosMTheta )*radCorrection2;
}

