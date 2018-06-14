#ifndef FORM_FACTORS_H
#define FORM_FACTORS_H

#include "GConstants.h"


double FA(const double qSq)///axial form factor
{
	return gA/std::pow(1.0 - qSq/MA2, 2);
}


double FP(const double qSq)///pseudoscalar form factor
{
	return 2.0*M2*FA(qSq)/(piCMass2 - qSq);
}


double ProtonKelly_GE(const double qSq)///proton electic form factor
{
    ///J.J. Kelly, Phys. Rev. C 70, 068202 (2004)
    const double tau( -0.25*qSq/pMass2 );
    const double gEp(     (1.0 - 0.24*tau)/(  1.0 + ( 10.98 + (12.82 + 21.97*tau)*tau )*tau  )     );

	return gEp;
}


double ProtonKelly_GM(const double qSq)///proton magnetic form factor
{
    ///J.J. Kelly, Phys. Rev. C 70, 068202 (2004)
    const double tau( -0.25*qSq/pMass2 );
    const double gMp(     (1.0 + 0.12*tau)/(  1.0 + ( 10.97 + (18.86 + 6.55*tau)*tau )*tau  )     );

	return mu_p*gMp;
}


double NeutronKelly_GE(const double qSq)///neutron electic form factor
{
    ///S. Riordan et al., Phys. Rev. Lett 105, 262302 (2010)
	const double tau( -0.25*qSq/nMass2 );
	const double gEn(   1.39*tau/( 1.0 + 2.00*tau )/std::pow(1 - qSq/MV2, 2)   );

    //J.J. Kelly, Phys. Rev. C 70, 068202 (2004)
    //const double gEn(   1.70*tau/( 1.0 + 3.30*tau )/std::pow(1 - qSq/MV2, 2)   );

	return gEn;
}


double NeutronKelly_GM(const double qSq)///neutron magnetic form factor
{
    ///J.J. Kelly, Phys. Rev. C 70, 068202 (2004)
    const double tau( -0.25*qSq/nMass2 );
    const double gMn(     (1.0 + 2.33*tau)/(  1.0 + ( 14.72 + (24.20 + 84.1*tau)*tau )*tau  )     );

	return mu_n*gMn;
}


#endif


