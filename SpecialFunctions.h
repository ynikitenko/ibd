#ifndef SPECIAL_FUNCTIONS_H
#define SPECIAL_FUNCTIONS_H

#include <cmath>

///accuracy verified with the book of Abramowitz and Stegun, p. 1005
///over the interval [0.5:1.0] every 0.01, the differences are below 1.0e-9 (rounding error)
/// spence(x) = - f^{AS}(1 - x) = - Li_2(x)
double spence(const double x)
{
    if ( x > 1.0 )
        return 0.0;

    if ( std::abs(x) < 1.0e-1 )
        return -x*(     1.0 + x*(    0.25 + x*(   1.0/9.0 + x*(  0.0625 + x*( 0.04 + x*(1.0/36.0 + x*(1.0/49.0 + x*0.015625)) )  )   )    )     );
        /// -( x + x^2/2^2 + x^3/3^2 + x^4/4^2 + ... + x^n/n^2 )

    if ( std::abs(1 - x) < 1.0e-10 )
        return -M_PI*M_PI/6.0 - spence(1 - x);

    if ( x > 0.5 )
        return -M_PI*M_PI/6.0 + log(x)*log(1 - x) - spence(1 - x);


    static const double xWeights[] = {0.1488743389, 0.4333953941, 0.6794095682, 0.8650633666, 0.9739065285};
    static const double fWeights[] = {0.2955242247, 0.2692667193, 0.2190863625, 0.1494513491, 0.0666713443};

    const double xr( 0.5*x );
    double sum(0.0), xm(0.0);

    for (int iCnt(0); iCnt < 5; ++iCnt)
    {
        xm = xr*xWeights[iCnt];

        sum += fWeights[iCnt]*( log(1.0 - xr - xm)/(xr + xm) + log(1.0 - xr + xm)/(xr - xm) );
    }

    return sum*0.5*x;

}


double j0(const double x)
{
    if (x < 0.4)
    {
        const double xSq(x*x);
        return 1.0 - xSq/6*(  1.0 - xSq/20*(1.0 - xSq/42*(1.0 - xSq/72*(1.0 - xSq/110)))   );
    }

    return sin(x)/x;
}

double slater(const double x)
{
    if (x < 0.4)
    {
        const double xSq(x*x);
        return 1.0 - 0.1*xSq*(  1.0 - xSq/28*(1.0 - xSq/54*(1.0 - xSq/88*(1.0 - xSq/130)))   );
    }

    return 3.0*( sin(x) - x*cos(x) )/(x*x*x);
}

//slater - j0
double L(const double x)
{
    if (x < 0.4)
    {
        const double xSq(x*x);
        return xSq/15*(  1.0 - xSq/14*(1.0 - xSq/36*(1.0 - xSq/66*(1.0 - xSq/104)))  );
    }

    return slater(x) - j0(x);
}

//slater - j0
double L0(const double x)
{
    return L(x);
}

//L0 - x^2/15*slater
double L1(const double x)
{
    const double xSq(x*x);

    if (x < 0.4)
        return xSq*xSq/525*(  1.0 - xSq/18*(1.0 - xSq/44*(1.0 - xSq/78))  );

    return L(x) - xSq/15.0*slater(x);
}

//L0 - x^2/15*slater - x^2/35*L
double L2(const double x)
{
    const double xSq(x*x);

    if (x < 0.4)
        return xSq*xSq*xSq/33075*(  1.0 - xSq/22*(1.0 - xSq/52)  );

    return L(x)*(1.0 - xSq/35.0) - xSq/15.0*slater(x);
}

//L0 - x^2/15*slater - 2*x^2/45*L + x^4/945*l
double L3(const double x)
{
    const double xSq(x*x);

    if (x < 0.4)
        return xSq*xSq*xSq*xSq/3274425*(1.0 - xSq/26);

    return L(x)*(1.0 - 2*xSq/45.0) - xSq/15.0*( 1.0 - xSq/63.0 )*slater(x);
}

double U(const double r)
{
    static const double kF( 1.333/fm );
    static const double eF( sqrt( M2 + kF*kF ) );

    static const double a1( M/( M + eF ) );
    static const double a2( M*(5*M + 6*eF)/pow(M + eF, 2) );
    static const double a3( M*(15*M2 + 34*M*eF + 20*eF*eF)/pow(M + eF, 3) );

    const double x( kF*r );
    const double reFSq( pow(r*eF, 2) );

    return slater(x) + (   1.5*a1*L0(x) + ( 3.75*a2*L1(x) + 39.375*a3*L2(x)/reFSq )/reFSq   )/reFSq;
}

double T(const double r)
{
    static const double kF( 1.333/fm );
    static const double eF( sqrt( M2 + kF*kF ) );

    static const double a1(   ( M + 2*eF )/( M + eF )   );
    static const double a2(   ( 3*M2 + 6*M*eF + 4*eF*eF )/pow(M + eF, 2)   );
    static const double a3(   ( 7*M2*M + 20*M*eF*(M + eF) + 8*eF*eF*eF )/pow(M + eF, 3)   );

    const double x( kF*r );
    const double reFSq( pow(r*eF, 2) );

    return slater(x) - (   1.5*a1*L0(x) + ( 3.75*a2*L1(x) + 39.375*a3*L2(x)/reFSq )/reFSq   )/reFSq;
}

double I(const double r)
{
    static const double kF( 1.333/fm );
    static const double eF( sqrt( M2 + kF*kF ) );

    static const double a1(   ( M + 2*eF )/( M + eF )   );
    static const double a2(   ( 5*M2 + 14*M*eF + 12*eF*eF )/pow(M + eF, 2)   );

    const double x( kF*r );
    const double reFSq( pow(r*eF, 2) );

    return L0(x) + ( 2.5*a1*L1(x) + 8.75*a2*L2(x)/reFSq )/reFSq;
}

double J(const double r)
{
    static const double kF( 1.333/fm );
    static const double eF( sqrt( M2 + kF*kF ) );

    static const double a1(   M/( M + eF )   );
    static const double a2(   M*( 5*M + 6*eF )/pow(M + eF, 2)   );

    const double x( kF*r );
    const double reFSq( pow(r*eF, 2) );

    return L0(x) + ( 2.5*a1*L1(x) + 8.75*a2*L2(x)/reFSq )/reFSq;
}

double K(const double r)
{
    static const double kF( 1.333/fm );
    static const double eF( sqrt( M2 + kF*kF ) );

    static const double a1(   M/( M + eF )   );
    static const double a2(   M*( M + 2*eF )/pow(M + eF, 2)   );

    const double x( kF*r );
    const double reFSq( pow(r*eF, 2) );

    return L0(x) - ( 2.5*a1*L1(x) + 26.25*a2*L2(x)/reFSq )/reFSq;
}

double R(const double r)
{
    static const double kF( 1.333/fm );
    static const double eF( sqrt( M2 + kF*kF ) );

    static const double a1(   ( M + 2*eF )/( M + eF )   );
    static const double a2(   ( 5*M2 + 14*M*eF + 12*eF*eF )/pow(M + eF, 2)   );

    const double x( kF*r );
    const double reFSq( pow(r*eF, 2) );

    return 5*L1(x) + ( 17.5*a1*L2(x) + 78.75*a2*L3(x)/reFSq )/reFSq;
}

#endif
