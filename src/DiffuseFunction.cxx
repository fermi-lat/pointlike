/** @file DiffuseFunction.cxx

$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/tools/src/DiffuseFunction.cxx,v 1.4 2006/08/15 21:26:27 burnett Exp $
*/

#include "pointlike/DiffuseFunction.h"
#include <cmath>

using namespace pointlike;

double DiffuseFunction::s_emin(10.);

double DiffuseFunction::energy_bin(int k)
{
    return s_emin*pow(2., k);
}

double DiffuseFunction::extraGal(double energy)const
{
    static double E0(100), flux(1.5E-5), alpha(1.1);
    return flux*alpha*pow(E0/energy, alpha)/energy;
}

int DiffuseFunction::layer(double e)const
{
    static double log2(log(2.0));
    double step ( log(e/s_emin)/log2 );
    return static_cast<int>(step);
}

void DiffuseFunction::setEnergy(double e)const
{
    static double log2(log(2.0));
    m_energy = e;
    double step ( log(e/s_emin)/log2 );
    m_layer = static_cast<int>(step);
    if( m_layer >= layers()-1 ){
        m_layer = layers()-1; // set for maximum
        m_fract=0;
    }else {
        m_fract = step - m_layer;
        m_fract = sqrt(m_fract); // interpolate according to e**-2?.
    }
}

double DiffuseFunction::h(double r, double alpha)
{
    return (1.-pow(r,-alpha))/alpha;
}

double DiffuseFunction::operator()(const astro::SkyDir& dir)const {
        double a ( m_data.pixelValue(dir, m_layer) );
        if( m_fract==0) return a;
        double b(m_data.pixelValue(dir, m_layer+1) );
        return b*m_fract + a*(1-m_fract) 
            + extraGal(m_energy);
    }

double DiffuseFunction::integral(const astro::SkyDir& dir, double a, double b)const
{
    static double log2(log(2.));
    ///@todo: generalize this for intervals larger than a factor of 3 
    int k(layer(a));   // nearest index
    double Ek( energy_bin(k) ); // energy for nearest value

    // flux*energy values for this and the next bin
    double Fk( Ek*m_data.pixelValue(dir,k) ) 
        ,  Fkp( 2.*Ek*m_data.pixelValue(dir,k+1) );

    // the power law index, and final integral
    double alpha( log(Fk/Fkp)/log2 )
        ,  Q( Fk * (h(b/Ek,alpha) - h(a/Ek,alpha)) );
    return Q;
}

std::vector<double> DiffuseFunction::integral(const astro::SkyDir& dir, const std::vector<double>&energies)const
{
    std::vector<double> result;
    static double infinity (300000);

    for( std::vector<double>::const_iterator it = energies.begin(); it!=energies.end(); ++it){
        std::vector<double>::const_iterator next(it+1);
        double a = *it;
        double b = next!=energies.end()? *next : infinity;
        result.push_back(integral(dir,a,b));
    }
    return result;
}
