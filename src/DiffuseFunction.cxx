/** @file DiffuseFunction.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/DiffuseFunction.cxx,v 1.4 2007/09/28 20:17:10 burnett Exp $
*/

#include "pointlike/DiffuseFunction.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "astro/Healpix.h"
#include "astro/HealPixel.h"
#include <cmath>
#include <map>

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

double DiffuseFunction::average(const astro::SkyDir& dir, double angle, double tolerance)const
{
    using astro::SkyDir;
    using astro::Healpix;

    static std::map<double, int> width_map;
    static bool map_built(false);

    int level, min_level = 6, max_level = 13;
    double result(0.0), previous(1e20);

    // Get value for one point at center
    previous = (*this) (dir);
    if (tolerance >= 0.5)  // If tolerance is higher than this, just return value at center.
        return previous;

    /* Build map of pixel widths by level, if not done yet.  Store twice the healpixel
       width for easy comparison. */
    if (!map_built)
    {
        width_map.clear();
        for (level = min_level; level <= max_level; ++level)
        {
            int nside(1 << level);
            int npix(12 * nside * nside);
            double width = sqrt(4 * M_PI / npix);  // Width of healpixel in radians
            width_map[2 * width] = level;
        }
        map_built = true;
    }

    // Use map to determine starting pixel level
    std::map<double, int>::iterator it = width_map.lower_bound(angle);
    if (it == width_map.end() || (it->first > angle && it != width_map.begin()))
        --it;
    level = it->second;

    // Get value for starting pixel level
    result = level_ave(dir, angle, level);

    // Iterate until result changes less than tolerance
    for(level += 1 ; fabs(result - previous) > tolerance && level < max_level; ++ level)
    {
        previous = result;
        result = level_ave(dir, angle, level);
    }

    return result;

}

// Calculate average for a given level
double DiffuseFunction::level_ave(const astro::SkyDir& dir, double angle, int level) const
{   
    using astro::Healpix;

    int nside(1 << level);
    std::vector<int> v;
    Healpix hpx(nside, astro::Healpix::NESTED, astro::SkyDir::GALACTIC);
    hpx.query_disc(dir, angle, v); 
    double av(0);

    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        astro::HealPixel hp(*it, level);
        av += (*this) (hp());
    }

    return av/v.size();
}
    
