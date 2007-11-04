/** @file SkySpectrum.cxx
    @brief implement class SkySpectrum

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SkySpectrum.cxx,v 1.1 2007/11/01 21:33:33 burnett Exp $

*/

#include "pointlike/SkySpectrum.h"
#include "astro/SkyDir.h"
#include "astro/Healpix.h"
#include "astro/HealPixel.h"

#include <vector>
#include <map>
#include <cassert>

using namespace pointlike;
SkySpectrum::SkySpectrum(double energy)
{
    setEnergy(energy);
}
void SkySpectrum::setEnergy(double e)const
{
    m_energy =e;
    m_use_range = false;
}

void SkySpectrum::setEnergyRange(double emin, double emax)const
{
    m_emin = emin; m_emax=emax;
    assert( emin>0 && emax/emin<3); // otherwise integral not (now) valid
    m_use_range=true;
}

double SkySpectrum::operator()(const astro::SkyDir& dir)const
{
    if( m_use_range){
        // evaluate the integral for the current energy range
        return integral(dir, m_emin, m_emax);
    }
    return value(dir, m_energy); //?? + extraGal(m_energy);

}

double SkySpectrum::operator()(const astro::SkyDir& dir, double energy)const
{
    setEnergy(energy); return (*this)(dir);
}

double SkySpectrum::operator()(const astro::SkyDir& dir, double emin, double emax)const
{
    setEnergyRange(emin, emax); return (*this)(dir);
}


double SkySpectrum::average(const astro::SkyDir& dir, double angle, double tolerance)const
{
    using astro::SkyDir;
    using astro::Healpix;

    static std::map<double, int> width_map;
    static bool map_built(false);

    int level, min_level = 6, max_level = 13;
    double result(0.0);

    // Get value for one point at center
    double previous = (*this) (dir);
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
    for(level += 1 ; fabs(result/previous -1.) > tolerance && level < max_level; ++ level)
    {
        previous = result;
        result = level_ave(dir, angle, level);
    }

    return result;

}

// Calculate average for a given level
double SkySpectrum::level_ave(const astro::SkyDir& dir, double angle, int level) const
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
    