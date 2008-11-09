/** @file BandBackground.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/BandBackground.cxx,v 1.3 2008/10/21 02:50:45 burnett Exp $

*/

#include "pointlike/BandBackground.h"

#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"

using namespace pointlike;
using astro::SkyDir;
namespace {
    bool verbose(false);
}
void BandBackground::set_verbose(bool q){verbose=q;} 

BandBackground::BandBackground(const skymaps::SkySpectrum& background, const skymaps::Band& band)
: m_background(background)
, m_band(band)
, m_emin(band.emin())
, m_emax(band.emax())
, m_event_class(band.event_class())
{}



double BandBackground::operator()(const astro::SkyDir& dir)const{
    return m_background.band_value(dir, m_band);
}


double BandBackground::average(const astro::SkyDir& dir, double angle, double tolerance)const
{
    using astro::SkyDir;
    using healpix::Healpix;

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
double BandBackground::level_ave(const astro::SkyDir& dir, double angle, int level) const
{   

    int nside(1 << level);
    std::vector<int> v;
    healpix::Healpix hpx(nside, healpix::Healpix::NESTED, astro::SkyDir::GALACTIC);
    hpx.query_disc(dir, angle, v); 
    double av(0);

    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {
        healpix::HealPixel hp(*it, level);
        av += (*this) (hp());
    }

    return av/v.size();
}
