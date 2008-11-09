/** @file BandBackground.h
    @brief declare BandBackground class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/BandBackground.h,v 1.2 2008/10/21 02:50:45 burnett Exp $
*/

#ifndef pointlike_BandBackground_h
#define pointlike_BandBackground_h

#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"

#include "skymaps/CompositeSkySpectrum.h"
#include "skymaps/Background.h"
#include "skymaps/Band.h"

namespace pointlike {


/** @class BandBackground
    @brief wrap a background, return rates for event class and energy range determined by a Band

    Also implements a function to integrate over a ROI.
*/

class BandBackground : public astro::SkyFunction {

public:

    /** @brief ctor
        @param background a composite; the first element should be a Background object
        @param band the band to use to select energy range and event class

    */
    BandBackground(const skymaps::SkySpectrum& background, const skymaps::Band& band);

    virtual ~BandBackground(){};

    double operator()(const astro::SkyDir& dir)const;

    double average(const astro::SkyDir& dir, double angle, double tolerance)const;
    static void set_verbose(bool q=true);
            
private:
    // Calculate average for a given level
    double level_ave(const astro::SkyDir& dir, double angle, int level) const;
#if 0
    const skymaps::CompositeSkySpectrum & m_background;
    const skymaps::Background* m_diffuse; ///< extract pointer to the diffuse cmponent
#else
    const skymaps::SkySpectrum& m_background;
#endif
    const skymaps::Band& m_band;
    double m_emin, m_emax;
    int m_event_class;

};

}

#endif
