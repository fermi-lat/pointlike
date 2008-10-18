/** @file BandBackground.h
    @brief declare BandBackground class

$Header$
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
    BandBackground(const skymaps::CompositeSkySpectrum& background, const skymaps::Band& band);

    virtual ~BandBackground(){};

    double operator()(const astro::SkyDir& dir)const;

    double average(const astro::SkyDir& dir, double angle, double tolerance)const;

            
private:
    // Calculate average for a given level
    double level_ave(const astro::SkyDir& dir, double angle, int level) const;

    const skymaps::CompositeSkySpectrum & m_background;
    const skymaps::Background* m_diffuse; ///< extract pointer to the diffuse cmponent
    
    double m_emin, m_emax;
    int m_event_class;
};

}

#endif
