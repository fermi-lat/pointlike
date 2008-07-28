/** @file Draw.h 
@brief declaration of the Draw wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Draw.h,v 1.5 2008/07/19 14:46:09 burnett Exp $
*/


#ifndef pointlike_Draw_h
#define pointlike_Draw_h

namespace astro { class SkyDir; }
#include <string>
#include <vector>
#include "embed_python/Module.h"
namespace skymaps{ class BinnedPhotonData; }

namespace pointlike {

    //! @class Draw
    //! @brief manage creating images to FITS files from a BinnedPhotonData object
    class Draw {
    public:


        Draw(const skymaps::BinnedPhotonData& map);

        //! create FITS image file using the data
        //! @param dir center
        //! @param outputFile file to write
        //! @param pixelsize in degrees
        //! @param fov  field of view (deg) if 180, use car

        void region(const astro::SkyDir& dir, std::string outputFile, double pixelsize,
                    double fov, bool smooth = false, int mincount = 0);

        void sky(std::string outputfile, double pixelsize, bool smooth = false, int mincount = 0);

        void galactic(){m_galactic = true;}      ///< set galactic
        void equatorial(){m_galactic=false;}     ///< set equatorial
        void projection(std::string p){m_proj = p;} ///< set the projection

    private:
        const skymaps::BinnedPhotonData& m_map;
        bool m_galactic;    ///< galactic or equatorial
        std::string m_proj; ///< projection (CAR, AIT, etc.)
    };


}// namespace pointline


#endif

