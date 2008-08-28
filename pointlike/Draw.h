/** @file Draw.h 
@brief declaration of the Draw wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Draw.h,v 1.7 2008/08/21 03:22:01 burnett Exp $
*/


#ifndef pointlike_Draw_h
#define pointlike_Draw_h

namespace astro { class SkyDir; }
#include <string>
#include <vector>
#include "embed_python/Module.h"
namespace skymaps{ class BinnedPhotonData; class SkySpectrum;}

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

        //! @brief all sky image, default AIT, galactic
        void sky(std::string outputfile, double pixelsize=0.1, bool smooth = false, int mincount = 0);

        //! @brief make a CAR projection for use by google sky
        void googleSky(std::string outfile, double pixelsize=0.1,bool smooth = false, int mincount = 0);

        void galactic(){m_galactic = true;}      ///< set galactic
        void equatorial(){m_galactic=false;}     ///< set equatorial
        void projection(std::string p){m_proj = p;} ///< set the projection

        void use_exposure(const skymaps::SkySpectrum* exp){m_exposure=exp;}

        ///! will add layers with count map values
        void set_layers(int n){m_layers=n;}

    private:
        const skymaps::BinnedPhotonData& m_map;
        bool m_galactic;    ///< galactic or equatorial
        std::string m_proj; ///< projection (CAR, AIT, etc.)
        const skymaps::SkySpectrum * m_exposure; ///< exposure to use for normalization, if present. (energy?)
        bool m_layers;
    };


}// namespace pointline


#endif

