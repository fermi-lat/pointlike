/** @file Draw.h 
@brief declaration of the Draw wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Draw.h,v 1.9 2008/09/05 23:40:10 funk Exp $
*/


#ifndef pointlike_Draw_h
#define pointlike_Draw_h

namespace astro { class SkyDir; }
#include <string>
#include <vector>
#include "embed_python/Module.h"
namespace skymaps{ class BinnedPhotonData; class SkySpectrum;}

namespace pointlike {

    class Data;  // forward declaration

    //! @class Draw
    //! @brief manage creating images to FITS files from a BinnedPhotonData object
    class Draw {
    public:

        //! @brief ctor sets data
      
        Draw(const skymaps::BinnedPhotonData& map, 
            const skymaps::SkySpectrum* background = 0,
            bool ts = false);

        Draw(const Data& data);

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
        int set_layers(int n){int t = m_layers; m_layers=n; return t;}

    private:
        const skymaps::BinnedPhotonData& m_map;
        const skymaps::SkySpectrum* m_background;
        bool m_galactic;    ///< galactic or equatorial
        std::string m_proj; ///< projection (CAR, AIT, etc.)
        const skymaps::SkySpectrum * m_exposure; ///< exposure to use for normalization, if present. (energy?)
        int m_layers;
        bool m_ts;
    };


}// namespace pointline


#endif

