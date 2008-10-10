/** @file Draw.h 
@brief declaration of the Draw wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Draw.h,v 1.11 2008/09/24 18:01:42 burnett Exp $
*/


#ifndef pointlike_Draw_h
#define pointlike_Draw_h
#include "astro/SkyFunction.h"

namespace astro { class SkyDir; }
#include <string>
#include <vector>
#include "embed_python/Module.h"
namespace skymaps{ class BinnedPhotonData; class SkySpectrum; class Band;}

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
        void zenith(){m_zenith = true;}

        void use_exposure(const astro::SkyFunction* exp){m_exposure=exp;}

        ///! will add layers with count map values
        int set_layers(int n){int t = m_layers; m_layers=n; return t;}

    private:
        const skymaps::BinnedPhotonData& m_map;
        const skymaps::SkySpectrum* m_background;
        bool m_galactic;    ///< galactic or equatorial
        bool m_zenith;      ///< set for zenith coords
        std::string m_proj; ///< projection (CAR, AIT, etc.)
        const astro::SkyFunction* m_exposure; ///< exposure to use for normalization, if present. (energy?)
        int m_layers;
        bool m_ts;
    };

    /** @class SkyDensity
        @brief adapt a BinnedPhotonData to give density
     
      This should go with skymaps,  
    */
    class SkyDensity : public astro::SkyFunction
    {
    public:
        SkyDensity(const skymaps::BinnedPhotonData& data, bool smooth=false, int mincount=0
            ,const astro::SkyFunction* exposure=0);

        //! @brief ctor to display the given band only
        SkyDensity(const skymaps::Band& band);



        double operator()(const astro::SkyDir & sd) const ;

    private:
        const skymaps::BinnedPhotonData* m_data;
        const skymaps::Band* m_band;
        bool m_smooth;
        int m_mincount;
        const astro::SkyFunction* m_exposure;
    };

}// namespace pointline


#endif

