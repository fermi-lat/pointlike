/** @file Draw.h 
@brief declaration of the Draw wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Draw.h,v 1.16 2009/03/06 21:41:47 markusa Exp $
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
    class SourceLikelihood;
    class PointSourceLikelihood;

    //! @class Draw
    //! @brief manage creating images to FITS files from a BinnedPhotonData object
    class Draw {
    public:

        //! @brief ctor sets data
      
      Draw(skymaps::BinnedPhotonData& map, 
	   const skymaps::SkySpectrum* background = 0,
	   bool ts = false, double emin = 100, double minalpha = 0.05);
      
      Draw(Data& data);

        //! create FITS image file using the data
        //! @param dir center
        //! @param outputFile file to write
        //! @param pixelsize in degrees
        //! @param fov  field of view (deg) if 180, use car

        void region(const astro::SkyDir& dir, std::string outputFile, double pixelsize,
                     double fov, bool smooth = false, int mincount = 0, int kernel = 0, double smooth_radius = 3);
        
	void density(const astro::SkyDir& dir, std::string outputFile, double pixelsize,
                    double fov, bool smooth = false, int mincount = 0, int kernel = 0, double smooth_radius = 3);

        //! @brief all sky image, default AIT, galactic
        void sky(std::string outputfile, double pixelsize=0.1, bool smooth = false, int mincount = 0);

        //! @brief make a CAR projection for use by google sky
        void googleSky(std::string outfile, double pixelsize=0.1,bool smooth = false, int mincount = 0);

        void galactic(){m_galactic = true;}      ///< set galactic
        void equatorial(){m_galactic=false;}     ///< set equatorial
        void projection(std::string p){m_proj = p;} ///< set the projection
        void zenith(){m_zenith = true;}

        void use_exposure(const astro::SkyFunction* exp){m_exposure=exp;}
        void use_exposure2(const astro::SkyFunction* exp){m_exposure2=exp;}
        void use_exposure(const skymaps::SkySpectrum* exp){m_exposure=(const astro::SkyFunction*)(exp);}
        void use_exposure2(const skymaps::SkySpectrum* exp){m_exposure2=(const astro::SkyFunction*)(exp);}

        ///! will add layers with count map values
        int set_layers(int n){int t = m_layers; m_layers=n; return t;}

    private:
      skymaps::BinnedPhotonData& m_map;
      const skymaps::SkySpectrum* m_background;
      bool m_galactic;    ///< galactic or equatorial
      bool m_zenith;      ///< set for zenith coords
      std::string m_proj; ///< projection (CAR, AIT, etc.)
      const astro::SkyFunction* m_exposure; ///< exposure to use for normalization, if present. (energy?)
      const astro::SkyFunction* m_exposure2; // possible entry for back exposure
      int m_layers;
      double m_emin;
      double m_minalpha;
      bool m_ts;
    };

    /** @class SkyDensity
        @brief adapt a BinnedPhotonData to give density
     
      This should go with skymaps,  
    */
    class SkyDensity : public astro::SkyFunction
    {
    public:
        //SkyDensity(const skymaps::BinnedPhotonData& data, bool smooth=false, int mincount=0
        //    ,const astro::SkyFunction* exposure=0);
        SkyDensity(const skymaps::BinnedPhotonData& data, bool smooth=false, int mincount=0
            ,const astro::SkyFunction* exposure=0,const astro::SkyFunction* exposure2=0,int kernel=0,double smooth_radius=3);

        //! @brief ctor to display the given band only
        SkyDensity(const skymaps::Band& band);



        double operator()(const astro::SkyDir & sd) const ;

    private:
        const skymaps::BinnedPhotonData* m_data;
        const skymaps::Band* m_band;
        bool m_smooth;
        int m_mincount, m_kernel;
        double m_radius;
        const astro::SkyFunction* m_exposure;
        const astro::SkyFunction* m_exposure2;
    };

}// namespace pointline


#endif

