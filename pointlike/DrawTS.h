/** @file DrawTS.h 
@brief declaration of the DrawTS wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/DrawTS.h,v 1.4 2008/05/28 21:39:32 burnett Exp $
*/


#ifndef pointlike_DrawTS_h
#define pointlike_DrawTS_h

namespace astro { class SkyDir; }
#include <string>
#include <vector>
#include "embed_python/Module.h"
#include "skymaps/SkySpectrum.h"
namespace skymaps{ class BinnedPhotonData; }

namespace pointlike {
  
  //! @class DrawTS
  //! @brief manage creating images to FITS files from a BinnedPhotonData object
  class DrawTS {
  public:
    
    
    DrawTS(const skymaps::BinnedPhotonData& map, 	   
	   const skymaps::SkySpectrum& background);
  
    //! create FITS image file using the data
  //! @param dir center
  //! @param outputFile file to write
  //! @param pixelsize in degrees
  //! @param fov  field of view (deg) if 180, use car
    
    void region(const astro::SkyDir& dir, std::string outputFile, 
		double pixelsize, double fov, bool fillts);
    
    void sky(std::string outputfile, double pixelsize);
    
    void galactic(){m_galactic = true;}      ///< set galactic
    void equatorial(){m_galactic=false;}     ///< set equatorial
    void projection(std::string p){m_proj = p;} ///< set the projection
    
  private:
    const skymaps::BinnedPhotonData& m_map;
    const skymaps::SkySpectrum& m_background;
    bool m_galactic;    ///< galactic or equatorial
    std::string m_proj; ///< projection (CAR, AIT, etc.)
  };
  

}// namespace pointline


#endif

