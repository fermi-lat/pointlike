/** @file Draw.cxx
    @brief implementation of Draw

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Draw.cxx,v 1.6 2008/01/27 02:31:33 burnett Exp $

*/


#include "pointlike/DrawTS.h"
#include "pointlike/SourceLikelihood.h"

#include "astro/SkyDir.h"

#include "skymaps/SkyImage.h"
#include "skymaps/PhotonMap.h"
#include "skymaps/BinnedPhotonData.h"

using namespace pointlike;
using astro::SkyDir;
using skymaps::BinnedPhotonData;
using skymaps::SkyImage;

DrawTS::DrawTS(const BinnedPhotonData& map,
	       const skymaps::SkySpectrum& background): 
  m_map(map),
  m_background(background),
  m_galactic(true),
  m_proj(""){
}

void DrawTS::region(const astro::SkyDir& dir, std::string outputFile, 
		    double pixel, double fov,
		    bool fillts){
  
  //     int layers((m_countType != NONE)? 9:1);
  //     if (fillts) layers =  ((m_countType != NONE)? 10:1);
  
  std::string proj (fov>90? "AIT":"ZEA");
  if( !m_proj.empty()){ proj = m_proj;}
  
  int layers(1);
  SkyImage image(dir, outputFile, pixel, fov, layers, proj,  m_galactic);

  /// @class SkyDensity
  /// @brief adapt a BinnedPhotonData to give density
  class SkyDensity : public astro::SkyFunction {
  public:
    SkyDensity(const BinnedPhotonData& data):
      m_data(data){}

    double operator()(const astro::SkyDir & sd) const {
      double  value = m_data.density(sd);
      return value;    
    }
  private:
    const BinnedPhotonData& m_data;
  };

  image.fill(SkyDensity(m_map), 0); // PhotonMap is a SkyFunction of the density 
  std::cout 
    <<   "\t minimum "<< image.minimum()
    << "\n\t maximum "<< image.maximum()
    << "\n\t average "<< (image.total()/image.count())
    << std::endl;

  class SkyTS: public astro::SkyFunction {
  public:
    SkyTS(const skymaps::PhotonMap& data, int level, 
	  const skymaps::SkySpectrum& background):
      m_data(data), m_level(level), m_background(background){}

    double operator()(const astro::SkyDir& sd) const {
      std::cout << "Fix the () operator " << std::endl;
//       pointlike::SourceLikelihood ps(m_data, "test", sd);
//       if (m_level > 0) ps.set_levels(m_level, m_level);
//       else ps.set_levels(m_data.minLevel(), 
// 			 m_data.minLevel() + m_data.levels()-1);
//       ps.set_verbose(true);
//       ps.set_diffuse(const_cast<skymaps::SkySpectrum*>(&m_background));
//       double ts = ps.maximize();
      // 	std::cout << "-------- Returning : " << ts << " for position: " << 
      // 	  sd.ra() << " " << sd.dec() << std::endl;
//       return ts;
      return -1;
    }
  private:
    const skymaps::PhotonMap& m_data;
    int m_level;
    const skymaps::SkySpectrum& m_background;
  };

  
//   int layer = 1, level = m_map.minLevel(), minLevel=level;
//   for (; level < minLevel + m_map.levels()
// // 	 && (sqrt(healpix::HealPixel(SkyDir(0,0), level).area())) * (180/M_PI) >= 1.15 * pixel;
//        ++level, ++ layer){
//     std::cout << "Filling image layer "<<layer<< " with TS on level "
// 	      <<level << std::endl;
//     image.fill(SkyTS(m_map, level, m_background), layer);  
//   }
  
  
  std::cout << "Writing image to file \""
	    << outputFile << "\""<<std::endl;
  
}

void DrawTS::sky(std::string outputfile, double pixel){

  region(SkyDir(0,0, SkyDir::GALACTIC), outputfile, pixel, 180., true);
}
