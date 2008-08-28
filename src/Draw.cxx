/** @file Draw.cxx
@brief implementation of Draw

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Draw.cxx,v 1.10 2008/08/21 03:23:31 burnett Exp $

*/


#include "pointlike/Draw.h"

#include "astro/SkyDir.h"

#include "skymaps/SkyImage.h"
#include "skymaps/BinnedPhotonData.h"
#include "skymaps/SkySpectrum.h"

using namespace pointlike;
using astro::SkyDir;
using skymaps::BinnedPhotonData;
using skymaps::SkyImage;

Draw::Draw(const BinnedPhotonData& map)
: m_map(map)
, m_galactic(true)
, m_proj("")
, m_exposure(0) // default: do not apply
, m_layers(1)
{}

void Draw::region(const astro::SkyDir& dir, std::string outputFile, double pixel, 
                  double fov, bool smooth, int mincount)
{

    std::string proj (fov>90? "AIT":"ZEA");
    if( !m_proj.empty()){ proj = m_proj;}
                
    SkyImage image(dir, outputFile, pixel, fov, m_layers, proj,  m_galactic);

    /// @class SkyDensity
    /// @brief adapt a BinnedPhotonData to give density
    class SkyDensity : public astro::SkyFunction
    {
        public:
            SkyDensity(const BinnedPhotonData& data, bool smooth, int mincount
                ,const skymaps::SkySpectrum* exposure):
              m_data(data),
              m_smooth(smooth),
              m_mincount(mincount)
              ,m_exposure(exposure)
              {}

              double operator()(const astro::SkyDir & sd) const 
              {
                  double  value;
                  if (m_smooth)
                    value = m_data.smoothDensity(sd, m_mincount);
                  else
                    value = m_data.density(sd);

                  if(m_exposure!=0){
                      // note we are not using energy dependence here
                      double exposure( (*m_exposure)(sd) );
                      if( exposure>0.) value /= exposure;
                  }
                  return value;    
              }
        private:
            const BinnedPhotonData& m_data;
            bool m_smooth;
            int m_mincount;
            const skymaps::SkySpectrum* m_exposure;
    };

    image.fill(SkyDensity(m_map, smooth, mincount, m_exposure), 0); // PhotonMap is a SkyFunction of the density 
    std::cout 
        <<   "\t minimum "<< image.minimum()
        << "\n\t maximum "<< image.maximum()
        << "\n\t average "<< (image.total()/image.count())
        << std::endl;


    /// @class SkyCount
    /// @brief adapt a BinnedPhotonData to give counts for bins with given energy
    class SkyCount : public astro::SkyFunction {
    public:
        SkyCount(const BinnedPhotonData& data, double energy): 
          m_data(data), m_energy(energy){}

          double operator()(const astro::SkyDir & sd) const {
              double  value = m_data.value(sd, m_energy); 
              return value;    
          }
    private:
        const BinnedPhotonData& m_data;
        double m_energy;
    };
#if 0

    std::cout << "Filling layers "<< layer << " and above with ... ";
    int points(0);
    for (BinnedPhotonData::const_iterator it = m_map.begin();
        it!= m_map.end(); ++it)
    {
        if (it->first.level() >= level) {
            int layer = it->first.level() - m_map.minLevel() + 1;
            if(image.addPoint((it->first)(), it->second, layer)) points+= it->second;
        }
    }
    std::cout <<  points << " hit display pixels" << std::endl;

    for(int layer=1; layer!=m_layers; ++m_layers){
       image.fill(SkyCount(m_map), 100.);
    }
#endif

    std::cout << "Writing image to file \""
        << outputFile << "\""<<std::endl;

}
void Draw::sky(std::string outputfile, double pixel, bool smooth, int mincount)
{
    region(SkyDir(0,0, SkyDir::GALACTIC), outputfile, pixel, 180., smooth, mincount);
}

void Draw::googleSky(std::string outfile, double pixelsize, bool smooth, int mincount)
{
    std::string proj = m_proj;    m_proj = "CAR";
    bool gal = m_galactic; m_galactic = false;
    region(SkyDir(180,0), outfile, pixelsize, 180);
    m_proj = proj;
    m_galactic = galactic;
}
