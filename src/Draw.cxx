/** @file Draw.cxx
@brief implementation of Draw

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Draw.cxx,v 1.7 2008/05/28 21:39:32 burnett Exp $

*/


#include "pointlike/Draw.h"

#include "astro/SkyDir.h"

#include "skymaps/SkyImage.h"
#include "skymaps/BinnedPhotonData.h"

using namespace pointlike;
using astro::SkyDir;
using skymaps::BinnedPhotonData;
using skymaps::SkyImage;

Draw::Draw(const BinnedPhotonData& map)
: m_map(map)
, m_galactic(true)
, m_proj("")
{}

void Draw::region(const astro::SkyDir& dir, std::string outputFile, double pixel, 
                  double fov, bool smooth, int mincount)
{

    std::string proj (fov>90? "AIT":"ZEA");
    if( !m_proj.empty()){ proj = m_proj;}
                
    int layers(1);
    SkyImage image(dir, outputFile, pixel, fov, layers, proj,  m_galactic);

    /// @class SkyDensity
    /// @brief adapt a BinnedPhotonData to give density
    class SkyDensity : public astro::SkyFunction
    {
        public:
            SkyDensity(const BinnedPhotonData& data, bool smooth, int mincount):
              m_data(data),
              m_smooth(smooth),
              m_mincount(mincount) {}

              double operator()(const astro::SkyDir & sd) const 
              {
                  double  value;
                  if (m_smooth)
                    value = m_data.smoothDensity(sd, m_mincount);
                  else
                    value = m_data.density(sd);
                  return value;    
              }
        private:
            const BinnedPhotonData& m_data;
            bool m_smooth;
            int m_mincount;
    };

    image.fill(SkyDensity(m_map, smooth, mincount), 0); // PhotonMap is a SkyFunction of the density 
    std::cout 
        <<   "\t minimum "<< image.minimum()
        << "\n\t maximum "<< image.maximum()
        << "\n\t average "<< (image.total()/image.count())
        << std::endl;

#if 0
    /// @class SkyCount
    /// @brief adapt a PhotonMap to give weighted value
    class SkyCount : public astro::SkyFunction {
    public:
        SkyCount(const PhotonMap& data, int level, CountType counts=WEIGHTED):
          m_data(data), m_level(level), m_counts(counts) {}

          double operator()(const astro::SkyDir & sd) const {
              bool includeChildren = (m_counts == CHILDREN || m_counts == WEIGHTED),
                  weighted        = (m_counts == WEIGHTED);
              double  value = m_data.photonCount(healpix::HealPixel(sd, m_level), includeChildren, weighted); 
              return value;    
          }
    private:
        const PhotonMap& m_data;
        int m_level;
        CountType m_counts;
    };

    // Where HealPixel width for level > display pixel width, use fill().
    int layer = 1, level = m_map.minLevel(), minLevel=level;
    for (; level < minLevel + m_map.levels()
        && (sqrt(healpix::HealPixel(SkyDir(0,0), level).area())) * (180/M_PI) >= 1.15 * pixel;
        ++level, ++ layer)
    {
        std::cout << "Filling image layer "<<layer<< " with  counts on level "<<level << std::endl;
        image.fill(SkyCount(m_map, level, m_countType), layer);  
    }
    // Where HealPixel width for level <= display pixel width, use addPoint().

    std::cout << "Filling layers "<< layer << " and above with ... ";
    int points(0);
    for (PhotonMap::const_iterator it = m_map.begin();
        it!= m_map.end(); ++it)
    {
        if (it->first.level() >= level) {
            int layer = it->first.level() - m_map.minLevel() + 1;
            if(image.addPoint((it->first)(), it->second, layer)) points+= it->second;
        }
    }
    std::cout <<  points << " hit display pixels" << std::endl;
#endif

    std::cout << "Writing image to file \""
        << outputFile << "\""<<std::endl;

}
void Draw::sky(std::string outputfile, double pixel)
{
    region(SkyDir(0,0, SkyDir::GALACTIC), outputfile, pixel, 180.);
}
