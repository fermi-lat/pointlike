/** @file Draw.cxx
@brief implementation of Draw

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Draw.cxx,v 1.1 2007/08/30 18:12:54 burnett Exp $

*/


#include "pointlike/Draw.h"

#include "astro/SkyDir.h"

#include "map_tools/PhotonMap.h"
#include "map_tools/SkyImage.h"

using namespace pointlike;
using astro::SkyDir;

Draw::Draw(const map_tools::PhotonMap& map)
: m_map(map)
, m_galactic(true)
, m_proj("")
, m_countType(WEIGHTED)
{}

void Draw::region(const astro::SkyDir& dir, std::string outputFile, double pixel, double fov )
{
    int layers((m_countType != NONE)? 9:1);

    std::string proj (fov>90? "AIT":"ZEA");
    if( !m_proj.empty()){ proj = m_proj;}
                
    map_tools::SkyImage image(dir, outputFile, pixel, fov, layers, proj,  m_galactic);
    std::cout << "Filling image layer 0 with density ..." << std::endl;
    image.fill(m_map, 0); // PhotonMap is a SkyFunction of the density 
    std::cout 
        <<   "\t minimum "<< image.minimum()
        << "\n\t maximum "<< image.maximum()
        << "\n\t average "<< (image.total()/image.count())
        << std::endl;

    if( layers==1) return;

    /// @class SkyCount
    /// @brief adapt a PhotonMap to give weighted value
    class SkyCount : public astro::SkyFunction {
    public:
        SkyCount(const map_tools::PhotonMap& data, int level, CountType counts=WEIGHTED):
          m_data(data), m_level(level), m_counts(counts) {}

          double operator()(const astro::SkyDir & sd) const {
              bool includeChildren = (m_counts == CHILDREN || m_counts == WEIGHTED),
                  weighted        = (m_counts == WEIGHTED);
              double  value = m_data.photonCount(astro::HealPixel(sd, m_level), includeChildren, weighted); 
              return value;    
          }
    private:
        const map_tools::PhotonMap& m_data;
        int m_level;
        CountType m_counts;
    };

    // Where HealPixel width for level > display pixel width, use fill().
    int layer = 1, level = m_map.minLevel(), minLevel=level;
    for (; level < minLevel + m_map.levels()
        && (sqrt(astro::HealPixel(SkyDir(0,0), level).area())) * (180/M_PI) >= 1.15 * pixel;
        ++level, ++ layer)
    {
        std::cout << "Filling image layer "<<layer<< " with  counts on level "<<level << std::endl;
        image.fill(SkyCount(m_map, level), layer);  
    }
    // Where HealPixel width for level <= display pixel width, use addPoint().

    std::cout << "Filling layers "<< layer << " and above with ... ";
    int points(0);
    for (map_tools::PhotonMap::const_iterator it = m_map.begin();
        it!= m_map.end(); ++it)
    {
        if (it->first.level() >= level) {
            int layer = it->first.level() - m_map.minLevel() + 1;
            if(image.addPoint((it->first)(), it->second, layer)) points+= it->second;
        }
    }
    std::cout <<  points << " hit display pixels" << std::endl;

    std::cout << "Writing image to file \""
        << outputFile << "\""<<std::endl;

}
void Draw::sky(std::string outputfile, double pixel)
{
    region(SkyDir(0,0, SkyDir::GALACTIC), outputfile, pixel, 180.);
}
