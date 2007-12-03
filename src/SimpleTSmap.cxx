/** @file SimpleTSmap.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SimpleTSmap.cxx,v 1.1 2007/11/27 04:35:33 burnett Exp $
*/

#include "pointlike/SimpleTSmap.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/PhotonMap.h"

#include "healpix/HealPixel.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <time.h>


namespace{  // anonymous namespace 

    void timer(std::string label="",std::ostream& out = std::cout)
    {
        time_t t1;
        time(&t1);
        out << "\n"<< ctime(&t1);
        if( !label.empty()) out << label<<std::endl;;
    }
    std::string names[]={"TS", "signal", "gradient"};

}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace pointlike;
using astro::SkyDir;
using healpix::HealPixel;

SimpleTSmap::SimpleTSmap(const pointlike::PhotonMap& pmap, 
                         const pointlike::SkySpectrum& background)
: m_pmap(pmap)
, m_background(background)
, m_level(8)
{
}

void SimpleTSmap::run( astro::SkyDir& center, double radius, int level, bool verbose)
{
    timer();
    m_level= level;

    std::vector<std::pair<healpix::HealPixel, int> > v;
    PointSourceLikelihood::set_levels(m_level, m_level);

    int num(0);
    // Extract the pixels to be examined
    m_pmap.extract_level(center, radius, v, m_level, true);
    std::cout << "Examining " << v.size() << " pixels at level " << m_level << std::endl;
    for(std::vector<std::pair<healpix::HealPixel, int> >::const_iterator it = v.begin();
        it != v.end(); ++it, ++num)
    {

        HealPixel hpix(it->first);
        size_t index = hpix.index();
        SkyDir sd;
        int count = static_cast<int>(m_pmap.photonCount(hpix, sd));

        PointSourceLikelihood ps(m_pmap, "test", sd);
        //ps.set_verbose(true);
        // todo: fix this need for a cast
        ps.set_diffuse(const_cast<pointlike::SkySpectrum*>(&m_background));
        double ts(ps.maximize());
        //ps.printSpectrum();
        if( ts>5 ){
            SimpleLikelihood& like(*ps[m_level]);
            m_tsmap[index].push_back(ts);
            m_tsmap[index].push_back(like.signal());
            m_tsmap[index].push_back(like.gradient().mag());
            if( verbose ){
                std::cout << hpix.index() << "\t" << ts << std::endl;
            }
        }

    }
    timer();
    std::cout << "Found solutions at " << m_tsmap.size() << " pixels." << std::endl;
}

SimpleTSmap::~SimpleTSmap()
{}

double SimpleTSmap::value(const astro::SkyDir& dir, double )const
{
    HealPixel hp(dir, m_level);
    std::map<int, std::vector<float> >::const_iterator it (m_tsmap.find(hp.index()) );
    return it != m_tsmap.end() ? it->second[0] : 0;    
}

double SimpleTSmap::integral(const astro::SkyDir& dir, double , double )const
{
    return value(dir, 1000);
}

std::string SimpleTSmap::name()const
{
    std::stringstream msg;
    msg <<"TS map, at level" << m_level;
    return msg.str();
}
void SimpleTSmap::save(std::string filename)
{
    std::ofstream out(filename.c_str());
    out << m_level << std::endl;
    for( std::map<int,std::vector<float> >::const_iterator it = m_tsmap.begin(); it!=m_tsmap.end(); ++it){
        out << it->first << "\t";
        std::copy(it->second.begin(), it->second.end() , 
                         std::ostream_iterator<float>(out, "\t") );
        out << std::endl;
    }
}

void SimpleTSmap::restore(std::string filename)
{
    std::ifstream in(filename.c_str());
    in >> m_level;
    while(!in.eof()){
        int index;
        float ts;
        in >> index >> ts;
        m_tsmap[index].push_back(ts);
    }
}

float SimpleTSmap::operator[] (int index)const
{
    std::map<int,std::vector<float> >::const_iterator it = m_tsmap.find(index);
    return it==m_tsmap.end()? 0 : it->second[0];
}


    
