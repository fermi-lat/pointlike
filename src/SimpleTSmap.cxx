/** @file SimpleTSmap.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SimpleTSmap.cxx,v 1.5 2008/05/02 23:31:04 burnett Exp $
*/

#include "pointlike/SimpleTSmap.h"
#include "pointlike/PointSourceLikelihood.h"
#include "skymaps/BinnedPhotonData.h"

#include "healpix/HealPixel.h"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <iterator>
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

    double tsmin(10.);

}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace pointlike;
using astro::SkyDir;
using healpix::HealPixel;

SimpleTSmap::SimpleTSmap(const skymaps::BinnedPhotonData& pmap)
: m_pmap(pmap)
{
}

void SimpleTSmap::run( astro::SkyDir& center, double radius, double emin, bool verbose)
{
    timer();

    std::vector<std::pair<int, int> > v;
    skymaps::BinnedPhotonData::const_iterator bpit(m_pmap.begin() );

    while( (*bpit).emin()< emin) { 
        ++bpit;
    }
    const skymaps::Band& b (* --bpit); 
    m_nside = b.nside();
    b.query_disk(center, radius, v);

    std::cout << "Examining " << v.size() << " pixels with nside= " << m_nside << std::endl;
    int num(0);
    for(std::vector<std::pair<int, int> >::const_iterator it = v.begin();
        it != v.end(); ++it, ++num)
    {
        int index(it->first);
        const SkyDir& sd(b.dir(index));
#ifdef OLD
        int count = static_cast<int>(m_pmap.photonCount(hpix, sd));
#else
        int count (it->second);
#endif
        PointSourceLikelihood ps(m_pmap, "test", sd);
        //ps.set_verbose(true);
        // todo: fix this need for a cast
        double ts(ps.maximize());
        //ps.printSpectrum();
        if( ts>tsmin ){

            m_tsmap[index].push_back(ts);
            m_tsmap[index].push_back(count);
#if 0 //fix this?
            SimpleLikelihood& like( *ps[m_level] );
            m_tsmap[index].push_back(like.signal());
            m_tsmap[index].push_back(like.gradient().mag());
#endif
        }

    }
    timer();
    std::cout << "Found solutions at " << m_tsmap.size() << " pixels." << std::endl;
}

SimpleTSmap::~SimpleTSmap()
{}

double SimpleTSmap::value(const astro::SkyDir& dir, double )const
{
#if 0
    HealPixel hp(dir, m_level);
    std::map<int, std::vector<float> >::const_iterator it (m_tsmap.find(hp.index()) );
    return it != m_tsmap.end() ? it->second[0] : 0;   
#else
    return 0;
#endif
}

double SimpleTSmap::integral(const astro::SkyDir& dir, double , double )const
{
    return value(dir, 1000);
}

std::string SimpleTSmap::name()const
{
    std::stringstream msg;
    msg <<"TS map, with nside=" << m_nside;
    return msg.str();
}
void SimpleTSmap::save(std::string filename)
{
    std::ofstream out(filename.c_str());
    out << m_nside << std::endl;
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
    in >> m_nside;
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


    
