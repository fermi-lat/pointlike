/** @file PhotonMap.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/PhotonMap.cxx,v 1.3 2007/11/20 23:14:29 burnett Exp $
*/

#include "pointlike/PhotonMap.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <cmath>
#include <utility>
#include <stdexcept>
//#include <errno.h>

using healpix::HealPixel;
using healpix::Healpix;
using astro::SkyDir;
#include "astro/Photon.h"

using astro::Photon;
using healpix::HealPixel;

using namespace pointlike;

#if 0
PhotonMap::PhotonMap(const std::string & inputFile,const std::string & tablename )
: map_tools::PhotonMap(inputFile, tablename)
, m_name(inputFile)
{
}

#endif

   ///@brief data value for bin with given energy 
    ///@param e energy in MeV
double PhotonMap::value(const SkyDir& dir, double energy)const
{
    Photon gam(dir, energy, 0); // just for the next interface :-)
    PhotonMap* self = const_cast<PhotonMap*>(this);  // pixel not const
    healpix::HealPixel pix( self->pixel(gam) );
    return photonCount(pix) ;
}

    ///@brief integral for the energy limits, in the given direction
    /// Assume that request is for an energy bin.
double PhotonMap::integral(const SkyDir& dir, double a, double b)const
{
    return value(dir, sqrt(a*b));
}


PhotonMap::PhotonMap(double emin, double eratio, int nlevels, int minlevel)
: m_emin(emin)
, m_logeratio(log(eratio))
, m_levels(nlevels)
, m_minlevel(minlevel)
, m_photons(0)
, m_pixels(0)
{}

PhotonMap::PhotonMap(const std::string & inputFile, const std::string & tablename)
: m_photons(0)
, m_pixels(0)
{
    const tip::Table & table=*tip::IFileSvc::instance().readTable(inputFile, tablename);
    const tip::Header& hdr = table.getHeader();
    double eratio;
    int stored_photons(0), stored_pixels(0);

	// Guard against headers not being found in fits file.  Set to default on error

	try	{hdr["EMIN"].get(m_emin);} catch (const std::exception& ) {m_emin = 100.;}
	try
	{
		hdr["ERATIO"].get(eratio);
		m_logeratio = log(eratio);
	}
	catch (const std::exception& ) {m_logeratio = log(2.35);}
    try	{hdr["LEVELS"].get(m_levels);} catch (const std::exception& ) {m_levels = 8;}
    try	{hdr["MINLEVEL"].get(m_minlevel);} catch (const std::exception& ) {m_minlevel = 6;}
	try
	{
		hdr["PHOTONS"].get(stored_photons);
		hdr["PIXELS"].get(stored_pixels);
	}
	catch (const std::exception& ) {}

    tip::Table::ConstIterator itor = table.begin();
    std::cout << "Creating PhotonMap from file " << inputFile << ", table " << tablename << std::endl;

    for(tip::Table::ConstIterator itor = table.begin(); itor != table.end(); ++itor)
    {
        long level, index, count;
        (*itor)["LEVEL"].get(level);
        (*itor)["INDEX"].get(index);
        (*itor)["COUNT"].get(count);
        HealPixel p(index, level);
        this->insert( value_type(p, count) );
        m_pixels ++;
        m_photons += count;
    }
    delete &table; 
    std::cout << "Photons available: " << stored_photons 
        << "  Pixels available: " << stored_pixels <<std::endl;
    std::cout << "Photons loaded: " << m_photons 
        << "  Pixels created: " << m_pixels <<std::endl;
}



void PhotonMap::addPhoton(const Photon& gamma)
{
    if( gamma.energy() < m_emin) return;
    HealPixel p = pixel(gamma);
    iterator it = this->find(p);
    if( it==this->end()){
        // create a new pixel with the entry
        this->insert( value_type(p,1) );
        ++m_pixels;
    }else{
        // just increment the count
        ++(it->second); 
    }
    ++m_photons;
}

void PhotonMap::addPixel(const healpix::HealPixel & px, int count)
{
    this->insert( value_type(px, count) );
    m_pixels ++;
    m_photons += count;
}


HealPixel PhotonMap::pixel(const Photon& gamma)
{
    int i( static_cast<int>(log(gamma.energy()/m_emin)/m_logeratio) );
    if( i>m_levels-1) i= m_levels-1;
    return HealPixel(gamma.dir(), i+m_minlevel);
}

// Return density for a given direction, in photons/area of base pixel.
double PhotonMap::density (const SkyDir & sd) const
{
    double result(0), weight(1.0);

    for (int ebin = 0; ebin < m_levels; ++ebin, weight*=4.0) { 
        HealPixel hpx(sd, ebin+m_minlevel);
        const_iterator mit = this->find(hpx);

        if (mit != this->end()) {
            result += weight* mit->second;;
        }
    }

    return result;
}

int PhotonMap::extract(const SkyDir& dir, double radius,
                       std::vector<std::pair<HealPixel, int> >& vec,
                       int summary_level, int select_level) const
{
    //unused bool allsky(radius>=180); // maybe use to simplify below, but seems fast
    radius *= (M_PI / 180); // convert to radians
    if (summary_level == -1)
        summary_level = m_minlevel; // default level to test
    int nside( 1<< summary_level);
    int npix( 12 * nside * nside);
    int total(0);
    vec.clear();

    // Get pixels in summary level that are within radius
    std::vector<int> v;
    healpix::Healpix hpx(nside, Healpix::NESTED, HealPixel::s_coordsys);
    hpx.query_disc(dir, radius, v);  
    int max_level = m_minlevel + m_levels - 1;

    // Add summary level pixels and all their children to return vector
    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {

        HealPixel hp(*it, summary_level);
        HealPixel boundary(hp.lastChildIndex(max_level), max_level);
        for(const_iterator it2 = lower_bound(hp);
            it2 != end() && it2->first <= boundary; ++it2)
        {
            if(select_level != -1 && it2->first.level() != select_level) continue;
            int count = it2->second;
            vec.push_back(std::make_pair(it2->first, count));
            total+=count;
        }
    }

    return total;
}

int PhotonMap::extract_level(const SkyDir& dir, double radius,
                       std::vector<std::pair<HealPixel, int> >& vec,
                       int select_level, bool include_all) const
{
    bool allsky(radius>=180); // maybe use to simplify below, but seems fast
    radius *= (M_PI / 180); // convert to radians
    if (select_level == -1)
        select_level = m_minlevel; // default level to select
    int nside( 1<< select_level);
    int total(0);
    vec.clear();

    // Get pixels in select level that are within radius
    std::vector<int> v;
    Healpix hpx(nside, Healpix::NESTED, HealPixel::s_coordsys);
    hpx.query_disc(dir, radius, v);  

    // Add select level pixels to return vector
    for (std::vector<int>::const_iterator it = v.begin(); it != v.end(); ++it)
    {

        HealPixel hp(*it, select_level);
        const_iterator it2 = find(hp);
        if (it2 == end()) // Not in PhotonMap
        {
            if (include_all) // Add anyway
                vec.push_back(std::make_pair(hp, 0));
        }
        else // In PhotonMap
        {
            int count = it2->second;
            vec.push_back(std::make_pair(it2->first, count));
            total += count;
        }
    }

    return total;
}
//! Count the photons, perhaps weighted, within a given pixel.
double PhotonMap::photonCount(const HealPixel & px, bool includeChildren,
                              bool weighted) const
{

    if (!includeChildren) // No children
    {
        const_iterator it = find(px);
        if (it != end()) {
            double weight = 1 << 2*(it->first.level() - m_minlevel);

            return weighted? it->second * weight  : it->second; 
        }else  return 0;
    }else{ // Include children

        double count = 0;
        int maxLevel = m_minlevel + m_levels - 1;
        HealPixel boundary(px.lastChildIndex(maxLevel), maxLevel);
        for (PhotonMap::const_iterator it =lower_bound(px);
            it != end() && it->first <= boundary; ++it)
        {
            double weight = 1 << 2*(it->first.level() - m_minlevel);

            count += weighted? it->second * weight
                : it->second; 
        }
        return count;
    }
}

//! Count the photons within a given pixel, weighted with children.  Also return weighted direction.
double PhotonMap::photonCount(const HealPixel & px, SkyDir & NewDir) const
{


    double count = 0;
    int maxLevel = m_minlevel + m_levels - 1;
    HealPixel boundary(px.lastChildIndex(maxLevel), maxLevel);
    CLHEP::Hep3Vector v(0,0,0);

    for (PhotonMap::const_iterator it = lower_bound(px);
        it != end() && it->first <= boundary; ++it)
    {
        double weight = 1 << 2*(it->first.level() - m_minlevel);

        count += it->second * weight; 
        v +=  weight * it->second * (it->first)().dir();
    }

    NewDir = SkyDir(v);
    return count;
}

std::vector<double> PhotonMap::energyBins()const
{
    std::vector<double> result; result.push_back(m_emin);
    double eratio(exp(m_logeratio));
    for(int i = 1; i< m_levels; ++i) result.push_back(result.back()*eratio);
    return result;
}

void PhotonMap::write(const std::string & outputFile,
                                           const std::string & tablename,
                                           bool clobber) const
{
    if (clobber)
    {
        int rc = std::remove(outputFile.c_str());
        if( rc == -1 && errno == EACCES ) 
            throw std::runtime_error(std::string(" Cannot remove file " + outputFile));
    }

    // now add a table to the file
    tip::IFileSvc::instance().appendTable(outputFile, tablename);
    tip::Table & table = *tip::IFileSvc::instance().editTable( outputFile, tablename);

    table.appendField("LEVEL", "1J");
    table.appendField("INDEX", "1J");
    table.appendField("COUNT", "1J");
    table.setNumRecords(m_pixels);

    // get iterators for the Table and the Map
    tip::Table::Iterator itor = table.begin();
    const_iterator pmitor = begin();

    // now just copy
    for( ; pmitor != end(); ++pmitor, ++itor)
    {
        (*itor)["LEVEL"].set(pmitor->first.level());
        (*itor)["INDEX"].set(pmitor->first.index());
        (*itor)["COUNT"].set(pmitor->second);
    }

    // set the headers (TODO: do the comments, too)
    tip::Header& hdr = table.getHeader();
    hdr["NAXIS1"].set(3 * sizeof(long));
    hdr["EMIN"].set(m_emin); 
    hdr["ERATIO"].set(exp(m_logeratio)); 
    hdr["LEVELS"].set(m_levels); 
    hdr["MINLEVEL"].set(m_minlevel); 
    hdr["PHOTONS"].set(m_photons); 
    hdr["PIXELS"].set(m_pixels);

    // close it?
    delete &table;
}

