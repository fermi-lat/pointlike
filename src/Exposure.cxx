/** @file Exposure.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Exposure.cxx,v 1.2 2007/11/22 02:38:49 burnett Exp $
*/

#include "pointlike/Exposure.h"
#include "healpix/HealpixArrayIO.h"

#include <stdexcept>

namespace{  // anonymous namespace for helper classes
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class IrfAeff
@brief function class implements effective area, as adapter to irfInterface::IAeff
*/
class IrfAeff {
public:
    /**
    @param cutoff limit for cos(theta)
    */
    IrfAeff( double cutoff=0.25)
        : m_cutoff(cutoff)
    {}

    double operator()(double costh) const
    {
           return costh<m_cutoff? 0 : (costh-m_cutoff)/(1.-m_cutoff);
    }
    double m_cutoff;
} aeff;


}  // anonymous namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

using namespace pointlike;

Exposure::Exposure(const std::string& fits_file, const std::string& tablename)
: m_filename(fits_file)
, m_exposure( healpix::HealpixArrayIO::instance().read(fits_file, tablename) )
{
    
}

Exposure::~Exposure()
{}


double Exposure::value(const astro::SkyDir& dir, double)const
{
    const healpix::CosineBinner& binner = m_exposure[dir];
    double val( binner( aeff ) );
    return val;
}

///@brief integral for the energy limits, in the given direction
double Exposure::integral(const astro::SkyDir& dir, double , double )const
{
    return value(dir, 1000.);
}

std::string Exposure::name()const
{
    return "Exposure: "+m_filename;
}
