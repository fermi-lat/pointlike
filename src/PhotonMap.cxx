/** @file PhotonMap.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/PhotonMap.cxx,v 1.9 2007/11/01 21:33:33 burnett Exp $
*/

#include "pointlike/PhotonMap.h"
#include <stdexcept>
#include <cmath>

using namespace pointlike;


PhotonMap::PhotonMap(const std::string & inputFile,const std::string & tablename )
: map_tools::PhotonMap(inputFile, tablename)
, m_name(inputFile)
{
}

   ///@brief data value for bin with given energy 
    ///@param e energy in MeV
double PhotonMap::value(const astro::SkyDir& dir, double energy)const
{
    astro::Photon gam(dir, energy, 0); // just for the next interface :-)
    PhotonMap* self = const_cast<PhotonMap*>(this);  // pixel not const
    astro::HealPixel pix( self->pixel(gam) );
    return photonCount(pix) ;
}

    ///@brief integral for the energy limits, in the given direction
    /// Assume that request is for an energy bin.
double PhotonMap::integral(const astro::SkyDir& dir, double a, double b)const
{
    return value(dir, sqrt(a*b));
}
