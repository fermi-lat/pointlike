/** @file CompositeSkySpectrum.cxx
    @brief implement class CompositeSkySpectrum

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/CompositeSkySpectrum.cxx,v 1.1 2007/11/04 22:11:32 burnett Exp $

*/

#include "pointlike/CompositeSkySpectrum.h"
#include <cassert>

using namespace pointlike;

CompositeSkySpectrum::CompositeSkySpectrum(const pointlike::SkySpectrum* component, double norm)
{
    if( component!=0 ) add(component, norm);
}

void CompositeSkySpectrum::add(const pointlike::SkySpectrum* component, double norm)
{
    push_back(std::make_pair(norm, component));
}


double CompositeSkySpectrum::value(const astro::SkyDir& dir, double energy) const
{
    double ret(0);
    std::vector< std::pair<double, const pointlike::SkySpectrum*> >::const_iterator it = begin();

    for( ; it!=end(); ++it){
        ret+= it->first * (it->second)->value(dir, energy);
    }
    return ret;
}



///@brief integral for the energy limits, in the given direction
double CompositeSkySpectrum::integral(const astro::SkyDir& dir, double emin, double emax)const
{
    double ret(0);
    std::vector< std::pair<double, const pointlike::SkySpectrum*> >::const_iterator it = begin();

    for( ; it!=end(); ++it){
        ret+= it->first * (it->second)->integral(dir, emin, emax);
    }
    return ret;

}

std::string CompositeSkySpectrum::name()const {
    std::string ret(m_name);
    if( ret.empty()){
        std::vector< std::pair<double, const pointlike::SkySpectrum*> >::const_iterator it = begin();
        for( ; it!=end(); ++it){
            ret+= ", " + (it->second)->name();
        }
    }
    return ret;
}
