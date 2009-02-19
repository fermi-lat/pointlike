/** @file ComplexSkySpectrum.cxx
    @brief implement class ComplexSkySpectrum

$Header: /nfs/slac/g/glast/ground/cvs/skymaps/src/ComplexSkySpectrum.cxx,v 1.2 2008/10/21 02:47:00 burnett Exp $

*/

#include "pointlike/ComplexSkySpectrum.h"
#include "skymaps/DiffuseFunction.h"
#include "skymaps/IsotropicPowerLaw.h"

#include <cassert>
#include <stdexcept>

using namespace pointlike;


ComplexSkySpectrum::ComplexSkySpectrum(const std::string& name, const std::string& relation,size_t elem)
   :m_name(name)
   ,m_relation(relation)
   ,m_rootFunc(m_name.c_str(),m_relation.c_str(),elem)
   ,m_nElements(elem)
{
   resize(elem,std::make_pair(0,static_cast<SkySpectrum*>(0)));
};

void ComplexSkySpectrum::insert(size_t pos, const skymaps::SkySpectrum* component, double norm)
{
    if(pos>=m_nElements) throw std::runtime_error("Index error in ComplexSkySpectrum.");
    at(pos)=std::make_pair(norm, component);
}


double ComplexSkySpectrum::value(const astro::SkyDir& dir, double energy) const
{
    std::vector< std::pair<double, const skymaps::SkySpectrum*> >::const_iterator it = begin();
    for( size_t i=0; it!=end(); ++it,++i){
        if(it->second==0) throw std::runtime_error("Element not set in ComplexSkySpectrum.");
        m_rootFunc.SetParameter(i,it->first * (it->second)->value(dir, energy));
    }
    double ret=m_rootFunc.Eval(0.);
    return ret;
}

double ComplexSkySpectrum::band_value(const astro::SkyDir& dir, const skymaps::Band& band)const
{
    std::vector< std::pair<double, const skymaps::SkySpectrum*> >::const_iterator it = begin();

    for( size_t i=0; it!=end(); ++it,++i){
        if(it->second==0) throw std::runtime_error("Element not set in ComplexSkySpectrum.");
        m_rootFunc.SetParameter(i,it->first * (it->second)->band_value(dir, band));
    }
    double ret=m_rootFunc.Eval(0.);
    return ret;

}



///@brief integral for the energy limits, in the given direction
double ComplexSkySpectrum::integral(const astro::SkyDir& dir, double emin, double emax)const
{
    // estimate integral by assuming power-law from start to end
    double fa(value(dir, emin)), fb(value(dir,emax));
    double q ( 1. - log(fa/fb)/log(emax/emin) );
//    std::cout<<"ComplexSkySpectrum::integral: l="<<dir.l()<<" b="<<dir.b()<<" fa="<<fa<<" fb="<<fb
//    <<" emin="<<emin<<" emax="<<emax<<" result="<<(fa* emin * (pow(emax/emin, q)-1)/q)<<std::endl;

    return fa* emin * (pow(emax/emin, q)-1)/q;
}

std::string ComplexSkySpectrum::name()const {
    std::string ret(m_name);
    return ret;
}
