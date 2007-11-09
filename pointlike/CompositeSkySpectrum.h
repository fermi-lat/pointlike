/** @file CompositeSkySpectrum.h
    @brief declare class CompositeSkySpectrum

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/CompositeSkySpectrum.h,v 1.1 2007/11/04 22:11:32 burnett Exp $

*/
#ifndef pointlike_CompositeSkySpectrum_h
#define pointlike_CompositeSkySpectrum_h

#include "pointlike/SkySpectrum.h"

#include <vector>
#include <algorithm> 
namespace pointlike {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class CompositeSkySpectrum
    @brief a SkySpectrum that combines weighted SkySpectrum objects

*/

class CompositeSkySpectrum: public pointlike::SkySpectrum, public  std::vector< std::pair<double, const pointlike::SkySpectrum*> >{
public:
    /// @brief ctor - default, or same as an add
    CompositeSkySpectrum(const pointlike::SkySpectrum* diffuse=0, double norm=1.);

    /// @brief add a new diffuse component
    /// @param component pointer to a SkySpectrum
    /// @param norm[1] the normalization factor
    void add(const pointlike::SkySpectrum* component, double norm=1.);

    ///@brief return differential value 
    ///@param energy energy in MeV
    virtual double value(const astro::SkyDir& dir, double energy)const;
   
    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    ///@brief name for identification
    /// default is to make a list of the names of the components
    virtual std::string name()const;

    ///@brief override name
    void setName(const std::string & name){m_name=name;}
private:
    std::string m_name;

};

}
#endif
