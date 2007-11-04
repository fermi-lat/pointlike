/** @file CompositeSkySpectrum.h
    @brief declare class CompositeSkySpectrum

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/CompositeSkySpectrum.h,v 1.1 2007/11/01 21:33:33 burnett Exp $

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
    /// @brief ctor 
    CompositeSkySpectrum(){}

    /// @brief add a new diffuse component
    /// @diffuse pointer to a SkySpectrum
    /// @norm[1] the normalization factor
    void add(const pointlike::SkySpectrum* diffuse, double norm=1.);

    ///@brief return differential value 
    ///@param energy energy in MeV
    virtual double value(const astro::SkyDir& dir, double energy)const;
   
    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    ///@brief name for identification
    virtual std::string name()const;
private:

};

}
#endif
