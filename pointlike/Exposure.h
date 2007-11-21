/** @file Exposure.h
    @brief declare class Exposure

$Header:$

*/
#ifndef pointlike_Exposure_h
#define pointlike_Exposure_h

#include "pointlike/SkySpectrum.h"


namespace pointlike {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Exposure
    @brief a SkyFunction that represents the exposure over the sky


    (dummy for now)
*/

class Exposure : public pointlike::SkySpectrum {
public:
    Exposure();
    ~Exposure();

    ///@brief a single energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const;


private:
};


}
#endif

