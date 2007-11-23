/** @file Exposure.h
    @brief declare class Exposure

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Exposure.h,v 1.1 2007/11/21 07:00:39 burnett Exp $

*/
#ifndef pointlike_Exposure_h
#define pointlike_Exposure_h

#include "pointlike/SkySpectrum.h"
#include <string>
#include "healpix/HealpixArrayIO.h"
#include "healpix/CosineBinner.h"


namespace pointlike {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class Exposure
    @brief a SkyFunction that represents the exposure over the sky


    (dummy for now)
*/

class Exposure : public pointlike::SkySpectrum {
public:
    Exposure(const std::string & fits_file, const std::string& tablename="Exposure");
    ~Exposure();

    ///@brief a single energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const;


private:
    std::string m_filename;
    healpix::HealpixArray<healpix::CosineBinner> m_exposure;
    
};


}
#endif

