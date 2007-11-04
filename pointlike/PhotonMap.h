/** @file PhotonMap.h
    @brief declare class PhotonMap

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/PhotonMap.h,v 1.9 2007/11/01 21:33:33 burnett Exp $

*/
#ifndef pointlike_PhotonMap_h
#define pointlike_PhotonMap_h

#include "pointlike/SkySpectrum.h"
#include "map_tools/PhotonMap.h"

#include <string>

namespace pointlike {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class PhotonMap
    @brief a SkySpectrum that adapts the class map_tools::PhotonMap. 

*/

class PhotonMap : public pointlike::SkySpectrum, public map_tools::PhotonMap {
public:
    PhotonMap(const std::string & inputFile,const std::string & tablename="PHOTONMAP" );
    ~PhotonMap(){};

    
    ///@brief data value for bin with given energy 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    /// Assume that request is for an energy bin.
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const{return m_name;}

    ///@ allow different name
    void setName(const std::string name){m_name = name;}

private:
    std::string m_name; ///< use name of file as descriptive name
};


}

#endif
