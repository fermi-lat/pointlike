/** @file ComplexSkySpectrum.h
    @brief declare class ComplexSkySpectrum

$Header: /nfs/slac/g/glast/ground/cvs/skymaps/skymaps/ComplexSkySpectrum.h,v 1.2 2008/10/21 02:47:00 burnett Exp $

*/
#ifndef skymaps_ComplexSkySpectrum_h
#define skymaps_ComplexSkySpectrum_h

#include "skymaps/SkySpectrum.h"
#include "TF1.h"

#include <vector>
#include <algorithm> 



namespace pointlike {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class ComplexSkySpectrum
    @brief a SkySpectrum that combines weighted SkySpectrum objects

*/

class ComplexSkySpectrum: public skymaps::SkySpectrum,
    public  std::vector< std::pair<double, const skymaps::SkySpectrum*> >{
public:
    /// @brief ctor - default, or same as an add

    ComplexSkySpectrum(const std::string& name, const std::string& relation, size_t elem);

    /// @brief add a new diffuse component
    /// @param component pointer to a SkySpectrum
    /// @param norm[1] the normalization factor
    void insert(size_t pos,const skymaps::SkySpectrum* component, double norm=1.);

    ///@brief return differential value 
    ///@param energy energy in MeV
    virtual double value(const astro::SkyDir& dir, double energy)const;

    ///@brief use a band to select interval. A subclass can use it to select subclass
    /// base class implemention integrates over band
    virtual double band_value(const astro::SkyDir& dir, const skymaps::Band& band)const;
   
    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    ///@brief name for identification
    /// default is to make a list of the names of the components
    virtual std::string name()const;

    ///@brief override name
    void setName(const std::string & name){m_name=name;}
private:
    std::string m_name;
    std::string m_relation;
    mutable TF1 m_rootFunc;
    size_t m_nElements;
};

}
#endif
