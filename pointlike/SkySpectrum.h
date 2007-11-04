/** @file SkySpectrum.h
    @brief declare class SkySpectrum

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SkySpectrum.h,v 1.2 2007/11/01 21:43:41 burnett Exp $

*/
#ifndef pointlike_SkySpectrum_h
#define pointlike_SkySpectrum_h

#include "astro/SkyFunction.h"

namespace pointlike {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class SkySpectrum
    @brief a SkyFunction that implements a spectrum at each point

    It is virtual: subclasses must implement the differential and integralfunctions value and integral
*/

    class SkySpectrum: public astro::SkyFunction {
public:
    SkySpectrum(double energy=1000);
    virtual ~SkySpectrum(){}

    /// @brief set an energy to evaluate the differential distribution with the basic operator()
    void setEnergy(double e)const;

    /// @brief set the energy range evaluation 
    /// change evaluation to a range of energies
    void setEnergyRange(double emin, double emax)const;

    ///@brief return differential value 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const=0;
   
    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const=0;

    ///@brief name for identification
    virtual std::string name()const=0;

    /// @brief Implement the SkyFunction interface, convenient for creating a SkyImage.
    /// @return value for given direction and currently selected energy, or integral over the energy range 
    double operator()(const astro::SkyDir& dir)const; 

    /// @brief functor that allows energy as well
    double operator()(const astro::SkyDir& dir, double energy)const;
    
    /// @brief functor that returns an integral over the energies as well
    double operator()(const astro::SkyDir& dir, double emin, double emax)const;

    ///! average, for the given energy or energy range, about the direction and cone angle(radians)
    double average(const astro::SkyDir& dir, double angle, double tolerance)const;


private:
    mutable double m_energy;
    mutable double m_emin, m_emax; ///< range for integral
    mutable bool m_use_range; ///< true: evaluate specified energy range; false: value at energy

    double level_ave(const astro::SkyDir& dir, double angle, int level) const;


};

}
#endif
