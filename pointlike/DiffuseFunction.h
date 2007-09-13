/** @file DiffuseFunction.h

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/DiffuseFunction.h,v 1.2 2007/09/09 19:50:06 burnett Exp $

*/
#ifndef pointlike_DiffuseFunction_h
#define pointlike_DiffuseFunction_h
#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"
#include "map_tools/SkyImage.h"
//#include "tools/Aeff.h"

#include <vector>
#include <cassert>
namespace pointlike {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class DiffuseFunction
    @brief a SkyFunction that adapts a diffuse map. also computes extragal diffuse

*/

class DiffuseFunction : public astro::SkyFunction {
public:
    DiffuseFunction(std::string diffuse_cube_file, double energy=1000.)
        : m_data(diffuse_cube_file)
        , m_fract(0)
        , m_emin(0), m_emax(0)
    {
        setEnergy(energy);
    }

    double extraGal(double energy) const; ///< access to the extra galactic component

    void setEnergy(double e)const;

    void setLayer(int layer){ 
        m_layer=layer;
        m_fract=0;
    }
    /// @brief set the energy range for integrals: note const
    void setEnergyRange(double emin, double emax)const{
        m_emin = emin; m_emax=emax;
        assert( emin>0 && emax/emin<3); // otherwise integral not (now) valid
    }
    
    /// Implement SkyFunction
    ///@return interpolation of the table for given direction and current energy 
    /// 
    double operator()(const astro::SkyDir& dir)const; 

    /// functor that allows energy as well
    double operator()(const astro::SkyDir& dir, double energy)const
    {
        setEnergy(energy); return (*this)(dir);
    }
    

    ///@return integral for the energy limits, in the given direction
    double integral(const astro::SkyDir& dir, double a, double b)const;

    ///! average, for the given energy, about the direction and cone angle (in radians).  level is healpix level for pixelization.
    double average(const astro::SkyDir& dir, double angle, int level = 9)const;

#if 0 // not implemented yet
    ///@return integral for the energy limits, over the function, in the given direction
    double integral(const astro::SkyDir& dir, const Aeff& f, double a, double b)const;
#endif


    //-------------------------------------------------------
    /** @brief set vector of values for set of energy bins
        @param dir  direction 
        @param energies vector of the bin edges.
        @return result vector of the values
    */
    std::vector<double> integral(const astro::SkyDir& dir, 
        const std::vector<double>&energies)const;

    /// @return number of layers
    /// @todo: get number from file
    int layers()const { return 17;}


private:
    static double s_emin;
    mutable double m_energy;
    int layer(double e)const;
    static double h(double r, double alpha);

    static double energy_bin(int k);
    map_tools::SkyImage m_data;
    mutable int m_layer;
    mutable double m_fract; ///< current fractional
    mutable double m_emin, m_emax; ///< range for integral
};
} // namespace pointlike
#endif

