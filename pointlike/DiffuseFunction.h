/** @file DiffuseFunction.h

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/DiffuseFunction.h,v 1.6 2007/10/26 00:27:19 burnett Exp $

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
    
    /// Implement the SkyFunction interface, convenient for creating a SkyImage.
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

    ///! average, for the given energy, about the direction and cone angle (in radians).  
    /// @param dir direction
    /// @param cone half-angle, radians
    /// @param tolerance relative error 
    /// Note that if the tolerance is > 0.5, the value returned will be that at the center
    double average(const astro::SkyDir& dir, double angle, double tolerance = 1e-3)const;


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
    int layers()const { return m_data.layers();}

    /// @brief access to the contained SkyImage
    const map_tools::SkyImage& image()const { return m_data;}



private:
    double level_ave(const astro::SkyDir& dir, double angle, int level) const;

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

