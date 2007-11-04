/** @file DiffuseFunction.h
    @brief declare class DiffuseFunction

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/DiffuseFunction.h,v 1.9 2007/11/01 21:33:33 burnett Exp $

*/
#ifndef pointlike_DiffuseFunction_h
#define pointlike_DiffuseFunction_h
#include "pointlike/SkySpectrum.h"
#include "astro/SkyFunction.h"
#include "astro/SkyDir.h"
#include "map_tools/SkyImage.h"
//#include "tools/Aeff.h"

#include <vector>
#include <cassert>

namespace pointlike {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class DiffuseFunction
    @brief a SkyFunction that adapts a diffuse map. also includes extragal diffuse

*/

class DiffuseFunction : public pointlike::SkySpectrum {
public:

    /** @brief ctor that reads a FITS cube represention of the diffuse, with multple layers for
         eneries. 
         @param diffuse_cube_file Name of file. It must have an extension ENERGIES with the corresponding energy valuse
         @param energy[1000] initial energy for the SkyFunction 

    */
    DiffuseFunction(std::string diffuse_cube_file, double energy=1000.);
 
    virtual ~DiffuseFunction();

    double extraGal(double energy) const; ///< access to the extra galactic component

    ///@brief interpolate table 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const{return m_name;}


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


    std::vector<double> m_energies; ///< list of energies
    int layer(double e)const;

    double energy_bin(int k) const;
    map_tools::SkyImage m_data;
    std::string m_name;
    double m_emin, m_emax;

};
} // namespace pointlike
#endif

