/** @file SpectralFunction.h
@brief declaration of class SpectralFunction

$Header$
*/

#ifndef pointlike_SpectralFunction_h
#define pointlike_SpectralFunction_h

#include "astro/SkyDir.h"

#include <vector>
#include <string>

namespace skymaps{
class Band;
class Exposure;

}
namespace pointlike{

    /** @class SpectralFunction
        @brief evaluate a function of energy


    */
    class SpectralFunction {
    public:
        typedef enum {PowerLaw, ExpCutoff, BrokenPowerLaw} Type;

        SpectralFunction(std::string name, const std::vector<double>& pars);

        /// @brief evaluate function itself
        double operator()(double energy)const{return value(energy);}
        /// @brief the function

        double value(double energy)const;

        /// @brief return expected counts in band
        ///
        /// @param dir the direction for exposure
        /// @param band the band object
        /// Note that this depends on exposure and effective area
        /// It will to the integral over the band
        double expected(const astro::SkyDir& dir,const skymaps::Band& band)const; 
    

        int type()const{return m_type;}
        std::string name()const{return m_name;}
        
        std::vector<double>pars()const{return m_pars;}

        void set_parameters(const std::vector<double>& pars){m_pars=pars;}

        static void set_exposures(const skymaps::Exposure* front, const skymaps::Exposure* back);

        /// @brief access to exposure object
        static const skymaps::Exposure* get_exposure(int n);
        static void set_simpson(int n);

    private:
        std::string m_name;
        int m_type;
        std::vector<double> m_pars;
        
        static int s_n;///< Simpson's rule count
        static double s_e0; 
        static double s_flux_scale;


        /// list of exposure objects
        static std::vector<const skymaps::Exposure*> s_exposures;

    };

}
#endif
