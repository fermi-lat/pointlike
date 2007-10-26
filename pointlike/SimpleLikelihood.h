/** @file SimpleLikelihood.h
    @brief declaration of class SimpleLikelihood

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SimpleLikelihood.h,v 1.8 2007/10/06 17:18:55 burnett Exp $

*/

#ifndef tools_Likelihood_h
#define tools_Likelihood_h

#include "astro/HealPixel.h"
#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"

#include "pointlike/PsfFunction.h"

#include <vector>
#include <utility>

namespace pointlike {
class DiffuseFunction;

/** @class SimpleLikelihood
@brief Simple SimpleLikelihood analysis

Manage likelihood analysis of a simple model: constant background plus a point source, for data in an
energy band such that the PSF is constant.

*/

class SimpleLikelihood  : public astro::SkyFunction{
public:
    /** ctor
    @param data   vector of directions, weights
    @param dir    initial direction
    @param gamma  power-law value for PsfFunction
    @param sigma  scale factor to define u
    @param background [-1] background density, events/solid angle (negative to not use)
    @param umax   [25] maximum value for the u variable
    */
    SimpleLikelihood(const std::vector<std::pair<astro::HealPixel,int> >& data,
        const astro::SkyDir& dir, 
        double gamma, double sigma, double background=-1, double umax=s_defaultUmax, double energy=1000);


    //! @return log likelihood for the signal fraction
    //! @param a value for signal fraction (default: use last fit)
    double operator()( double a=-1)const;

    //! @brief maximize the likelihood (minimize -log L)
    //! @return (signal_fraction, error)
    std::pair<double,double> maximize();

           
    int photons()const { return m_photon_count; } ///< number of photons used

    /// @brief First derivitive: gradient of function just evaluated  
    Hep3Vector gradient() const;

    /// @return Second derivative along arbitrary direction.
    double curvature() const;

    /// @return value of likelihood determined by gradient calculation
    double value() const{return m_w;}

    /// @brief update the direction
    void setDir(const astro::SkyDir& dir);

    double alpha()const { return m_alpha; }

    double sigma_alpha() const { return m_sigma_alpha;}

    /// @brief return the test statistic
    double TS(double alpha = -1)const;

    /// @brief set a known background density, which enables use of Poisson. Negative for not used
    void setBackgroundDensity(double b){m_background=b;}

    /// @return the  negative log likelihood for the poisson, or zero if 
    /// no background estimate
    double poissonLikelihood(double a)const;

    /// @return signal
    double signal(double a=-1)const{return (a<0?m_alpha:a)*m_photon_count;}

    /// @return background estimte
    double background()const{return m_background*solidAngle();}

    /// @return the solid angle in sr for the circular aperature used for analysis
    double solidAngle()const;

    /// set/retrieve the umax parameter
    double umax()const {return m_umax;}

    /// check average u
    double average_u()const{ return m_avu;}

    /// average normalized background: should be 1
    double average_b()const{ return m_avb;}
    static double s_defaultUmax;

    double feval(double k);

    double kcurvature(double k);

    void setEnergy(double e){m_energy = e;}

    /// @brief implement the SkyFunction interface
    /// @return the events/pixel corresponding to the solution
    double operator()(const astro::SkyDir& dir)const;

    static DiffuseFunction* s_diffuse;
    static double s_tolerance; // for integral
private:

    //! @brief a quick estimate of the signal fraction
    //! @return the value of of the signal fraction
    double estimate() const;

    /// @brief derivatives of log poisson with respect to alpha for case with known background
    /// If background not specified, return (0,0)
    std::pair<double,double>poissonDerivatives(double x);

    astro::SkyDir m_dir;
    mutable int    m_photon_count; // total number of photons
    double m_fint, m_fint2; //integral of f, f^2


    mutable double m_w;      // likelihood from gradient

    //! vector of healpixels and the number of photons in each
    const std::vector<std::pair<astro::HealPixel,int> >& m_vec;
    //! simplified set with function or distances from m_dir 
    std::vector<std::pair<double, int> > m_vec2;
    std::vector<double> m_vec3; //storage of u values for fast Likelihood recalculation
    double m_averageF;
    pointlike::PsfFunction m_psf;
    double m_sigma;
    double m_alpha, m_sigma_alpha; ///< current fit value, error
    mutable double m_curv;  // saved curvature from gradient calculation
    double m_background;  ///< expected background (negative: no estimate)

    double m_umax; ///< maximum value of u, for selection of data, fits
    double m_avu, m_avb;
    double m_energy; ///< median energy
    
};
}
#endif

