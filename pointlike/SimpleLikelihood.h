/** @file SimpleLikelihood.h
    @brief declaration of class SimpleLikelihood

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SimpleLikelihood.h,v 1.26 2008/05/08 05:42:17 mar0 Exp $

*/

#ifndef tools_Likelihood_h
#define tools_Likelihood_h
#ifdef OLD
#include "healpix/HealPixel.h"
#else
#include "skymaps/Band.h"
#endif
#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"

#include "skymaps/PsfFunction.h"
#include "skymaps/SkySpectrum.h"

#include <vector>
#include <utility>

namespace pointlike {
class SkySpectrum;

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
    @param umax   [25] maximum value for the u variable
    @param diffuse

    */
    SimpleLikelihood(const skymaps::Band& data,
        const astro::SkyDir& dir, 
        double umax 
        ,const skymaps::SkySpectrum* diffuse);

    ~SimpleLikelihood();

    const skymaps::Band& band()const {return m_band;}

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
    void setDir(const astro::SkyDir& dir, bool subset=false);

    double alpha()const { return m_alpha; }
 
    void setalpha(double alpha) {m_alpha=alpha;}

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

    double feval(double k);
    double geval(double k);

    void changepsf(){}; // note not implemented
    void setgamma(double gamma) {m_psf=skymaps::PsfFunction(gamma);}

    /// @brief implement the SkyFunction interface
    /// @return the events/sr corresponding to the solution
    double operator()(const astro::SkyDir& dir)const;

    /// @brief special display function
    /// @param dir direciton
    /// @mode 0: same as the operator; 1: data; 2: background; 3:fit; 4:residual
    ///     
    double display(const astro::SkyDir& dir, int mode)const;

    /// @brief access to the effective sigma (radians)  used for the fits
    double sigma()const{ return m_sigma;}
    void setsigma(double sigma) {m_sigma=sigma;}
    double gamma()const{ return m_psf.gamma();} 

    void recalc(bool subset=true);
    void reload(bool subset=true);

    /// @brief access to the diffuse background component 
    const skymaps::SkySpectrum* diffuse() const;

    /// @brief set the diffuse component
    void setDiffuse(skymaps::SkySpectrum* diff);

    static double tolerance();
    static void setTolerance(double tol);
    static double defaultUmax();
    static void setDefaultUmax(double umax);


private:


    static double s_defaultUmax;
    static double s_tolerance; // for integral

    
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
    const skymaps::Band& m_band;
    typedef std::vector<std::pair<astro::SkyDir, int> > PixelList;
    PixelList m_vec;
    
    //! simplified set with function or distances from m_dir 
    std::vector<std::pair<double, int> > m_vec2;  //stores <log-like,nphotons>
    std::vector<int> m_vec4; //stores subset healpix indices
    double m_averageF;
    skymaps::PsfFunction m_psf;
    double m_sigma;
    double m_alpha, m_sigma_alpha; ///< current fit value, error
    mutable double m_curv;  // saved curvature from gradient calculation
    double m_background;  ///< expected background (negative: no estimate)
    double m_umax; ///< maximum value of u, for selection of data, fits
    double m_avu, m_avb;
    double m_emin, m_emax; ///< energy range for this object
    double m_F;       ///< integral of PSF over u

    const skymaps::SkySpectrum * m_diffuse; ///< background distribution to use

    class NormalizedBackground; // forward declaration of helper class
    NormalizedBackground* m_back; ///< instance of private helper class set when direction is
    
};
}
#endif

