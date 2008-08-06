/** @file ExtendedLikelihood.h
    @brief declaration of class ExtendedLikelihood

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/ExtendedLikelihood.h,v 1.1 2008/06/18 01:19:23 funk Exp $

*/

#ifndef tools_ExtendedLikelihood_h
#define tools_ExtendedLikelihood_h
#ifdef OLD
#include "healpix/HealPixel.h"
#else
#include "skymaps/Band.h"
#endif
#include "astro/SkyDir.h"
#include "astro/SkyFunction.h"

#include "skymaps/PsfFunction.h"
#include "skymaps/SkySpectrum.h"
#include <TMath.h>
#include "pointlike/ExtendedSourcePDF.h"

#include <vector>
#include <utility>

namespace pointlike {
  class SkySpectrum;
  
  /** @class ExtendedLikelihood
      @brief Simple ExtendedLikelihood analysis
      
      Manage likelihood analysis of a simple model: constant background plus a point source, for data in an
      energy band such that the PSF is constant.
      
  */
  
  class ExtendedLikelihood  : public astro::SkyFunction{
    
  public:
    
    
    /** ctor
    @param data   vector of directions, weights
    @param dir    initial direction
    @param profile source profile (POINTSRC, DISK, GAUSS)
    @param src_param initial parameters of the source profile
    @param umax   [25] maximum value for the u variable
    @param diffuse
    */
    
    ExtendedLikelihood(const skymaps::Band& data,
		       const astro::SkyDir& dir, 
		       const std::string& name, const std::vector<double>& src_param, 
		       double umax, const skymaps::SkySpectrum* diffuse);
    
    ExtendedLikelihood(const skymaps::Band& data,
		       const astro::SkyDir& dir, 
		       double umax, const skymaps::SkySpectrum* diffuse);
    
    ~ExtendedLikelihood();
    
    const skymaps::Band& band()const {return m_band;}
    
    //! @return log likelihood for the signal fraction
    //! @param a value for signal fraction (default: use last fit)
    double operator()( double a=-1)const;
    
    //! @brief maximize the likelihood (minimize -log L)
    //! @return (signal_fraction, error)
    std::pair<double,double> maximize();
    std::pair<double,double> maximizeMinuit();
    std::pair<double,double> maximize(bool newStyle);
    
    int photons()const { return m_photon_count; } ///< number of photons used
    
    /// @brief First derivitive: gradient of function just evaluated  
    Hep3Vector ps_gradient() const;
    const std::vector<double>& gradient(const CLHEP::Hep3Vector& ex,
					const CLHEP::Hep3Vector& ey) const;
    
    /// @return Second derivative along arbitrary direction.
    double ps_curvature() const;
    
    /// @return value of likelihood determined by gradient calculation
    double value() const{return m_w;}
    
    /// @brief update the direction
    void setDir(const astro::SkyDir& dir, bool subset=false);
    void setDir(const astro::SkyDir& dir, const std::vector<double>& src_param,
		bool subset=false);
    
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
    
    /// check average f
    double average_F()const{ return m_averageF;}

    /// average normalized background: should be 1
    double average_b()const{ return m_avb;}

    void changepsf(){}; // note not implemented
    double gamma()const { return m_psf.gamma();}
    void   setgamma(double gamma) {m_psf.setGamma(gamma);}

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

    double psfIntegral()const{ return m_fint;}

    void recalc(bool subset=true);
    void reload(bool subset=true);

    const std::vector<double>& residual() const{ return m_vloglike; };
    const std::vector< std::vector<double> >& jacobian() const{ return m_vjacobian; };


    /// @brief access to the diffuse component
    const skymaps::SkySpectrum* diffuse() const;

    /// @brief set the diffuse component
    void setDiffuse(skymaps::SkySpectrum* diff);

    static double tolerance();
    static void setTolerance(double tol);
    static double defaultUmax();
    static void setDefaultUmax(double umax);
    static double defaultRoI();
    static void setDefaultRoI(double roi);


  private: 
    
    static double s_defaultUmax;
    static double s_defaultRoI;
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
    
    ExtendedSourcePseudoPSF m_psf;
    std::vector<double> m_src_param;
    
    double m_sigma;
    double m_alpha, m_sigma_alpha; ///< current fit value, error
    mutable double m_curv;  // saved curvature from gradient calculation
    double m_background;  ///< expected background (negative: no estimate)
    double m_umax; ///< maximum value of u, for selection of data, fits
    double m_avu, m_avb;
    double m_emin, m_emax; ///< energy range for this object
    double m_F;       ///< integral of PSF over u
    
    const skymaps::SkySpectrum * m_diffuse; ///< background distribution to use
    mutable std::vector<double> m_gradient;
    
    mutable std::vector<double> m_vloglike;
    mutable std::vector< std::vector<double> > m_vjacobian;
    
    class NormalizedBackground; // forward declaration of helper class
    NormalizedBackground* m_back; ///< instance of private helper class set when direction is
    
  };
}
#endif

