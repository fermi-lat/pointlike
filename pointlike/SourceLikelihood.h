/** @file SourceLikelihood.h

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SourceLikelihood.h,v 1.22 2008/02/19 21:00:32 burnett Exp $
*/

#ifndef tools_SourceLikelihood_h
#define tools_SourceLikelihood_h

#include "pointlike/ExtendedLikelihood.h"

#include "skymaps/SkySpectrum.h"

#include "astro/SkyDir.h"
#ifdef OLD
#include "healpix/HealPixel.h"
#else
#include "skymaps/Band.h"
#include "skymaps/BinnedPhotonData.h"
#endif
#include <iostream>
#include <map>
#include <vector>

namespace embed_python { class Module; }
namespace skymaps { 
  class BinnedPhotonData;
  class CompositeSkySpectrum;
}

namespace pointlike {
  
  
  /** @class SourceLikelihood
      @brief manage a set of ExtendedLikelihood objects, one for each energy band
      
      Note that it is a map of the ExtendedLikelihood objects, with the key an index to a corresponding Band 
      object, as maintained by the corresponding BinnedPhotonData.
      
      
  */
  class SourceLikelihood : public  std::vector<ExtendedLikelihood*>, 
			   public skymaps::SkySpectrum {
  public:
    

    /** ctor
	@param data   source of the data to fit
	@param name   source name for printout
	@param dir    initial direction
    */
    SourceLikelihood(const skymaps::BinnedPhotonData& data,
		     std::string name,
		     const astro::SkyDir& dir,
		     const std::string type="point",
		     const std::vector<double>& src_param= std::vector<double>(0));
    
    ~SourceLikelihood();

    //! fit to signal fractions 
    /// @return total TS
    double  maximize();

    //! change the current direction -- resets data and refits
    void setDir(const astro::SkyDir& dir,bool subset=false);
    void setDir(const astro::SkyDir& dir, const std::vector<double>& srcparam, 
		bool subset=false);

    /// @return the gradient, summed over all bands
    const CLHEP::Hep3Vector& ps_gradient() const;

    ///@return the curvature, summed over all bands
    double ps_curvature() const;

    const std::vector<double> gradient() const;

    /// @brief 
    void printSpectrum();

    /// @brief perform localization fit, maximizing joint likelihood
    /// @param skip [0] number of bands to skip
    /// @return error circle radius (deg) or large number corresponding to error condition
    double localize(int skip);
    double localizeMinuit(int skip);
    double localizeMinpack(int skip);

    /// @brief invoke localize with skip values from skip1 to skip 2 or until good fit
    double localize(int skip1, int skip2);

    /// @brief localate with iteration to refit the levels, using parameters set in ctor
    double localize();

    std::string name()const{return m_name;}

    const astro::SkyDir& dir()const{return m_dir;}
    
    double TS()const { return m_TS; } 
    
    double logL(){ return m_loglike;}

    double errorCircle()const{return  sqrt(1./ps_curvature())*180/M_PI;}

    void set_ostream(std::ostream& out){m_out=&out;}

    static void set_verbose(bool verbosity=true); //{s_verbose=verbosity;}

    static bool verbose(); //{return s_verbose;}

    ///! implement the SkyFunction interface
    ///@brief return differential value 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    /// @brief set all parameters using the embedded python module
    static void SourceLikelihood::setParameters(const embed_python::Module& par);
    
    /// @brief set radius for individual fits
    static void setDefaultUmax(double umax);
    
    //     /// @brief access to the sigma (radians) used for the individual ExtendedLikelihood objects
    //     double sigma(int level)const;

    std::string type()const {return m_type;};
    
    std::vector<double> sourceParameters() const { return m_sourceParameters; };
    std::vector<double> sourceParErrors() const   { return m_sourceParErrors; };
    double errorX() const {return m_errorX;};
    double errorY() const {return m_errorY;};
    
    //     static double set_gamma_level(int level, double v);
    
    //     static double set_sigma_level(int level, double v);
    
    ///! Set the global diffuse background function, return current value 
    static  skymaps::SkySpectrum* set_diffuse( skymaps::SkySpectrum* diffuse, 
					       double exposure = 1.0);
    
    ///! add a point source fit to the background for subsequent fits
    void addBackgroundPointSource(const SourceLikelihood* fit);
    
    ///! remove all such
    void clearBackgroundPointSource();
    
//     ///! @brief recalculate likelihoods using any static changes made to parameters
//     void recalc(int level);

    int& useMinuit(){ return s_useMinuit;};
    int  useMinuit() const { return s_useMinuit;};

    ///! access to background model 
    const skymaps::SkySpectrum * background()const;

//     /// @brief get the starting, ending levels used
//     static int minlevel();
//     static int maxlevel();
//     static void set_levels(int minlevel, int maxlevel=13);

    void setMinuitLevelSkip(int skip){s_minuitLevelSkip=skip;};
    int  minuitLevelsToSkip() const {return s_minuitLevelSkip;};
    
    /// @brief set the integration tolerance for the background, return present value
    static double set_tolerance(double tol);
    
    /// @brief set the range of energy to fit
    static void set_energy_range(double emin, double emax=1e6);
    
    /// @brief special display function
    /// @param dir direction
    /// @param energy selects energy band
    /// @mode 0: same as the operator; 1: data; 2: background; 3:fit; 4:residual
    ///     
    double display(const astro::SkyDir& dir, double energy, int mode)const;
    std::vector<double> energyList() const;

  private:
    void setup(const skymaps::BinnedPhotonData& data);
    std::string m_name;
    astro::SkyDir m_dir; ///< common direction
    double m_dir_sigma;  ///< error circle from fit (radians)
    double m_loglike;    ///< total loglike
    double m_TS;         ///< total TS value

    std::string m_type;
    unsigned int m_npar;
    std::vector<double> m_sourceParameters;
    std::vector<double> m_sourceParErrors;
    double m_errorX;
    double m_errorY;
    
    std::ostream * m_out;
    std::ostream& out()const{return *m_out;}
    mutable CLHEP::Hep3Vector m_gradient; ///< current gradient
    
    static skymaps::SkySpectrum* s_diffuse; ///< global diffuse used by all PSL objects
    
    skymaps::CompositeSkySpectrum * m_background;  ///< background spectrum to use
    
    static double s_emin, s_emax, s_minalpha, s_TSmin, s_tolerance, 
      s_maxstep; //
    static int s_useMinuit;
    static int s_minuitLevelSkip;
    static int s_skip1, s_skip2, s_itermax, s_verbose;
    
  };
  
    /** @class SLdisplay

*/
  class SLdisplay :  public skymaps::SkySpectrum {
  public:
    SLdisplay(const SourceLikelihood & psl, int mode);
    virtual double value(const astro::SkyDir& dir, double e)const;
    
    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;
    
    std::string name()const;
  private:
    const SourceLikelihood& m_psl;
    int m_mode;
    
    
  };
    
};
#endif

