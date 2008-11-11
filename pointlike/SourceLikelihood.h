/** @file SourceLikelihood.h

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SourceLikelihood.h,v 1.9 2008/10/20 23:40:12 markusa Exp $
*/

#ifndef tools_SourceLikelihood_h
#define tools_SourceLikelihood_h

#include "pointlike/ExtendedLikelihood.h"

#include "skymaps/SkySpectrum.h"

#include "astro/SkyDir.h"
#include "skymaps/Band.h"
#include "skymaps/BinnedPhotonData.h"

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
    SourceLikelihood(skymaps::BinnedPhotonData& data,
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

    const std::vector<double> gradient() const;

    /// @brief 
    void printSpectrum();

    /// @brief perform localization fit, maximizing joint likelihood
    /// @param skip [0] number of bands to skip
    /// @return error circle radius (deg) or large number corresponding to error condition
    double localizeMinuit();

    /// @brief invoke localize until good fit
    double localize();

    /// @brief localate with iteration to refit the levels, using parameters set in ctor
    double fit();

    std::string name()const{return m_name;}

    const astro::SkyDir& dir()const{return m_dir;}
    
    double TS()const { return m_TS; } 
    double TS(int band) const;
    double alpha(int band) const;

    double logL(){ return m_loglike;}

    void set_ostream(std::ostream& out){m_out=&out;}

   ///@static get/set functions for SourceLikelihood parameters to be used in the python interface 
    ///@param e energy in MeV

    static void setVerbose(bool verbosity=true) {s_verbose=verbosity;};
    static bool verbose() {return s_verbose;};

    ///! implement the SkyFunction interface
    ///@brief return differential value 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    /// @brief set all parameters using the embedded python module
    static void SourceLikelihood::setParameters(const embed_python::Module& par);
    
    /// @brief set radius for individual fits
    static void setDefaultUmax(double umax){pointlike::ExtendedLikelihood::setDefaultUmax(umax); };
    static double defaultUmax(double umax) { return pointlike::ExtendedLikelihood::defaultUmax(); };

    /// @brief set the integration tolerance for the background, return present value
    static void setTolerance(double tol){  pointlike::ExtendedLikelihood::setTolerance(tol);};
    static double tolerance()   { return pointlike::ExtendedLikelihood::tolerance();};

    static void setDefaultRoI(double roi){  pointlike::ExtendedLikelihood::setDefaultRoI(roi);};
    static double defaultRoI()  { return pointlike::ExtendedLikelihood::defaultRoI();};
    
    /// @brief set the range of energy to fit
    static void setEnergyRange(double emin,double emax=1e6){ s_emin=emin; s_emax=emax;};
    static void set_energy_range(double emin, double emax=1e6) {s_emin=emin; s_emax=emax;};///compatibility function
    static double emin() { return s_emin;};
    static double emax() { return s_emax;};

    static int nParameter() { return s_npar;};

    /// @brief set minimum alpha
    static void setMinAlpha(double alpha){ s_minalpha=alpha; };
    static void set_min_alpha(double alpha){ s_minalpha=alpha; };
    static double minAlpha() { return s_minalpha; };

    /// @brief set minimum TS
    static void setMinTS(double ts){ s_TSmin=ts; };
    static double minTS() { return s_TSmin; };

    /// @brief set fit accuracy
    static void setFitAccuracy(double acc){ s_accuracy=acc; };
    static double fitAccuracy() { return s_accuracy; };

    void fixPosition();
    void fixPosition(const astro::SkyDir& dir);
    void freePosition();

    void fixExtension();
    void fixExtension(std::vector<double> radii);
    void freeExtension();

    /// @brief set minuit mode to be used by vector of strings: 
    /// @param modes can be (MIGRAD|SIMPLEX),(MINOS|HESSE),(GRAD,NOGRAD)
    static void setMinuitMode(const std::vector<std::string> modeVec);
    static std::vector<std::string> minuitMode() ;

    /// @brief set gamma parameter of psf for energy bins
    /// @param selection: 'front' or 'back'
    /// @param gamma: vector with gamma constants 
    static void setGamma(const std::string selection,const std::vector<double> gamma);
    static std::vector<double> gamma(const std::string selection) ;

    /// @brief set sigma parameter of psf for energy bins
    /// @param selection: 'front' or 'back'
    /// @param sigma vector with gamma constants 
    static void setSigma(const std::string selection,const std::vector<double> gamma);
    static std::vector<double> sigma(const std::string selection) ;

    /// @brief set scale factor for extension used in minimization (debug use only)
    static void setExtscale(double es){ s_extscale=es; };
    static double extscale() { return s_extscale; };

    /// @brief set maximum size of source allowed in fit (in rad):
    static void setMaxSize(double ms){ s_maxsize=ms; };
    static double maxSize() { return s_maxsize; };

    std::string type()const {return m_type;};
    
    std::vector<double> sourceParameters() const { return m_sourceParameters; };
    std::vector<double> sourceParErrors() const   { return m_sourceParErrors; };
    double errorX() const {return m_errorX;};
    double errorY() const {return m_errorY;};
    std::map< std::string,std::vector<double> > errorsMINOS() const ;


    
    ///! Set the global diffuse background function, return current value 
    static  skymaps::SkySpectrum* set_diffuse( skymaps::SkySpectrum* diffuse, 
					       double exposure = 1.0);

    ///! Set exposure maps (for each event class, return current value 
    static  skymaps::SkySpectrum* set_exposure( skymaps::SkySpectrum* exposure, 
					       int event_class);
    
    ///! add a point source fit to the background for subsequent fits
    void addBackgroundPointSource(const SourceLikelihood* fit);
    
    ///! remove all such
    void clearBackgroundPointSource();

    ///! access to background model 
    const skymaps::SkySpectrum * background()const;

//     /// @brief set the integration tolerance for the background, return present value
//     static double set_tolerance(double tol);
    
    /// @brief special display function
    /// @param dir direction
    /// @param energy selects energy band
    /// @mode 0: same as the operator; 1: data; 2: background; 3:fit; 4:residual
    ///     
    double display(const astro::SkyDir& dir, double energy, int mode)const;
    std::vector<double> energyList() const;

  private:
    void setup(skymaps::BinnedPhotonData& data);
    std::string m_name;
    astro::SkyDir m_dir; ///< common direction
    double m_loglike;    ///< total loglike
    double m_TS;         ///< total TS value

    std::string m_type;
    unsigned int m_npar;
    std::vector<double> m_sourceParameters;
    std::vector<double> m_sourceParErrors;
    std::vector<bool>   m_sourceParametersFixMask;
    std::vector<double> m_errMINOSParabolic;
    std::vector<double> m_errMINOSPlus;
    std::vector<double> m_errMINOSMinus;
    std::vector<double> m_errMINOSGlobalCorr;
    
    double m_errorX;
    double m_errorY;
    
    std::ostream * m_out;
    std::ostream& out()const{return *m_out;}
    mutable CLHEP::Hep3Vector m_gradient; ///< current gradient
    
    static skymaps::SkySpectrum* s_diffuse; ///< global diffuse used by all SL objects
    static std::vector<skymaps::SkySpectrum*> s_exposure; ///< global exposure maps used by all SL objects
    
    skymaps::CompositeSkySpectrum * m_background;  ///< background spectrum to use
    
    static double s_emin, s_emax, s_minalpha, s_TSmin, s_tolerance, 
      s_maxstep,s_accuracy,s_maxsize,s_extscale; //
    static int s_useMinuit;
    static int s_useSimplex;
    static int s_itermax, s_verbose;
    static int s_useMinos;
    static int s_useGradient;
    static int s_npar;
    static std::vector<double> s_gamma_front; 
    static std::vector<double> s_sigma_front; 
    static std::vector<double> s_gamma_back; 
    static std::vector<double> s_sigma_back; 
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

