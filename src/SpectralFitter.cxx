#include "pointlike/SpectralFitter.h"
#include "pointlike/SourceLikelihood.h"
#include "pointlike/HypergeometricFunction.h"
#include "st_facilities/GaussianQuadrature.h"

#include "embed_python/Module.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <stdlib.h>
#include <time.h>

#include <TMinuit.h>

namespace pointlike{

  int gFitCounter=0;
  int gUseExposurePDF=1;
  int gUseMultiwavelengthData=0;

  // Global pointers
  pointlike::SourceLikelihood* gSourcePointer = NULL;
  pointlike::SpectralModel* gModelPointer = NULL;
  pointlike::MWData* gMWDataPointer = NULL;

  // Initialize static variables
  double SpectralFitter::s_lower_bound(200.);
  double SpectralFitter::s_upper_bound(6.e4);

  int    SpectralFitter::s_useDefaultParams(0);
  int    SpectralFitter::s_useSimplex(1);
  int    SpectralFitter::s_useGradient(0);
  int    SpectralFitter::s_useMinos(0);
  
  int SpectralFitter::s_useExposurePDF(1);

  int SpectralFitter::s_useMultiwavelengthData(0);

  double SpectralFitter::s_stepIncrement(0.01);

  double SpectralFitter::s_accuracy(0.0001);

  int    SpectralFitter::s_useUpperLimit(1);
  double SpectralFitter::s_TS_threshold(25.);
  double SpectralFitter::s_index(2.0);
  double SpectralFitter::s_upper_limit_lower_bound(100.);
  double SpectralFitter::s_upper_limit_upper_bound(3.e5);
  double SpectralFitter::s_band_confidence_limit(0.9);

  //-------------------------------------------------------------------------
  //  Define the wrapper function input to Minuit
  //-------------------------------------------------------------------------

  void minuit_spectral_wrapper(Int_t &,
			       Double_t* derivative,
			       Double_t & Loglike,
			       Double_t par[],
			       Int_t iflag){
    if(!gSourcePointer)
      throw std::invalid_argument("SpectralFitter::gSourcePointer not set");
    if(!gModelPointer)
      throw std::invalid_argument("SpectralFitter::gModelPointer not set");

    // Number of energy bins
    int numBins=gSourcePointer->size();

    double E_min,E_max;
    int n;
    double mu_alpha,sigma_alpha,mu_exp,sigma_exp;

    int npar=gModelPointer->get_npar();
    std::vector<double> specParams(npar);
    
    // Working in log space
    for(int i=0;i<npar;i++){
      specParams[i]=exp(par[i]);
    }
      
    gModelPointer->set_params(specParams);
    
    PDF pdf(gSourcePointer,gModelPointer,gUseExposurePDF);

    double logLikeSum=0.;
    
    for(int bin=0;bin<numBins;bin++){
      E_min=gSourcePointer->at(bin)->band().emin();
      E_max=gSourcePointer->at(bin)->band().emax();
      
      // Enforce energy range for fitting - bin must fully contained
      if(E_min < gModelPointer->get_lower_bound() || 
	 E_max > gModelPointer->get_upper_bound()) continue;

      pdf.set_bin(bin);
      logLikeSum-=pdf.get_loglikelihood(); // Minuit minimizes -> negative sign
    }

    // Incorporate multiwavelength data
    if(gUseMultiwavelengthData){
      
      if(!gMWDataPointer)
	 throw std::invalid_argument("SpectralFitter::gMWDataPointer not set");
      
      PDFMultiwavelength MWpdf(gModelPointer);

      double E,dNdE,dNdE_err_lo,dNdE_err_hi;

      int MWbins=gMWDataPointer->get_E().size();
      
      for(int MWbin=0;MWbin<MWbins;MWbin++){
	E           = gMWDataPointer->get_E()[MWbin];
	dNdE        = gMWDataPointer->get_dNdE()[MWbin];
	dNdE_err_lo = gMWDataPointer->get_dNdE_err_lo()[MWbin];
	dNdE_err_hi = gMWDataPointer->get_dNdE_err_hi()[MWbin];
	
	MWpdf.set_bin(E,dNdE,dNdE_err_lo,dNdE_err_hi);
	
	std::cout << -1.*MWpdf.get_loglikelihood() << std::endl;

	logLikeSum-=MWpdf.get_loglikelihood(); // Minuit minimizes -> negative sign
      }
    }

    Loglike=logLikeSum;

    std::cout << "Iteration " << gFitCounter << ":   ";
    gModelPointer->show();
    std::cout << "   -LogLikeSum = " 
	      << std::setprecision(4) << Loglike << std::endl;
    
    ++gFitCounter;
  }

  //--------------------------------------------------------------------------
  //  Constructor for spectral fitting
  //--------------------------------------------------------------------------

  SpectralFitter::SpectralFitter(SourceLikelihood& source,
				 SpectralModel& model){
    
    // Set the global pointers for access by the Minuit wrapper function
    gSourcePointer = &source;
    gModelPointer = &model;
    m_model_pointer = &model;

    // Set the energy range for spectral fitting
    this->setFitRange(s_lower_bound,s_upper_bound);

    m_npar=model.get_npar();

    this->setCombined();

    this->initialize();

  }

  //--------------------------------------------------------------------------
  //  Initialize Minuit parameters
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::initialize(){

    gFitCounter=0;

    m_Minuit=new TMinuit(m_npar);
    //m_Minuit.SetPrintLevel(verbose()-1);
    
    m_Minuit->SetFCN(minuit_spectral_wrapper);

    // The parameters will be optimized in log space
    m_parName  = gModelPointer->get_parName();
    m_stepSize = gModelPointer->get_stepSize(s_stepIncrement);
    m_minVal   = gModelPointer->get_minVal();
    m_maxVal   = gModelPointer->get_maxVal();

    // Use the default starting values or those provided by the model?
    if(s_useDefaultParams){
      m_par = gModelPointer->get_default_par();
      for(int i=0;i<m_npar;i++){
	m_par[i]=log(m_par[i]);
      }
    }
    else{
      m_par = gModelPointer->get_params();
      for(int i=0;i<m_npar;i++){
	m_par[i]=log(m_par[i]);
      }
    }

    for(int i=0;i<m_npar;i++){
      m_Minuit->DefineParameter(i,m_parName[i].c_str(),m_par[i],
			      m_stepSize[i],m_minVal[i],m_maxVal[i]);
      m_specParams.push_back(m_par[i]);
      m_specParamErrors.push_back(0.);
    }

  }

  //--------------------------------------------------------------------------
  //  Destructor
  //--------------------------------------------------------------------------

  SpectralFitter::~SpectralFitter(){
    if(m_Minuit) delete m_Minuit;
    if(m_FeldmanCousins) delete m_FeldmanCousins;
    if(m_model_pointer){
      m_model_pointer=NULL;
      delete m_model_pointer;
    }
    if(m_MWData) delete m_MWData;
  }

  //--------------------------------------------------------------------------
  //  Spectral fitting method
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::specfitMinuit(double scale){
    // Initialize timer
    time_t start,end;
    time(&start);

    // Check if using decorrelation energy
    if(gModelPointer->get_spec_type()=="POWER_LAW" && scale!=-1.){
      gModelPointer->set_scale(scale);
    }

    // Check if using the full exposure PDF convolution
    gUseExposurePDF = s_useExposurePDF;

    // Check if incorporating multiwavelength data
    gUseMultiwavelengthData = s_useMultiwavelengthData;
    if(s_useMultiwavelengthData)
      gMWDataPointer=m_MWData;

    // The Minuit output flag - query after each command
    Int_t ierflag=0;

    // Errors for likelihood fitting
    Double_t arglist[20];
    arglist[0]=0.5;
    int nargs=1;
    m_Minuit->mnexcm("SET ERR",arglist,nargs,ierflag);
    
    // Force Minuit to use the user-defined gradient function
    if(s_useGradient){ 
      nargs=1; arglist[0] = 1; 
      m_Minuit->mnexcm("SET GRA", arglist, nargs, ierflag);
    }

    // Set Minuit strategy
    nargs=1; arglist[0] = 1; 
    m_Minuit->mnexcm("SET STR", arglist, nargs, ierflag);
    
    // arglist[0] corresponds to the maximum number of function calls
    arglist[0] = 2000; 
    // arglist[1] controls fit tolerance
    if (s_useSimplex==1) arglist[1]=s_accuracy;
    else arglist[1] = s_accuracy*1000.;

    std::cout << std::endl << "===== BEGIN SPECTRAL FITTING ======" 
	      << std::endl << std::endl;
	   
    nargs = 2;
    if(s_useSimplex==1)
      m_Minuit->mnexcm("SIMPLEX", arglist, nargs, ierflag);
    else
      m_Minuit->mnexcm("MIGRAD", arglist, nargs, ierflag);  
    if((ierflag == 4 || s_useSimplex==1)&&(!s_useMinos)){
      m_Minuit->mnexcm("HESSE", arglist, nargs, ierflag);
    }
  
    // Get the spectral parameters from Minuit
    for(int i=0;i<m_npar;i++){
      m_Minuit->GetParameter(i,m_specParams[i],m_specParamErrors[i]);
    }
    
    std::cout << std::endl << "===== SPECTRAL FITTING RESULTS ====="
	      << std::endl;

    std::cout << std::endl << "Source: " << gSourcePointer->name()
	      << std::endl;

    std::cout << std::endl << "Spectral fitting in range "
	      << std::fixed
	      << std::setprecision(0)
	      << gModelPointer->get_lower_bound()
	      << " MeV < E < "
	      << gModelPointer->get_upper_bound()
	      << " MeV"
	      << std::endl;

    gModelPointer->show_legend();

    for(int i=0;i<m_npar;i++){
      m_specParams[i]=exp(m_specParams[i]); // Return from log space
      m_specParamErrors[i]=m_specParams[i]*(exp(m_specParamErrors[i])-1);
      std::cout << "PARAMETER: " 
		<< std::scientific
		<< std::setprecision(4)
		<< std::setw(15) << m_parName[i] << " = "
		<< m_specParams[i]
		<< " +/- " << m_specParamErrors[i]
		<< std::endl;
    }

    // Record the optimized model parameters and errors
    gModelPointer->set_params(m_specParams);
    gModelPointer->set_param_errors(m_specParamErrors);
    
    // Set the model dependent energies and exposures
    int numBins=gSourcePointer->size();

    double E_min,E_max;
    for(int bin=0; bin<numBins; bin++){
      E_min=gSourcePointer->at(bin)->band().emin();
      E_max=gSourcePointer->at(bin)->band().emax();
      
      // Enforce energy range for spectral fitting
      if(E_min < gModelPointer->get_lower_bound() || 
	 E_max > gModelPointer->get_upper_bound()) continue;

      else{
	m_model_energies.push_back(gModelPointer->get_model_E(E_min,E_max));
	m_exposures.push_back(gModelPointer->get_model_exposure(gSourcePointer,bin));
      }
    }
    
    gModelPointer->set_model_energies(m_model_energies);
    gModelPointer->set_exposures(m_exposures);

    double fmin,fedm,errdef;
    int npari,nparx,istat;
    m_Minuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
    gModelPointer->set_logLikeSum(fmin);

    // Pick up the covariance matrix entries
    double C[m_npar][m_npar];
    m_Minuit->mnemat(&C[0][0],m_npar);
    for(int i=0;i<m_npar;i++){
      for(int j=0;j<m_npar;j++){

	// Return from log space
	//m_covar_entries.push_back(m_specParams[i]*m_specParams[j]*(exp(C[i][j])-1));

	m_covar_entries.push_back(C[i][j]);
      }
    }
    gModelPointer->set_covar_entries(m_covar_entries);

    std::vector<double> flux_errors=gModelPointer->get_flux_errors(100.,3.e5);

    std::cout << std::endl 
	      << "Average integrated photon flux (100 MeV < E < 300 GeV) = " 
	      << gModelPointer->get_flux(100.,3.e5)
	      << " +/- "
	      << flux_errors[0]
      	      << "/"
	      << flux_errors[1]
	      << " ph cm^-2 s^-1"
	      << std::endl;

    std::cout << std::endl
	      << "Minimized -logLikeSum = " << fmin << std::endl
	      << std::endl;

    if(ierflag==4){
      std::cerr<<"WARNING: Minuit returned ierflag=4: Fit did not converge."<<std::endl;
    }

    time(&end);
    double SFtime=difftime(end,start);
    std::cout << std::fixed << std::setprecision(0)
	      << "Spectral fitting computation time = " 
	      << SFtime 
	      << " seconds" 
	      << std::endl;

    // Calculate pivot energy and repeat the fit
    if(gModelPointer->get_spec_type()=="POWER_LAW" && scale==-1.){
      double log_covar[2][2];
      m_Minuit->mnemat(&log_covar[0][0],2);
      double log_covar01=log_covar[0][1];
      double covar01=m_specParams[0]*m_specParams[1]*(exp(log_covar01)-1);
      double log_covar11=log_covar[1][1];
      double covar11=pow(m_specParams[1],2)*(exp(log_covar11)-1);
      double E_decorrelation=100.*exp(covar01/(m_specParams[0]*covar11));
      std::cout << "Calculate pivot energy using covariance matrix"
		<< std::endl << std::endl
		<< "Pivot energy = " 
		<< std::fixed << std::setprecision(1)
		<< E_decorrelation << " MeV"
		<< std::endl << std::endl
		<< "Repeat spectral fitting" 
		<< std::endl << std::endl;
      this->initialize();
      this->specfitMinuit(E_decorrelation);
    }
    else{
      if(s_useUpperLimit)
	this->setFluxUpperLimits();
    }

  }

  //--------------------------------------------------------------------------
  //  Integral flux upper limit calculation 
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::setFluxUpperLimits(std::vector<double> confidence_limits){
    // Initialize timer
    time_t start,end;
    time(&start);

    if(!gSourcePointer)
      throw std::invalid_argument("SpectralFitter::gSourcePointer not set");

    if(!confidence_limits.empty()) m_confidence_limits=confidence_limits;
    else if(m_confidence_limits.empty()) this->setConfidenceLimits();

    std::cout << std::endl
	      << "===== INTEGRAL PHOTON FLUX UPPER LIMITS ====="
	      << std::endl;

    // Choose fitted spectral model or default power law model
    if(gSourcePointer->TS()>s_TS_threshold && gModelPointer && gFitCounter>0){
      m_model_pointer=gModelPointer;
      std::cout << std::endl 
		<< "Using fitted spectral model to compute upper limits"
		<< std::endl;
    }
    else{
      std::vector<double> pl_params(2);
      pl_params[0]=1.e-9;
      pl_params[1]=s_index;
      m_model_pointer=new PowerLaw(pl_params);
      std::cout << std::endl
		<< "Assume a power law spectral model with photon index "
		<< std::fixed << std::setprecision(2)
		<< s_index
		<< std::endl;
    }
    
    // Number of energy bins
    int numBins=gSourcePointer->size();

    double E_min,E_max;

    int n;
    double alpha, alpha_error; 
    double band_exposure;

    double psf_correction;
    double exposure_sum=0;
    double weighted_exposure_sum=0;

    double fit_lower_bound=s_upper_limit_upper_bound;
    double fit_upper_bound=s_upper_limit_lower_bound;

    double N_total=0;
    double N_signal=0;
    double N_background=0;

    m_FeldmanCousins=new FeldmanCousins();

    double FC_upper_limit;

    m_band_upper_limits.clear();

    std::cout << std::endl
	      << std::fixed << std::setprecision(2)
	      << "Integral photon fluxes within energy bands [ ph cm^-2 s^-1 ] ( Energies in MeV ; Confidence Limit = " << s_band_confidence_limit << ")"
	      << std::endl << std::endl
	      << std::setw(10) << "emin" << std::setw(10) << "emax" << std::setw(10) << "evclass" << std::setw(12) << "flux"
	      << std::endl;

    // Total photon counts
    for(int bin=0;bin<numBins;bin++){
      
      E_min=gSourcePointer->at(bin)->band().emin();
      E_max=gSourcePointer->at(bin)->band().emax();
      
      // Enforce range of upper limit calculation
      if(E_min < s_upper_limit_lower_bound || 
	 E_max > s_upper_limit_upper_bound) continue;

      if(E_min < fit_lower_bound) fit_lower_bound=E_min;
      if(E_max > fit_upper_bound) fit_upper_bound=E_max;

      n=gSourcePointer->at(bin)->photons();
      alpha=gSourcePointer->at(bin)->alpha();
      alpha_error=gSourcePointer->at(bin)->sigma_alpha();

      N_total=N_total+n;
      N_signal=N_signal+alpha*n;

      psf_correction=gSourcePointer->at(bin)->psfIntegral();

      // Model dependent exposure calculation
      band_exposure=m_model_pointer->get_model_exposure(gSourcePointer,bin);

      exposure_sum+=band_exposure;
      weighted_exposure_sum+=psf_correction*band_exposure;
      
      // Set flux upper limits for individual bands including signal fraction error
      FC_upper_limit=m_FeldmanCousins->get_upper_limit(gSourcePointer,bin,s_band_confidence_limit);

      m_band_upper_limits.push_back(FC_upper_limit/psf_correction);
      m_energy_upper_limits.push_back(m_model_pointer->get_model_E(E_min,E_max));
      m_exposure_upper_limits.push_back(band_exposure);

      std::cout << std::fixed << std::setprecision(0)
		<< std::setw(10) << E_min
		<< std::setw(10) << E_max
		<< std::setw(10) << gSourcePointer->at(bin)->band().event_class()
		<< std::scientific << std::setprecision(2) 
		<< std::setw(12) << FC_upper_limit/(psf_correction*band_exposure)
		<< std::endl;
    }

    // Get the range of energy used for fitting and use it for calculating the flux ratio
    double flux_ratio=m_model_pointer->get_flux(s_upper_limit_lower_bound,s_upper_limit_upper_bound)/m_model_pointer->get_flux(fit_lower_bound,fit_upper_bound);

    // Integrated exposure for combined front and back
    double integrated_exposure;
    if(m_combined)
      integrated_exposure=m_model_pointer->get_integrated_exposure(gSourcePointer,s_upper_limit_lower_bound,s_upper_limit_upper_bound);
    else
      integrated_exposure=m_model_pointer->get_integrated_exposure(gSourcePointer,s_upper_limit_lower_bound,s_upper_limit_upper_bound,0)+
	m_model_pointer->get_integrated_exposure(gSourcePointer,s_upper_limit_lower_bound,s_upper_limit_upper_bound,1);
    
    // Fraction of photons contained within region of interest
    double containment_fraction=weighted_exposure_sum/exposure_sum;

    N_background=N_total-N_signal;

    std::cout << std::endl
	      << std::fixed << std::setprecision(0)
	      << "Total photons  = "
	      << N_total
	      << std::endl
	      << "Signal photons = "
	      << N_signal
	      << std::endl << std::endl
	      << std::scientific << std::setprecision(4)
	      << "Integrated exposure = "
	      << integrated_exposure
	      << " cm^2 s"
	      << std::endl << std::endl << std::fixed << std::setprecision(0)
	      << "Confidence limit       Flux upper limit ( " << s_upper_limit_lower_bound << " MeV < E < " << s_upper_limit_upper_bound << " MeV ) [ ph cm^-2 s^-1 ]"
	      << std::endl;

    m_FeldmanCousins->set_counts(N_total,N_background);

    m_flux_upper_limits.clear();

    for(int i=0;i<m_confidence_limits.size();i++){

      FC_upper_limit=m_FeldmanCousins->get_upper_limit(m_confidence_limits[i]);

      m_flux_upper_limits.push_back((flux_ratio*FC_upper_limit)/(integrated_exposure*containment_fraction));

      std::cout << std::endl
		<< std::setw(23) << std::left
		<< std::fixed
		<< std::setprecision(2)
		<< m_confidence_limits[i]
		<< std::setw(15) << std::left
		<< std::scientific
		<< m_flux_upper_limits[i]
		<< std::endl;
    }

    // Restore the original spectral model
    m_model_pointer=gModelPointer;

    time(&end);
    double ULtime=difftime(end,start);
    std::cout << std::fixed << std::setprecision(0)
	      << "Upper limit computation time = " 
	      << ULtime 
	      << " seconds" 
	      << std::endl;

  }

  //--------------------------------------------------------------------------
  //  Functions to set spectral fitting parameters
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::setFitRange(double lower_bound,double upper_bound){
        gModelPointer->set_energy_range(lower_bound,upper_bound);
  }

  void pointlike::SpectralFitter::setFluxUpperLimitRange(double lower_bound,double upper_bound){
    s_upper_limit_lower_bound=lower_bound;
    s_upper_limit_upper_bound=upper_bound;
  }

  void pointlike::SpectralFitter::setConfidenceLimits(std::vector<double> confidence_limits){
    if(confidence_limits.size()==0){
      std::vector<double> cl(4);
      cl[0]=0.99;
      cl[1]=0.95;
      cl[2]=0.90;
      cl[3]=0.68;
      m_confidence_limits=cl;
    }
    else m_confidence_limits=confidence_limits;
  }

  //--------------------------------------------------------------------------
  //  Function for using combined front and back energy binning
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::setCombined(){
    if(!gSourcePointer)
      throw std::invalid_argument("SpectralFitter::gSourcePointer not set");
    
    m_combined=1;

    int numBins=gSourcePointer->size();

    for(int bin=0;bin<numBins;bin++){
      if(gSourcePointer->at(bin)->band().event_class()==0)
	m_combined=0;
    }
      
    if(!m_combined)
      std::cout << std::endl 
		<< "Treating front and back conversion events separately" 
		<< std::endl << std::endl;
    else
      std::cout << std::endl 
		<< "Using combined front and back energy bins" 
		<< std::endl << std::endl;
  }

}
