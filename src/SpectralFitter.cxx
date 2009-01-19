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

#include <TMinuit.h>

namespace pointlike{

  int gFitCounter=0;

  // Global pointers
  pointlike::SourceLikelihood* gSourcePointer = NULL;
  pointlike::SpectralModel* gModelPointer = NULL;
  
  // Initialize static variables
  double SpectralFitter::s_lower_bound(200.);
  double SpectralFitter::s_upper_bound(3.e4);

  int    SpectralFitter::s_useDefaultParams(0);
  int    SpectralFitter::s_useSimplex(1);
  int    SpectralFitter::s_useGradient(0);
  int    SpectralFitter::s_useMinos(0);
  
  double SpectralFitter::s_stepIncrement(0.01);

  double SpectralFitter::s_accuracy(0.0001);

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
    
    PDF pdf(gModelPointer);

    double exposure_error_frac;

    double logLikeSum=0.;
    
    for(int bin=0;bin<numBins;bin++){
      	
      n=gSourcePointer->at(bin)->photons();

      if(n==0) continue; // Skip empty energy bins

      E_min=gSourcePointer->at(bin)->band().emin();
      E_max=gSourcePointer->at(bin)->band().emax();
      
      // Enforce energy range for fitting - bin must fully contained
      if(E_min < gModelPointer->get_lower_bound() || 
	 E_max > gModelPointer->get_upper_bound()) continue;
      

      mu_alpha=gSourcePointer->at(bin)->alpha();
      sigma_alpha=gSourcePointer->at(bin)->sigma_alpha();

      // Model dependent exposure calculation
      mu_exp=gModelPointer->get_model_exposure(gSourcePointer,bin);
      exposure_error_frac=gModelPointer->get_exposure_uncertainty(E_min,E_max);
      sigma_exp=exposure_error_frac*mu_exp; 
      
      pdf.set_bin(n,E_min,E_max,mu_exp,sigma_exp,mu_alpha,sigma_alpha);

      logLikeSum-=pdf.get_loglikelihood(); // Negative sign -> Minuit minimizes
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
    
    // Set the energy range for spectral fitting
    this->setFitRange(s_lower_bound,s_upper_bound);

    m_npar=model.get_npar();

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
  }

  //--------------------------------------------------------------------------
  //  Spectral fitting method
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::specfitMinuit(double scale){

    // Check if using decorrelation energy
    if(gModelPointer->get_spec_type()=="POWER_LAW" && scale!=-1.){
      gModelPointer->set_scale(scale);
    }

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

    if(ierflag==4){
      std::cerr<<"WARNING: Minuit returned ierflag=4: Fit did not converge."<<std::endl;
    }

  }

  //--------------------------------------------------------------------------
  //  Integral flux upper limit calculation 
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::getFluxUpperLimit(std::vector<double> confidence_limit){
    
    if(!gSourcePointer)
      throw std::invalid_argument("SpectralFitter::gSourcePointer not set");

    if(confidence_limit.size()==0)
      confidence_limit.push_back(0.9);

    double index=2.;

    std::vector<double> pl_params(2);
    pl_params[0]=1.e-9;
    pl_params[1]=index;
    m_pl_model=new PowerLaw(pl_params);

    // Number of energy bins
    int numBins=gSourcePointer->size();

    double E_min,E_max;

    // Get the lowest energy used for fitting
    double fit_lower_bound=gSourcePointer->at(0)->band().emin();
    double flux_ratio=m_pl_model->get_flux(100.,3.e5)/m_pl_model->get_flux(fit_lower_bound,3.e5);

    // Fraction of photons contained
    double containment_fraction=0.68;

    int n;
    double alpha,exposure;

    double N_total=0;
    double N_signal=0;
    double N_background=0;

    double integrated_exposure=0;

    for(int bin=0;bin<numBins;bin++){
      
      // Photons counts 
      n=gSourcePointer->at(bin)->photons();
      alpha=gSourcePointer->at(bin)->alpha();

      N_total=N_total+n;
      N_signal=N_signal+floor(alpha*n);

      E_min=gSourcePointer->at(bin)->band().emin();
      E_max=gSourcePointer->at(bin)->band().emax();
      
      // Model dependent exposure calculation
      exposure=m_pl_model->get_model_exposure(gSourcePointer,bin);
      integrated_exposure=integrated_exposure+exposure;  

    }

    integrated_exposure=m_pl_model->get_integrated_exposure(gSourcePointer,100.,3.e5);

    std::cout << std::endl
	      << "===== INTEGRAL PHOTON FLUX UPPER LIMITS ====="
	      << std::endl;

    std::cout << std::endl
	      << "Assume a power law spectral model with spectral index "
	      << std::fixed << std::setprecision(2)
	      << index
	      << std::endl;

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
	      << std::endl << std::endl
	      << "Confidence limit       Flux upper limit ( 100 MeV < E < 300 GeV ) [ ph cm^-2 s^-1 ]"
	      << std::endl;

    m_FeldmanCousins=new FeldmanCousins(N_signal,N_background);
    
    double FC_upper_limit,flux_upper_limit;

    for(int i=0;i<confidence_limit.size();i++){

      FC_upper_limit=m_FeldmanCousins->get_upper_limit(confidence_limit[i]);

      flux_upper_limit=(flux_ratio*FC_upper_limit)/(integrated_exposure*containment_fraction);

      std::cout << std::endl
		<< std::setw(23) << std::left
		<< std::fixed
		<< std::setprecision(2)
		<< confidence_limit[i]
		<< std::setw(15) << std::left
		<< std::scientific
		<< flux_upper_limit
		<< std::endl;
    }

  }

  //--------------------------------------------------------------------------
  //  Functions to set spectral fitting parameters
  //--------------------------------------------------------------------------

  void pointlike::SpectralFitter::setFitRange(double lower_bound,double upper_bound){
        gModelPointer->set_energy_range(lower_bound,upper_bound);
  }

}
