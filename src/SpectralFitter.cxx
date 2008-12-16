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
  double SpectralFitter::s_emin(200.);
  double SpectralFitter::s_emax(3.e5);

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

    double logLikeSum=0.;
    
    for(int bin=0;bin<numBins;bin++){
      	
      n=gSourcePointer->at(bin)->photons();
      
      if(n==0) continue; // Skip empty energy bins

      E_min=gSourcePointer->at(bin)->band().emin();
      E_max=gSourcePointer->at(bin)->band().emax();

      mu_alpha=gSourcePointer->at(bin)->alpha();
      sigma_alpha=gSourcePointer->at(bin)->sigma_alpha();
      
      mu_exp=gSourcePointer->at(bin)->exposure();
      sigma_exp=0.01*mu_exp; // Until we get a better estimate from C&A group
      
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

    m_npar=model.get_npar();
    
    // Set the energy range for spectral fitting
    // Do we want a user specified energy range for spectral fitting?
    s_emin=source.emin();
    s_emax=source.emax();

    m_Minuit=new TMinuit(m_npar);
    //m_Minuit.SetPrintLevel(verbose()-1);
    
    m_Minuit->SetFCN(minuit_spectral_wrapper);
  
    // The parameters will be optimized in log space
    m_parName  = model.get_parName();
    m_stepSize = model.get_stepSize(s_stepIncrement);
    m_minVal   = model.get_minVal();
    m_maxVal   = model.get_maxVal();
    
    // Use the default starting values or those provided by the model
    if(s_useDefaultParams){
      m_par = model.get_default_par();
      for(int i=0;i<m_npar;i++){
	m_par[i]=log(m_par[i]);
      }
    }
    else{
      m_par = model.get_params();
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

  void pointlike::SpectralFitter::specfitMinuit(){

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

    gModelPointer->show_legend();

    for(int i=0;i<m_npar;i++){
      m_specParams[i]=exp(m_specParams[i]); // Return from log space
      m_specParamErrors[i]=m_specParams[i]*(exp(m_specParamErrors[i])-1);
      std::cout << "PARAMETER: " 
		<< std::scientific
		<< std::setw(15) << m_parName[i] << " = "
		<< m_specParams[i]
		<< " +/- " << m_specParamErrors[i]
		<< std::endl;
    }

    // Record the optimized model parameters and errors
    gModelPointer->set_params(m_specParams);
    gModelPointer->set_param_errors(m_specParamErrors);

    double fmin,fedm,errdef;
    int npari,nparx,istat;
    m_Minuit->mnstat(fmin,fedm,errdef,npari,nparx,istat);
    gModelPointer->set_logLikeSum(fmin);

    std::cout << std::endl 
	      << "Integrated photon flux (300 GeV > E > 100 MeV) = " 
	      << gModelPointer->get_flux(100.,3.e5) 
	      << " ph/cm^2/s"
	      << std::endl;

    std::cout << std::endl
	      << "Minimized -logLikeSum = " << fmin << std::endl
	      << std::endl;

    if(ierflag==4){
      std::cerr<<"WARNING: Minuit returned ierflag=4: Fit did not converge."<<std::endl;
    }

  }

}
