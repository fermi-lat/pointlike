#ifndef tools_SpectralFitter_h
#define tools_SpectralFitter_h

#include "pointlike/HypergeometricFunction.h"
#include "st_facilities/GaussianQuadrature.h"
#include "pointlike/SourceLikelihood.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <stdlib.h>

#include <cstdlib>
#include <time.h>

#include <TMinuit.h>

namespace embed_python{ class Module; }

namespace pointlike{

  class SourceLikelihood;

  //---------------------------------------------------------------------------
  //  Gaussian random number generation centered on zero with a standard
  //  deviation of one (Box-Muller transformation).
  //---------------------------------------------------------------------------

  class Random{
  public:
    Random(){
      srand(time(NULL));
    }

    double rand_unif(){
      return (rand()%10000)/10000.;
    }

    double rand_gauss(){
      double x1, x2, w, y1, y2;

      do {
	x1 = 2.0 * this->rand_unif() - 1.0;
	x2 = 2.0 * this->rand_unif() - 1.0;
	w = x1 * x1 + x2 * x2;
      } while ( w >= 1.0 );

      w = sqrt( (-2.0 * log( w ) ) / w );
      y1 = x1 * w;
      y2 = x2 * w;

      return y1;
    }
  };

  //---------------------------------------------------------------------------
  //  Base class for spectral models
  //---------------------------------------------------------------------------

  class SpectralModel{
  protected:
    std::string m_spec_type;
    double m_logLikeSum;
    std::vector<double> m_model_energies;
    std::vector<double> m_exposures;
    
    double m_scale;

    double m_lower_bound,m_upper_bound;

    std::vector<double> m_covar_entries;

  public:
    SpectralModel(std::string spec_type=""){
      m_spec_type=spec_type;
    }

    std::string get_spec_type(){
      return m_spec_type;
    }
    
    void set_logLikeSum(double logLikeSum){
      m_logLikeSum=logLikeSum; 
    }

    double get_logLikeSum(){
      return m_logLikeSum;
    }

    //virtual ~SpectralModel()=0;

    virtual void set_params(const std::vector<double> params)=0;
    virtual std::vector<double> get_params()=0;

    virtual void set_param_errors(const std::vector<double> param_errors)=0;
    virtual std::vector<double> get_param_errors()=0;
    
    virtual double get_flux(double E_min,double E_max)=0;

    std::vector<double> get_flux_errors(double E_min,double E_max){
      std::vector<double> best_fit_params=this->get_params();
      std::vector<double> param_errors=this->get_param_errors();
      double best_fit_flux=this->get_flux(E_min,E_max);

      int npar=this->get_npar();
      std::vector<double> new_params(npar);
            
      Random generator;

      std::vector<double> new_fluxes;
      
      double new_flux=1;

      int maxit=400;
      for(int i=0;i<maxit;i++){
	do{
	  for(int j=0;j<npar;j++){
	    new_params[j]=best_fit_params[j]+generator.rand_gauss()*param_errors[j];
	  }
	  
	  this->set_params(new_params);
	  new_flux=this->get_flux(100.,3.e5);
	} while( new_flux <= 0.0 );
	
	new_fluxes.push_back(this->get_flux(100.,3.e5));
      }      
      
      std::sort(new_fluxes.begin(),new_fluxes.end());

      // Reset original parameter values
      this->set_params(best_fit_params);

      std::vector<double> flux_errors(2);
      flux_errors[0]=new_fluxes[ceil(0.84*new_fluxes.size())]-best_fit_flux;
      flux_errors[1]=best_fit_flux-new_fluxes[floor(0.16*new_fluxes.size())];	

      return flux_errors;
    }

    virtual double get_dNdE(double E)=0;

    double get_model_exposure(pointlike::SourceLikelihood* SourcePointer,int bin){
      double E_min=SourcePointer->at(bin)->band().emin();
      double E_max=SourcePointer->at(bin)->band().emax();

      double weighted_exposure_sum=0;
      double normalization=0;
      
      double frac,mid_frac,next_frac;
      double E,mid_E,next_E;
      double dNdE,dE;
      double exposure;
      
      double max=100.;

      for(double i=0.;i<max;i+=1.){
	frac=i/max;
	mid_frac=(i+0.5)/max;
	next_frac=(i+1.)/max;

	E=exp((1.-frac)*log(E_min)+frac*log(E_max));
	mid_E=exp((1.-mid_frac)*log(E_min)+mid_frac*log(E_max));
	next_E=exp((1.-next_frac)*log(E_min)+next_frac*log(E_max));

	dE=next_E-E;

	dNdE=this->get_dNdE(mid_E);
	exposure=SourcePointer->at(bin)->exposure(mid_E);
	weighted_exposure_sum+=dNdE*exposure*dE;
	normalization+=dNdE*dE;
      }
      
      return weighted_exposure_sum/normalization;
    }

    double get_integrated_exposure(pointlike::SourceLikelihood* SourcePointer,double E_min,double E_max,int evclass=-1){
      double weighted_exposure_sum=0;
      double normalization=0;
      
      double frac,mid_frac,next_frac;
      double E,mid_E,next_E;
      double dNdE,dE;
      double exposure;
      
      double max=1000.;

      for(double i=0.;i<max;i+=1.){
	
	frac=i/max;
	mid_frac=(i+0.5)/max;
	next_frac=(i+1.)/max;
	
	E=exp((1.-frac)*log(E_min)+frac*log(E_max));
	mid_E=exp((1.-mid_frac)*log(E_min)+mid_frac*log(E_max));
	next_E=exp((1.-next_frac)*log(E_min)+next_frac*log(E_max));
	
	dE=next_E-E;

	dNdE=this->get_dNdE(mid_E);

	// We may select any band for the purpose of calculation
	if(evclass==0 || evclass==1)
	  exposure=SourcePointer->at(0)->exposure(mid_E,evclass);
	else
	  exposure=SourcePointer->at(0)->exposure(mid_E);
	
	weighted_exposure_sum+=dNdE*exposure*dE;
	normalization+=dNdE*dE;
      }
      
      return weighted_exposure_sum/normalization;
    }

    double get_model_E(double E_min,double E_max){
      double weighted_E_sum=0;
      double normalization=0;

      double frac,mid_frac,next_frac;
      double E,mid_E,next_E;
      double dNdE,dE;

      double max=100.;
      for(double i=0.;i<max;i+=1.){
	frac=i/max;
	mid_frac=(i+0.5)/max;
	next_frac=(i+1.)/max;

	E=exp((1.-frac)*log(E_min)+frac*log(E_max));
	mid_E=exp((1.-mid_frac)*log(E_min)+mid_frac*log(E_max));
	next_E=exp((1.-next_frac)*log(E_min)+next_frac*log(E_max));

	dE=next_E-E;

	dNdE=this->get_dNdE(mid_E);

	weighted_E_sum+=dNdE*mid_E*dE;
	normalization+=dNdE*dE;
      }
      
      return weighted_E_sum/normalization;
    }

    double get_exposure_uncertainty(double E_min,double E_max){
      // Until we get a better estimate from the C&A group
      // See Vela I paper
      double E_log_center=exp((log(E_max)+log(E_min))/2.);
      if(E_log_center<1.e2) return 0.2;
      if(E_log_center<1.e4) return 0.1;
      else return 0.3;
    }

    virtual int get_npar()=0;

    virtual void show()=0;
  
    virtual void show_legend()=0;

    virtual std::vector<std::string> get_parName()=0;
    virtual std::vector<double> get_default_par()=0;
    virtual std::vector<double> get_stepSize(double increment)=0;
    virtual std::vector<double> get_minVal()=0;
    virtual std::vector<double> get_maxVal()=0;

    void set_model_energies(std::vector<double> model_energies){ m_model_energies=model_energies; }
    std::vector<double> get_model_energies(){ return m_model_energies; }
    
    void set_exposures(std::vector<double> exposures){ m_exposures=exposures; }
    std::vector<double> get_exposures(){ return m_exposures; } 
    
    void set_energy_range(double lower_bound,double upper_bound){
      m_lower_bound=lower_bound;
      m_upper_bound=upper_bound;
    }
    double get_lower_bound(){ return m_lower_bound; }
    double get_upper_bound(){ return m_upper_bound; }

    double set_scale(double scale){ m_scale=scale; }
    double get_scale(){ return m_scale; }

    void set_covar_entries(std::vector<double> covar_entries){ m_covar_entries=covar_entries; }
    std::vector<double> get_covar_entries(){ return m_covar_entries; }

  };

  //---------------------------------------------------------------------------
  //  Power law spectral model
  //
  //  dN/dE = N_0*(E/E_0)^(-gamma)
  //
  //  N_0   = prefactor        [ ph/cm^2/s/MeV ]
  //  E_0   = energy scale     [ default = 100 MeV ]
  //  gamma = photon index
  //---------------------------------------------------------------------------

  class PowerLaw: public SpectralModel{
    double m_prefactor,m_index;
    double m_prefactor_error,m_index_error;

    double function(double E){
      double result=m_prefactor*E*pow(E/m_scale,-1.*m_index)/(1.-m_index);
      return result; 
    }
  public:
    PowerLaw(std::vector<double> params):
      SpectralModel("POWER_LAW"),
      m_prefactor(params[0]),
      m_index(params[1]){
      this->set_scale(100.);
    }

    PowerLaw():
      SpectralModel("POWER_LAW"){
      this->set_params(this->get_default_par());
      this->set_scale(100.);
    }

    void set_params(const std::vector<double> params){
      m_prefactor=params[0];
      m_index=params[1];
    }
    
    std::vector<double> get_params(){
      std::vector<double> params(2);
      params[0]=m_prefactor;
      params[1]=m_index;
      return params;
    }

    void set_param_errors(const std::vector<double> param_errors){
      m_prefactor_error=param_errors[0];
      m_index_error=param_errors[1];
    }

    std::vector<double> get_param_errors(){
      std::vector<double> param_errors(2);
      param_errors[0]=m_prefactor_error;
      param_errors[1]=m_index_error;
      return param_errors;
    }

    double get_flux(double E_min,double E_max){
      double flux=function(E_max)-function(E_min);
      return flux;
    }

    virtual double get_dNdE(double E){
      return m_prefactor*pow(E/m_scale,-1*m_index);
    }

    int get_npar(){
      return 2;
    }

    virtual ~PowerLaw(){};

    virtual void show(){
      std::cout << std::setprecision(4)
		<< std::scientific 
	//<< "N_0 = " << m_prefactor
		<< "flux_100 = " << this->get_flux(100.,3.e5) 
		<< std::fixed 
		<< "   gamma =  " << m_index;
    }

    virtual void show_legend(){
      std::cout << std::endl
		<< "  Power law spectral model" << std::endl
		<< std::endl
		<< "  dN/dE = N_0*(E/E_0)^(-gamma)" << std::endl
		<< std::endl
		<< "  N_0   = prefactor      [ ph/cm^2/s/MeV ]" << std::endl
		<< "  E_0   = pivot energy   [ "
		<< std::fixed
		<< std::setprecision(1)
		<< this->get_scale() 
		<< " MeV ]" << std::endl
		<< "  gamma = photon index" << std::endl
		<< std::endl;
    }

    virtual std::vector<std::string> get_parName(){
      std::vector<std::string> parName(2);
      parName[0]="N_0";
      parName[1]="gamma";
      return parName;
    }

    virtual std::vector<double> get_default_par(){
      std::vector<double> par(2);
      par[0]=1.e-9;
      par[1]=2.1;
      return par;
    }

    virtual std::vector<double> get_stepSize(double increment){
      std::vector<double> stepSize(2);
      std::vector<double> minVal=this->get_minVal();
      std::vector<double> maxVal=this->get_maxVal();
      stepSize[0]=increment*(maxVal[0]-minVal[0]);
      stepSize[1]=increment*(maxVal[1]-minVal[1]);
      return stepSize;
    }

    virtual std::vector<double> get_minVal(){
      std::vector<double> minVal(2);
      minVal[0]=log(1.e-20);
      minVal[1]=log(0.1);
      return minVal;
    }
  
    virtual std::vector<double> get_maxVal(){
      std::vector<double> maxVal(2);
      maxVal[0]=log(1.e2);
      maxVal[1]=log(5.);
      return maxVal;
    }
  };

  //---------------------------------------------------------------------------
  //  Broken power law spectral model
  //
  //  dN/dE = N_0*(E/E_b)^(-gamma_1) , E < E_b
  //        = N_0*(E/E_b)^(-gamma_2) , E > E_b
  //
  //  N_0     = prefactor                              [ ph/cm^2/s/MeV ]
  //  gamma_1 = photon index in lower energy regime
  //  gamma_2 = photon index in higher energy regime
  //  E_b     = break energy                           [ MeV ]
  //---------------------------------------------------------------------------

  class BrokenPowerLaw: public SpectralModel{
    double m_prefactor,m_index_1,m_index_2,m_E_break;
    double m_prefactor_error,m_index_1_error,m_index_2_error,m_E_break_error;

    double function(double E){
      double result=0.;
      if(E<=m_E_break){
	result=m_prefactor*E*pow(E/m_E_break,-1.*m_index_1)/(1.-m_index_1);
      }else{
	result=m_prefactor*E*pow(E/m_E_break,-1.*m_index_2)/(1.-m_index_2);
      }
      return result; 
    }

  public:
    BrokenPowerLaw(std::vector<double> params):
      SpectralModel("BROKEN_POWER_LAW"),
      m_prefactor(params[0]),
      m_index_1(params[1]),
      m_index_2(params[2]),
      m_E_break(params[3]){
    }

    BrokenPowerLaw():
      SpectralModel("BROKEN_POWER_LAW"){
      this->set_params(this->get_default_par());
    }

    void set_params(const std::vector<double> params){
      m_prefactor=params[0];
      m_index_1=params[1];
      m_index_2=params[2];
      m_E_break=params[3];
    }

    std::vector<double> get_params(){
      std::vector<double> params(4);
      params[0]=m_prefactor;
      params[1]=m_index_1;
      params[2]=m_index_2;
      params[3]=m_E_break;
      return params;
    }

    void set_param_errors(const std::vector<double> param_errors){
      m_prefactor_error=param_errors[0];
      m_index_1_error=param_errors[1];
      m_index_2_error=param_errors[2];
      m_E_break_error=param_errors[3];
    }

    std::vector<double> get_param_errors(){
      std::vector<double> param_errors(4);
      param_errors[0]=m_prefactor_error;
      param_errors[1]=m_index_1_error;
      param_errors[2]=m_index_2_error;
      param_errors[3]=m_E_break_error;
      return param_errors;
    }

    double get_flux(double E_min,double E_max){
      double flux=0.;
      double epsilon=1e-10;
      if(E_max<=m_E_break){
	flux=function(E_max-epsilon)-function(E_min);
      }
      else if(E_min>=m_E_break){
	flux=function(E_max)-function(E_min+epsilon);
      }
      else{
	flux=(function(m_E_break-epsilon)-function(E_min))+(function(E_max)-function(m_E_break+epsilon));
      }
      return flux;
    }

    virtual double get_dNdE(double E){
      if(E<m_E_break) 
	return m_prefactor*pow(E/m_E_break,-1*m_index_1);
      else 
	return m_prefactor*pow(E/m_E_break,-1*m_index_2);
    }

    int get_npar(){
      return 4;
    }

    virtual ~BrokenPowerLaw(){};

    virtual void show(){
      std::cout << std::setprecision(4)
		<< std::scientific 
	//<< "N_0 = " << m_prefactor
		<< "flux_100 = " << this->get_flux(100.,3.e5)
		<< std::fixed
		<< "   gamma_1 = " << m_index_1
		<< "   gamma_2 = " << m_index_2
		<< "   E_b = " << m_E_break;
    }

    virtual void show_legend(){
      std::cout << std::endl 
		<< "  Broken power law spectral model" << std::endl
		<< std::endl
		<< "  dN/dE = N_0*(E/E_b)^(-gamma_1) , E < E_b" << std::endl
		<< "        = N_0*(E/E_b)^(-gamma_2) , E > E_b" << std::endl
		<< std::endl
		<< "  N_0     = prefactor                              [ ph/cm^2/s/MeV ]" << std::endl
		<< "  gamma_1 = photon index in lower energy regime" << std::endl
		<< "  gamma_2 = photon index in higher energy regime" << std::endl
		<< "  E_b     = break energy                           [ MeV ]" << std::endl
		<< std::endl;
    }

    virtual std::vector<std::string> get_parName(){
      std::vector<std::string> parName(4);
      parName[0]="N_0";
      parName[1]="gamma_1";
      parName[2]="gamma_2";
      parName[3]="E_b";
      return parName;
    }

    virtual std::vector<double> get_default_par(){
      std::vector<double> par(4);
      par[0]=1.e-9;
      par[1]=2.1;
      par[2]=2.5;
      par[3]=1000.;
      return par;
    }

    virtual std::vector<double> get_stepSize(double increment){
      std::vector<double> stepSize(4);
      std::vector<double> minVal=this->get_minVal();
      std::vector<double> maxVal=this->get_maxVal();
      stepSize[0]=increment*(maxVal[0]-minVal[0]);
      stepSize[1]=increment*(maxVal[1]-minVal[1]);
      stepSize[2]=increment*(maxVal[2]-minVal[2]);
      stepSize[3]=increment*(maxVal[3]-minVal[3]);
      return stepSize;
    }

    virtual std::vector<double> get_minVal(){
      std::vector<double> minVal(4);
      minVal[0]=log(1.e-20);
      minVal[1]=log(0.1);
      minVal[2]=log(0.1);
      minVal[3]=log(50.);
      return minVal;
    }
  
    virtual std::vector<double> get_maxVal(){
      std::vector<double> maxVal(4);
      maxVal[0]=log(1.e2);
      maxVal[1]=log(5.);
      maxVal[2]=log(5.);
      maxVal[3]=log(3.e5);
      return maxVal;
    }
  };

  //---------------------------------------------------------------------------
  //  Exponential cutoff spectral model
  //
  //  dN/dE = N_0*exp(-E/E_c)*(E/E_0)^(-gamma)
  //
  //  N_0   = prefactor        [ ph/cm^2/s/MeV ]
  //  gamma = photon index
  //  E_c   = cutoff energy    [ MeV ]
  //  E_0   = energy scale     [ default = 100 MeV ]
  //---------------------------------------------------------------------------

  class ExpCutoff: public SpectralModel{
    double m_prefactor,m_index,m_cutoff;
    double m_prefactor_error,m_index_error,m_cutoff_error;

    class functor{
    private:
      double m_prefactor,m_index,m_cutoff;
    public:
      functor(double prefactor,double index,double cutoff):
	m_prefactor(prefactor),
	m_index(index),
	m_cutoff(cutoff){
      }
      double operator()(double E) const {
	double scale=100.;
	double result=m_prefactor*exp(-1.*E/m_cutoff)*pow(E/scale,-1.*m_index);
	return result;
      }
    };
    
  public:
    ExpCutoff(std::vector<double> params):
      SpectralModel("EXP_CUTOFF"),
      m_prefactor(params[0]),
      m_index(params[1]),
      m_cutoff(params[2]){
    }
    
    ExpCutoff():
      SpectralModel("EXP_CUTOFF"){
      this->set_params(this->get_default_par());
    }

    void set_params(const std::vector<double> params){
      m_prefactor=params[0];
      m_index=params[1];
      m_cutoff=params[2];
    }
    
    std::vector<double> get_params(){
      std::vector<double> params(3);
      params[0]=m_prefactor;
      params[1]=m_index;
      params[2]=m_cutoff;
      return params;
    }
    
    void set_param_errors(const std::vector<double> param_errors){
      m_prefactor_error=param_errors[0];
      m_index_error=param_errors[1];
      m_cutoff_error=param_errors[2];
    }

    std::vector<double> get_param_errors(){
      std::vector<double> param_errors(3);
      param_errors[0]=m_prefactor_error;
      param_errors[1]=m_index_error;
      param_errors[2]=m_cutoff_error;
      return param_errors;
    }

    double get_flux(double E_min,double E_max){
      double error=1.e-5, flux=0;
      int ierr = -1;
      functor arg(m_prefactor,m_index,m_cutoff);    
      flux=st_facilities::GaussianQuadrature::dgaus8<const ExpCutoff::functor>(arg,E_min,E_max,error,ierr);
      if(ierr!=1) throw std::runtime_error("WARNING. Error in dgaus8 integration.");
      return flux;
    }

    virtual double get_dNdE(double E){
      double scale=100.;
      return m_prefactor*exp(-1.*E/m_cutoff)*pow(E/scale,-1*m_index);
    }

    int get_npar(){
      return 3;
    }

    virtual ~ExpCutoff(){};

    virtual void show(){
      std::cout << std::setprecision(4)
		<< std::scientific 
	//<< "N_0 = " << m_prefactor
		<< "flux_100 = " << this->get_flux(100.,3.e5)
		<< std::fixed
		<< "   gamma = " << m_index
	        << "   E_c = " << m_cutoff;
    }

    virtual void show_legend(){
      std::cout << std::endl
		<< "  Exponential cutoff spectral model" << std::endl
		<< std::endl
		<< "  dN/dE = N_0*exp(-E/E_c)*(E/E_0)^(-gamma)" << std::endl
		<< std::endl
		<< "  N_0   = prefactor        [ ph/cm^2/s/MeV ]" << std::endl
		<< "  gamma = photon index" << std::endl
		<< "  E_c   = cutoff energy    [ MeV ]" << std::endl
		<< "  E_0   = energy scale     [ default = 100 MeV ]" << std::endl
		<< std::endl;
    }

    virtual std::vector<std::string> get_parName(){
      std::vector<std::string> parName(3);
      parName[0]="N_0";
      parName[1]="gamma";
      parName[2]="E_c";
      return parName;
    }
    
    virtual std::vector<double> get_default_par(){
      std::vector<double> par(3);
      par[0]=1.e-9;
      par[1]=2.1;
      par[2]=1000.;
      return par;
    }
    
    virtual std::vector<double> get_stepSize(double increment){
      std::vector<double> stepSize(3);
      std::vector<double> minVal=this->get_minVal();
      std::vector<double> maxVal=this->get_maxVal();
      stepSize[0]=increment*(maxVal[0]-minVal[0]);
      stepSize[1]=increment*(maxVal[1]-minVal[1]);
      stepSize[2]=increment*(maxVal[2]-minVal[2]);
      return stepSize;
    }

    virtual std::vector<double> get_minVal(){
      std::vector<double> minVal(3);
      minVal[0]=log(1.e-20);
      minVal[1]=log(1.);
      minVal[2]=log(200.);
      return minVal;
    }
    
    virtual std::vector<double> get_maxVal(){
      std::vector<double> maxVal(3);
      maxVal[0]=log(1.2);
      maxVal[1]=log(5.);
      maxVal[2]=log(5.12e5);
      return maxVal;
    }
  };
  
  class SpectralModelCollection{
  private:
    std::vector<SpectralModel*> m_spectral_pointer_vector;
  public:
    SpectralModelCollection(){
      m_spectral_pointer_vector=std::vector<SpectralModel*>();
    }
    
    void push_back(SpectralModel& spectral_model){
      m_spectral_pointer_vector.push_back(&spectral_model);
    }

    SpectralModel* get(int i){
      return m_spectral_pointer_vector[i];
    }
    
    int size(){
      return m_spectral_pointer_vector.size();
    }
  };

  //---------------------------------------------------------------------------
  //  Probability distribution function of the factors depending on exposure
  //---------------------------------------------------------------------------
  
  class PDFExposure{
  private:
    int m_n;
    double m_mu,m_sigma,m_A;
  public:
    PDFExposure(int n, double mu, double sigma, double A)
      :m_n(n),m_mu(mu),m_sigma(sigma),m_A(A){
    }
    
    void set(int n, double mu, double sigma, double A){
      m_n=n;
      m_mu=mu;
      m_sigma=sigma;
      m_A=A;
    }
    
    double get_log_pdf_exp(){
      double branch=m_mu-(m_A*m_sigma*m_sigma);
      
      double l_pdf=0;
      if(branch<=0){
	double l_num_c=-(m_n*log(2.)/2.)-(m_mu*m_mu/(2.*m_sigma*m_sigma))+(m_n*log(m_A*m_sigma));
	
	double l_den_c=(1/2.)*log(M_PI)+log(1+erf(m_mu/(sqrt(2.)*m_sigma)));
	
	double a=(1+m_n)/2.;
	double b=1/2.;
	double z=branch*branch/(2.*m_sigma*m_sigma);
	
	HypergeometricFunctionU hU(a,b,1e-15,l_num_c);
	
	l_pdf=log(hU(z))-l_den_c;
      }else{
	double l_num_1=(m_n*log(2.)/2.)-(m_mu*m_mu/(2.*m_sigma*m_sigma))+(m_n*log(m_A*m_sigma))+GammaFunction::lngamma((1+m_n)/2.);
	
	double l_num_2=((m_n+1)*log(2.)/2.)-(m_mu*m_mu/(2.*m_sigma*m_sigma))+(m_n*log(m_A*m_sigma))+log(branch)-log(m_sigma)+GammaFunction::lngamma(1+(m_n/2.));
	double l_den_c=(1/2.)*log(M_PI)+log(1+erf(m_mu/(sqrt(2.)*m_sigma)))+GammaFunction::lngamma(1+m_n);
	
	double a_1=(1+m_n)/2.;
	double b_1=1/2.;
	double z_1=branch*branch/(2.*m_sigma*m_sigma);
	
	double a_2=(2.+m_n)/2.;
	double b_2=3./2.;
	double z_2=branch*branch/(2.*m_sigma*m_sigma);
	
	HypergeometricFunction1F1 h1F1_1(a_1,b_1,1e-15,l_num_1-l_den_c);
	HypergeometricFunction1F1 h1F1_2(a_2,b_2,1e-15,l_num_2-l_den_c);
	
	l_pdf=log(h1F1_1(z_1)+h1F1_2(z_2));
      };
      return l_pdf;
    }

    double get_log_simple_pdf(){
      double n_model=m_A*m_mu;
      double l_pdf=-n_model+m_n*log(n_model)-GammaFunction::lngamma(1+m_n);
      return l_pdf;
    }

    double get_simple_pdf(){
      return exp(get_log_simple_pdf());
    }
  
    double get_pdf_exp(){
      return exp(get_log_pdf_exp());
    }
  };

  //---------------------------------------------------------------------------
  //  Total probability distribution function
  //---------------------------------------------------------------------------

  class PDF{
  private:
    
    int m_bin;

    int m_n;
    int m_useExposurePDF;
    double m_E_min,m_E_max,m_mu_exp,m_sigma_exp,m_mu_alpha,m_sigma_alpha;
    
    SpectralModel* m_model_pointer;

    SourceLikelihood* m_source_pointer;

    //-------------------------------------------------------------------------
    //  Helper function to calculate the argument going to numeric integrator
    //-------------------------------------------------------------------------

    class functor{
    private:
      int m_bin;

      int m_n;
      double m_E_min,m_E_max,m_mu_exp,m_sigma_exp,m_mu_alpha,m_sigma_alpha;
      double m_flux;
      
      SpectralModel* m_model_pointer;

      SourceLikelihood* m_source_pointer;

      int m_useExposurePDF;

      int m_normalize;

    public:
      functor(SourceLikelihood* SourcePointer,SpectralModel* ModelPointer,int bin,int useExposurePDF=0){
	m_source_pointer=SourcePointer;
	m_model_pointer=ModelPointer;

	m_bin=bin;

	m_n=m_source_pointer->at(m_bin)->photons();

	m_E_min=m_source_pointer->at(m_bin)->band().emin();
	m_E_max=m_source_pointer->at(m_bin)->band().emax();

	m_mu_alpha=m_source_pointer->at(m_bin)->alpha();
	m_sigma_alpha=m_source_pointer->at(m_bin)->sigma_alpha();

	// Model dependent exposure calculation
	m_mu_exp=m_model_pointer->get_model_exposure(m_source_pointer,m_bin);
	double exposure_error_frac=m_model_pointer->get_exposure_uncertainty(m_E_min,m_E_max);
	m_sigma_exp=exposure_error_frac*m_mu_exp;

	m_flux=m_model_pointer->get_flux(m_E_min,m_E_max);

	m_useExposurePDF=useExposurePDF;
      }

      void normalize(int normalize){ m_normalize=normalize; }

      double operator()(double alpha) const {
	double loglike_best = m_source_pointer->at(m_bin)->operator()(m_mu_alpha);
	double loglike      = m_source_pointer->at(m_bin)->operator()(alpha);

	// Normalization condition
	if(m_normalize) 
	  return exp(loglike_best-loglike);

	double A=m_flux/alpha;

	PDFExposure arg_exp(m_n,m_mu_exp,m_sigma_exp,A);

	// Case 1: No photon counts in bin
	if(m_n<1.){
	  if(m_useExposurePDF)
	    return arg_exp.get_pdf_exp();
	  else
	    return arg_exp.get_simple_pdf();
	}

	// Case 2: Photon counts in bin
	if(m_useExposurePDF)
	  return exp(loglike_best-loglike)*arg_exp.get_pdf_exp();
	else
	  return exp(loglike_best-loglike)*arg_exp.get_simple_pdf();
      }
    };

  public:
    PDF(SourceLikelihood* SourcePointer,SpectralModel* ModelPointer,int useExposurePDF=0):
      m_source_pointer(SourcePointer),
      m_model_pointer(ModelPointer),
      m_useExposurePDF(useExposurePDF){
    }

    void set_bin(int bin){
      m_bin=bin;

      m_n=m_source_pointer->at(m_bin)->photons();

      m_E_min=m_source_pointer->at(m_bin)->band().emin();
      m_E_max=m_source_pointer->at(m_bin)->band().emax();

      m_mu_alpha=m_source_pointer->at(m_bin)->alpha();
      m_sigma_alpha=m_source_pointer->at(m_bin)->sigma_alpha();

      // Model dependent exposure calculation
      m_mu_exp=m_model_pointer->get_model_exposure(m_source_pointer,m_bin);
      double exposure_error_frac=m_model_pointer->get_exposure_uncertainty(m_E_min,m_E_max);
      m_sigma_exp=exposure_error_frac*m_mu_exp;
    }

    double get_likelihood(){
      double error=1.e-8;
      int ierr = -1;
      double norm= 1.;
      double result=0;

      functor arg(m_source_pointer,m_model_pointer,m_bin,m_useExposurePDF);

      // If there are photon counts convolve internal PDF with gaussian distribution of the range alpha values 
      int exit_loop=0;
      while(exit_loop!=1){
	try{

	  // Integration limits
	  double alpha_min,alpha_max;
	  if(m_n>0){
	    alpha_min=m_mu_alpha-5.*m_sigma_alpha;
	    if(alpha_min < 1.e-5) alpha_min=1.e-5;
	    alpha_max=m_mu_alpha+5.*m_sigma_alpha;
	    if(alpha_max > 1.) alpha_max=1.;
	  }
	  else{
	    return arg(1.);
	    //alpha_min=1.e-5;
	    //alpha_max=1.;
	  }

	  // Compute the normalization factor
	  arg.normalize(1);
	  norm=1./st_facilities::GaussianQuadrature::dgaus8<const PDF::functor>(arg,alpha_min,alpha_max,error,ierr);

	  // Compute PDF
	  arg.normalize(0);
	  result=norm*st_facilities::GaussianQuadrature::dgaus8<const PDF::functor>(arg,alpha_min,alpha_max,error,ierr);
      
	  if(ierr!=1) throw std::runtime_error("WARNING. Error in dgaus8 integration.");
	}
	catch(std::runtime_error message){
	  error=error*10.; // Decrease the numeric accuracy and try again
	  if(error>1.e-4) throw std::runtime_error("WARNING. Upper bound on numeric integration accuracy exceeded.");
	  continue;
	}
	exit_loop=1;
      }

      return result;
    }

    double get_loglikelihood(){
      double likelihood=this->get_likelihood();
      if(likelihood < 1.e-280) 
	return -650.; // Limits of numerical accuracy
      else
	return log(likelihood);
    }

  };

  //--------------------------------------------------------------------------
  //  Multiwavelength probability distribution function 
  //--------------------------------------------------------------------------

  class PDFMultiwavelength{
  private:
    double m_E;
    double m_dNdE;
    double m_dNdE_err_lo;
    double m_dNdE_err_hi;

    SpectralModel* m_model_pointer;

    double m_model_dNdE;

  public:
    PDFMultiwavelength(SpectralModel* ModelPointer):
      m_model_pointer(ModelPointer){
    }

    void set_bin(double E,double dNdE,double dNdE_err_lo,double dNdE_err_hi){
      m_E           = E;
      m_dNdE        = dNdE;
      m_dNdE_err_lo = dNdE_err_lo;
      m_dNdE_err_hi = dNdE_err_hi;
      m_model_dNdE  = m_model_pointer->get_dNdE(m_E);
    }

    double get_likelihood(){
      double err_ratio=fabs(m_model_dNdE-m_dNdE)/m_dNdE_err_hi;
      if(err_ratio<5.)
	return 1.-erf(err_ratio/sqrt(2.));
      else
	return 5.733e-7;
    }

    double get_loglikelihood(){
      return log(this->get_likelihood());
    }

  };

  //--------------------------------------------------------------------------
  //  Feldman-Cousins object 
  //--------------------------------------------------------------------------

  class FeldmanCousins{
  private:
    double m_n_total;
    double m_n_signal,m_n_background;
    double m_lngamma;
    double m_step;
    double m_normalization;

  public:
    FeldmanCousins(double n_total=0.,double n_background=0.):
      m_n_signal(n_total-n_background),
      m_n_background(n_background),
      m_n_total(n_total),
      m_lngamma(GammaFunction::lngamma(m_n_total+1.)),
      m_step(0.01),
      m_normalization(0.){
    }

    void set_counts(double n_total,double n_background){
      m_n_signal=n_total-n_background;
      m_n_background=n_background;
      m_n_total=n_total;
      m_lngamma=GammaFunction::lngamma(m_n_total+1.);
    }

    double get_upper_limit(double confidence_limit){		         
      // Correct numeric calculation 
      // Feldman Cousins PDF becomes un-normalized at low counts
      this->set_normalization();
      
      double upper_limit=0;
      double sum_pdf=0.;

      while(sum_pdf<confidence_limit){
	sum_pdf=sum_pdf+m_step*this->pdf(upper_limit)*m_normalization;
	upper_limit+=m_step;
	if(upper_limit>10.*(m_n_total+1)) return -1.;
      }
      return upper_limit;
    }

    double pdf(double mu){
      double log_pdf=m_n_total*log(mu+m_n_background)-(mu+m_n_background)-m_lngamma;
      return exp(log_pdf);
    }

    void set_normalization(){
      double pdf;
      double sum_pdf=0.;
      double i=0.;
      double max_pdf=0.;

      double max_iterations;
      if(m_n_signal/m_n_total<0.1 && m_n_signal>0.)
	max_iterations=10.*m_n_signal;
      else if(m_n_signal<0.)
	max_iterations=10*m_n_background;
      else
	max_iterations=2*m_n_total;

      do{
	pdf=this->pdf(i);
	sum_pdf=sum_pdf+m_step*this->pdf(i);
	i+=m_step;
      }while(i<max_iterations && sum_pdf<1.);

      m_normalization=1/sum_pdf;
    }

    double get_upper_limit(SourceLikelihood* source_pointer,int bin,double confidence_limit){
      double alpha_min=1.e-5;
      double alpha_max=1.;

      double weighted_alpha_sum=0;
      double normalization=0;
      
      double frac,mid_frac,next_frac;
      double alpha,mid_alpha,next_alpha;
      double weight,dalpha;
      double exposure;
      
      double loglike_best = source_pointer->at(bin)->operator()();
      double loglike;    

      double max=100.;

      for(double i=0.;i<max;i+=1.){
	frac=i/max;
	mid_frac=(i+0.5)/max;
	next_frac=(i+1.)/max;

	alpha=exp((1.-frac)*log(alpha_min)+frac*log(alpha_max));
	mid_alpha=exp((1.-mid_frac)*log(alpha_min)+mid_frac*log(alpha_max));
	next_alpha=exp((1.-next_frac)*log(alpha_min)+next_frac*log(alpha_max));

	dalpha=next_alpha-alpha;

	loglike=source_pointer->at(bin)->operator()(mid_alpha);
	weight=exp(loglike_best-loglike);

	weighted_alpha_sum+=weight*mid_alpha*dalpha;
	normalization+=weight*dalpha;
      }      
      double best_alpha=weighted_alpha_sum/normalization;

      double n_total=source_pointer->at(bin)->photons();

      this->set_counts(n_total,(1.-best_alpha)*n_total);
      return this->get_upper_limit(confidence_limit);
    }

  };

  //--------------------------------------------------------------------------
  //  Read a textfile of multiwavelength data
  //--------------------------------------------------------------------------

  class MWData{
  private:
    std::vector<double> m_E;
    std::vector<double> m_dNdE;
    std::vector<double> m_dNdE_err_lo;
    std::vector<double> m_dNdE_err_hi;

  public:
    MWData(char filename[100],double scale=1.){
      this->read(filename,scale);
    }

    void read(char filename[100],double scale=1.){
      
      std::cout << "Opening " << filename << "..." << std::endl;

      char line[100];
      std::ifstream inFile;
      inFile.open(filename);

      if (!inFile) {
        std::cout << "WARNING: Unable to open multiwavelength data file." << std::endl;
        exit(1);
      }

      std::stringstream linestream;
      std::string linestring;

      std::vector<double> numbers;

      while(!inFile.eof()){
	inFile.getline(line,100);
	linestring.assign(line,strlen(line));
	std::cout << linestring << std::endl;
	if(linestring.find("#")==std::string::npos && !linestring.empty()){
	   numbers=numberize(line);
	   m_E.push_back(numbers[0]*scale);
	   m_dNdE.push_back(numbers[1]/scale);
	   m_dNdE_err_lo.push_back(numbers[2]/scale);
	   if(numbers.size()<4)
	     m_dNdE_err_hi.push_back(numbers[2]/scale);
	   else
	     m_dNdE_err_hi.push_back(numbers[3]/scale);
	}
      }

      inFile.close();
    }

    std::vector<double> numberize(std::string line){
      
      std::string whitespaces(" \t\f\v\n\r");
      size_t found;
      found=line.find_last_not_of(whitespaces);
      if (found!=std::string::npos)
	line.erase(found+1);
      else
	line.clear(); 

      std::vector<double> numbers;

      size_t begin=0;
      size_t end=0;
      std::string number_str;

      while(end<line.length() && begin<line.length()){
	begin=line.find_first_not_of(" ",end);
	end=line.find_first_of(" ",begin+1);
	number_str=line.substr(begin,end-begin);
	numbers.push_back(atof(number_str.c_str()));
      }

      return numbers;
    }

    std::vector<double> get_E(){ return m_E; }
    std::vector<double> get_dNdE(){ return m_dNdE; }
    std::vector<double> get_dNdE_err_lo(){ return m_dNdE_err_lo; }
    std::vector<double> get_dNdE_err_hi(){ return m_dNdE_err_hi; }

  };

  //--------------------------------------------------------------------------
  //  Class definition
  //--------------------------------------------------------------------------

  class SpectralFitter{
  private:

    static double s_lower_bound,s_upper_bound;

    static int s_useDefaultParams;
    static int s_useSimplex;
    static int s_useGradient;
    static int s_useMinos;

    static int s_useExposurePDF;

    static int s_useMultiwavelengthData;

    static double s_stepIncrement;

    static double s_accuracy;

    static int s_useUpperLimit;
    static double s_TS_threshold;
    static double s_index;
    static double s_upper_limit_lower_bound;
    static double s_upper_limit_upper_bound;
    static double s_band_confidence_limit;

    int m_combined;

    int m_npar;

    std::vector<double> m_specParams;
    std::vector<double> m_specParamErrors;
    
    std::vector<double> m_covar_entries;

    std::vector<std::string> m_parName;
    std::vector<double> m_par;
    std::vector<double> m_stepSize;
    std::vector<double> m_minVal;
    std::vector<double> m_maxVal;

    TMinuit* m_Minuit;

    std::vector<double> m_model_energies;
    std::vector<double> m_exposures;
    
    FeldmanCousins* m_FeldmanCousins;

    SpectralModel* m_model_pointer;

    MWData* m_MWData;

    std::vector<double> m_confidence_limits;
    std::vector<double> m_flux_upper_limits;

    std::vector<double> m_band_upper_limits;
    std::vector<double> m_energy_upper_limits;
    std::vector<double> m_exposure_upper_limits;

  public:
    SpectralFitter(SourceLikelihood& source, SpectralModel& model);

    ~SpectralFitter();

    void initialize();
    void specfitMinuit(double scale=-1.);
    void setFluxUpperLimits(std::vector<double> confidence_limits=NULL);

    // Functions for setting and checking spectral fitting parameters
    
    void setFitRange(double lower_bound,double upper_bound);
    void setFluxUpperLimitRange(double lower_bound,double upper_bound);
    void setConfidenceLimits(std::vector<double> confidence_limits=NULL);
    void useExposurePDF(int useExposurePDF=1) { s_useExposurePDF=useExposurePDF; }; 

    // Function for using combined front and back energy binning

    void setCombined();

    // Function deciding whether to perform the upper limit calculation

    void useUpperLimit(int useUpperLimit=1) { s_useUpperLimit=useUpperLimit; };

    // Function to get covariance matrix entries
    
    std::vector<double> get_covar_entries() { return m_covar_entries; };

    // Functions for exporting upper limit results

    std::vector<double> getUpperLimitRange(){
      std::vector<double> range(2);
      range[0]=s_upper_limit_lower_bound;
      range[1]=s_upper_limit_upper_bound;
      return range;
    }

    std::vector<double> getConfidenceLimits() { return m_confidence_limits; };
    std::vector<double> getFluxUpperLimits() { return m_flux_upper_limits; };

    std::vector<double> getBandUpperLimits() { return m_band_upper_limits; };
    std::vector<double> getEnergyUpperLimits() { return m_energy_upper_limits; };
    std::vector<double> getExposureUpperLimits() { return m_exposure_upper_limits; };

    // Get the model pointer from the spectral fitter

    SpectralModel* getModel() { return m_model_pointer; };

    // Functions to add multiwavelength data to the spectral fitting

    void setMWData(char filename[100],double scale=1.) { 
      m_MWData = new MWData(filename,scale); 
      s_useMultiwavelengthData = 1;
    };
    MWData* getMWData() { return m_MWData; };

  };

}

#endif
