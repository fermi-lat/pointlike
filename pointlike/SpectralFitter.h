#ifndef tools_SpectralFitter_h
#define tools_SpectralFitter_h

#include "pointlike/HypergeometricFunction.h"
#include "st_facilities/GaussianQuadrature.h"
#include "pointlike/SourceLikelihood.h"

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <cmath>
#include <vector>
#include <string>
#include <cstring>
#include <stdlib.h>

#include <cstdlib>
#include <time.h>

#include <TMinuit.h>
#include <TF1.h>

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
    double m_E_min,m_E_max;

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
      
      double frac,E,dNdE,exposure;
      
      double max=100.;

      for(double i=0.;i<max;i+=1.){
	frac=i/max;
	E=frac*E_max+(1-frac)*E_min;
	dNdE=this->get_dNdE(E);
	exposure=SourcePointer->at(bin)->exposure(E);
	weighted_exposure_sum+=dNdE*exposure;
	normalization+=dNdE;
      }
      
      return weighted_exposure_sum/normalization;
    }

    double get_model_E(double E_min,double E_max){
      if(E_min<this->get_lower_bound() || E_max>this->get_upper_bound()){
	std::cout << "WARNING. Exposure calculation out of range." << std::endl;
	return -1.;
      }

      double weighted_E_sum=0;
      double normalization=0;

      double frac,E,dNdE;

      double max=100.;
      for(double i=0.;i<max;i+=1.){
	frac=i/max;
	E=frac*E_max+(1.-frac)*E_min;
	dNdE=this->get_dNdE(E);
	weighted_E_sum+=dNdE*E;
	normalization+=dNdE;
      }
      
      return weighted_E_sum/normalization;
    }

    double get_exposure_uncertainty(double E_min,double E_max){
      // Until we get a better estimate from the C&A group
      double E_log_center=exp((log(E_max)+log(E_min))/2.);
      if(E_log_center<1.e2) return 0.2;
      if(E_log_center<1.e4) return 0.1;
      else return 0.3;
    }

    virtual void get_E2dNdE(TF1& func,std::vector<double> params)=0;

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

    void set_energy_range(double E_min,double E_max){
      m_E_min=E_min;
      m_E_max=E_max;
    }
    double get_lower_bound(){ return m_E_min; }
    double get_upper_bound(){ return m_E_max; }
  };

  //---------------------------------------------------------------------------
  //  Power law spectral model
  //
  //  dN/dE = N_0*(E/E_0)^(-gamma)
  //
  //  N_0   = prefactor        [ ph/cm^2/s/MeV ]
  //  E_0   = energy scale     [ default = 100 MeV ]
  //  gamma = spectral index
  //---------------------------------------------------------------------------

  class PowerLaw: public SpectralModel{
    double m_prefactor,m_index;
    double m_prefactor_error,m_index_error;

    double function(double E){
      double scale=100.;
      double result=m_prefactor*E*pow(E/scale,-1.*m_index)/(1.-m_index);
      return result; 
    }
  public:
    PowerLaw(std::vector<double> params):
      SpectralModel("POWER_LAW"),
      m_prefactor(params[0]),
      m_index(params[1]){
    }

    PowerLaw():
      SpectralModel("POWER_LAW"){
      this->set_params(this->get_default_par());
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
      double scale=100.;
      return m_prefactor*pow(E/scale,-1*m_index);
    }

    virtual void get_E2dNdE(TF1& func,std::vector<double> params){
      double scale=100.;
      func=TF1("f1","pow(x,2)*[0]*pow(x/[2],-1*[1])",10,5e5);
      func.SetParameters(params[0],params[1],scale);
    }

    int get_npar(){
      return 2;
    }

    virtual ~PowerLaw(){};

    virtual void show(){
      std::cout << std::setprecision(4)
		<< std::scientific 
		<< "N_0 = " << m_prefactor
		<< std::fixed 
		<< "   gamma =  " << m_index;
    }

    virtual void show_legend(){
      std::cout << std::endl
		<< "  Power law spectral model" << std::endl
		<< std::endl
		<< "  dN/dE = N_0*(E/E_0)^(-gamma)" << std::endl
		<< std::endl
		<< "  N_0   = prefactor        [ ph/cm^2/s/MeV ]" << std::endl
		<< "  E_0   = energy scale     [ default = 100 MeV ]" << std::endl
		<< "  gamma = spectral index" << std::endl
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
      minVal[1]=log(1.);
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
  //  gamma_1 = spectral index in lower energy regime
  //  gamma_2 = spectral index in higher energy regime
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

    virtual void get_E2dNdE(TF1& func,std::vector<double> params){
      func=TF1("f1","(x<[3])*pow(x,2)*[0]*pow(x/[3],-1*[1])+(x>[3])*pow(x,2)*[0]*pow(x/[3],-1*[2])",10,5e5);
      func.SetParameters(params[0],params[1],params[2],params[3]);
    }

    int get_npar(){
      return 4;
    }

    virtual ~BrokenPowerLaw(){};

    virtual void show(){
      std::cout << std::setprecision(4)
		<< std::scientific 
		<< "N_0 = " << m_prefactor
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
		<< "  gamma_1 = spectral index in lower energy regime" << std::endl
		<< "  gamma_2 = spectral index in higher energy regime" << std::endl
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
      minVal[1]=log(1.);
      minVal[2]=log(1.);
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
  //  gamma = spectral index
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

    virtual void get_E2dNdE(TF1& func,std::vector<double> params){
      double scale=100.;
      func=TF1("f1","pow(x,2)*[0]*exp(-1.*x/[2])*pow(x/[3],-1*[1])",10,5e5);
      func.SetParameters(params[0],params[1],params[2],scale);
    }

    int get_npar(){
      return 3;
    }

    virtual ~ExpCutoff(){};

    virtual void show(){
      std::cout << std::setprecision(4)
		<< std::scientific 
		<< "N_0 = " << m_prefactor
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
		<< "  gamma = spectral index" << std::endl
		<< "  E_c   = cutoff energy    [ MeV ]" << std::endl
		<< "  E_0   = energy scale     [ default = 100 MeV ]" << std::endl
		<< std::endl;
    }

    virtual std::vector<std::string> get_parName(){
      std::vector<std::string> parName(3);
      parName[0]="N";
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
  
    double get_pdf_exp(){
      return exp(get_log_pdf_exp());
    }
  };

  //---------------------------------------------------------------------------
  //  Total probability distribution function
  //---------------------------------------------------------------------------

  class PDF{
  private:
    int m_n;
    double m_E_min,m_E_max,m_mu_exp,m_sigma_exp,m_mu_alpha,m_sigma_alpha;
    
    SpectralModel* m_model_pointer;

    double m_flux;

    //-------------------------------------------------------------------------
    //  Helper function to calculate the argument going to numeric integrator
    //-------------------------------------------------------------------------

    class functor{
    private:
      int m_n;
      double m_E_min,m_E_max,m_mu_exp,m_sigma_exp,m_mu_alpha,m_sigma_alpha;
      double m_flux;
      
    public:
      functor(int n,double E_min,double E_max,double mu_exp,double sigma_exp,double mu_alpha,double sigma_alpha,double flux=0):
	m_n(n),
	m_E_min(E_min),
	m_E_max(E_max),
	m_mu_exp(mu_exp),
	m_sigma_exp(sigma_exp),
	m_mu_alpha(mu_alpha),
	m_sigma_alpha(sigma_alpha),
	m_flux(flux){
      }

      double operator()(double alpha) const {
	
	double result=0;

	double A=m_flux/alpha;

	PDFExposure arg_exp(m_n,m_mu_exp,m_sigma_exp,A);

	double g_alpha_num=sqrt(2./M_PI)*exp(-1.*(alpha-m_mu_alpha)*(alpha-m_mu_alpha)/(2.*m_sigma_alpha*m_sigma_alpha));
	
	double g_alpha_den=m_sigma_alpha*(-erf((-1.+m_mu_alpha)/(sqrt(2.)*m_sigma_alpha))+erf(m_mu_alpha/(sqrt(2.)*m_sigma_alpha)));

	result=g_alpha_num*arg_exp.get_pdf_exp()/g_alpha_den;

	return result;
      }
    };

  public:
    PDF(SpectralModel* ModelPointer):
      m_model_pointer(ModelPointer){
    }

    void set_bin(int n,double E_min,double E_max,double mu_exp,double sigma_exp,double mu_alpha,double sigma_alpha){
      m_n=n;
      m_E_min=E_min;
      m_E_max=E_max;
      m_mu_exp=mu_exp;
      m_sigma_exp=sigma_exp;
      m_mu_alpha=mu_alpha;
      m_sigma_alpha=sigma_alpha;
    }

    double get_likelihood(){
      double error=1.e-8, result=0;
      int ierr = -1;
      
      m_flux=m_model_pointer->get_flux(m_E_min,m_E_max);
      
      functor arg(m_n,m_E_min,m_E_max,m_mu_exp,m_sigma_exp,m_mu_alpha,m_sigma_alpha,m_flux);
    
      int exit_loop=0;
      while(exit_loop!=1){
	try{
	  result=st_facilities::GaussianQuadrature::dgaus8<const PDF::functor>(arg,1e-15,1.,error,ierr);
      
	  if(ierr!=1) throw std::runtime_error("WARNING. Error in dgaus8 integration.");
	}
	catch(std::runtime_error message){
	  error=error*10.; // Decrease the numeric accuracy and try again
	  if(error>100.) throw std::runtime_error("WARNING. Lower bound on numeric integration accuracy exceeded.");
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
  //  Class definition
  //--------------------------------------------------------------------------

  class SpectralFitter{
  private:

    static double s_emin,s_emax;

    static int s_useDefaultParams;
    static int s_useSimplex;
    static int s_useGradient;
    static int s_useMinos;

    static double s_stepIncrement;

    static double s_accuracy;

    int m_npar;

    std::vector<double> m_specParams;
    std::vector<double> m_specParamErrors;

    std::vector<std::string> m_parName;
    std::vector<double> m_par;
    std::vector<double> m_stepSize;
    std::vector<double> m_minVal;
    std::vector<double> m_maxVal;

    TMinuit* m_Minuit;

    std::vector<double> m_model_energies;
    std::vector<double> m_exposures;
    
  public:
    SpectralFitter(SourceLikelihood& source, SpectralModel& model);

    ~SpectralFitter();

    void specfitMinuit();

  };

}

#endif
