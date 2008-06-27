#ifndef ExtendedSourcePDF_h
#define ExtendedSourcePDF_h


#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

#include "pointlike/HypergeometricFunction.h"

#include "skymaps/PsfFunction.h"

#include "astro/SkyDir.h"

namespace pointlike{
    namespace extendedSources { 

    //--------------------------------------------------------------------------------------------------
    //  iSource
    //   
    //  Abstract base class for a source shape to be used in extendedSourcePDF
    //  operator() returns a radial profile in scaled space v=0.5 r**2/s**2
    //  min() returns a lower (integration) limit for the source radius
    //  max() returns an upper (integration) limit for the source radius
    //  set() is used to set the parameters for the radial profile of the source
    //--------------------------------------------------------------------------------------------------

	class iSource {
	   private:
	     double m_sigma;
             double m_vmax,m_vmin;

	   public:
             iSource(double sgm, double vmax=0.,double vmin=0.): m_sigma(sgm),m_vmax(vmax) , m_vmin(vmin) {};

             virtual void set(const std::vector<double>& param)=0;
             virtual double operator() (double v) const  =0;
             virtual double grad(double v,int icomp) const =0;
	     virtual double get(int iparam) const {  return NAN; iparam=0;};
	     virtual bool parameterShiftsMax(int iparam) const { return false; iparam=0; };
	     virtual bool parameterShiftsMin(int iparam) const { return false; iparam=0; };
             double min() const {return m_vmin;};
             double max() const {return m_vmax;};
             double sigma() const {return m_sigma;};
             double sigma2() const {return m_sigma*m_sigma;};
             void   min(double vmin){m_vmin=vmin;};
             void   max(double vmax){m_vmax=vmax;};
             void   sigma(double sgm){m_sigma=sgm;};
	};  

    //--------------------------------------------------------------------------------------------------
    //  diskSource
    //   
    //  a disk-like source, i.e. a constant-brightness disk/ring between rmin and rmax
    //--------------------------------------------------------------------------------------------------

	class diskSource: public iSource {
	   double m_rmax;
	   double m_rmin;

	   public:  
	     diskSource(double sigma, double rmax, double rmin=0.):iSource(sigma),m_rmax(rmax),m_rmin(rmin){
	        if(rmax>=0 && rmin>=0){
	  	   max(0.5*rmax*rmax/sigma2());
	           min(0.5*rmin*rmin/sigma2());
	          std::cout<<sigma<<" "<<min()<<" "<<max()<<std::endl;
	        } else max(rmax);
		   
	     };
	     void set(const std::vector<double>& param){ 
		if (param.size()>0) { max(0.5*param[0]*param[0]/sigma2()); m_rmax=param[0];};
		if (param.size()>1) { min(0.5*param[1]*param[1]/sigma2()); m_rmin=param[1]; };
//                std::cout<<sigma()<<" "<<min()<<" "<<max()<<" "<<param[0]<<" "<<param[1]<<std::endl;
	     };

	     bool parameterShiftsMax(int iparam) const {if(iparam==0) return true; return false; };
	     bool parameterShiftsMin(int iparam) const {if(iparam==1) return true; return false; };

	     double get(int ipar) const { 
		if (ipar==0) return m_rmax;
		if (ipar==1) return m_rmin;
		return NAN;
	     };

	    double operator()(double v) const { 
		return 1./(2*M_PI*(max()-min()));
		v=0;
	    };

	    double grad(double v,int icomp) const { 
	        double dm=max()-min();
		if (icomp==0)  return  -1./(2*M_PI*dm*dm)*m_rmax/sigma2();
		if (icomp==1) return   1./(2*M_PI*dm*dm)*m_rmin/sigma2();
		return NAN;
		v=0;
	    };

	};     

    //--------------------------------------------------------------------------------------------------
    //  gaussSource
    //   
    //  an extended source with a gaussian radial intensity profile and radius r 
    //--------------------------------------------------------------------------------------------------

	class gaussSource: public iSource {
	   mutable double m_radius;
	   double m_intLimit;
	   double m_intLimitNorm;
           mutable double m_scaledRadius;
	   public:  
	     gaussSource(double sigma, double rad): iSource(sigma),m_radius(rad),m_intLimit(4.6052) {  //int limit 4.6052 equals to 99% boundary..,.
		m_scaledRadius=0.5*m_radius*m_radius/(sigma2());
		max(m_intLimit*m_scaledRadius);
		m_intLimitNorm=1.-exp(-m_intLimit); 
	     };
	     void set(const std::vector<double>& param){ 
		if (param.size()>0) m_radius=param[0]; 
		m_scaledRadius=0.5*m_radius*m_radius/sigma2();
		max(m_intLimit*m_scaledRadius); 
	     };
	     
	     double get(int ipar) const {return m_radius; ipar=0;};
	     
	     double operator() (double v) const { return exp(-v/m_scaledRadius)/(2*M_PI*m_scaledRadius*m_intLimitNorm);};
	     
	     double grad(double v,int icomp) const { 
	         if (icomp==0) return (v/m_scaledRadius - 1) * exp(-v/m_scaledRadius)  /
		                              (2*M_PI*m_scaledRadius*m_scaledRadius*m_intLimitNorm) *m_radius/sigma2();
                 return NAN;
	     };
	};     

    //--------------------------------------------------------------------------------------------------
    //  generator
    //   
    //  a factory method to generate the different source types by name
    //--------------------------------------------------------------------------------------------------

	static iSource& generator(const std::string & name, double sigma){
	     if (name=="point") return *(new diskSource(sigma,-1.));
	     if (name=="pseudopoint") return *(new diskSource(sigma,1.e-4));
	     if (name=="disk") return *(new diskSource(sigma,1.e-2));
	     if (name=="gauss") return *(new gaussSource(sigma,1.e-2));
	     std::cerr<<"WARNING in extendedSources::generator : Source name "<<name<<" invalid. Returning point source."<<std::endl;
             return *(new diskSource(sigma,1.e-5));	 
	 };	 

/*
	static iSource& generator(const std::string & name, double sigma, const std::vector<double> params){
             iSource& src=generator(name,sigma);
	     src.set(params);
	     return src;
	 };	 
*/
    };

    //--------------------------------------------------------------------------------------------------
    //  ExtendedSourcePDF
    //   
    //  Semi-analytic calculation of the PDF for an extended source given the GLAST PSF parametrization 
    //  v corresponds to the radial coordinate of the source extension in scaled space v= 0.5 * r**2 / sigma**2
    //  u is the deviation from the source center u= 0.5 * r**2 / sigma**2
    //  the integral over v is performed numerically with dgaus8
    //  the integral over azimuth angles was done analytically
    //--------------------------------------------------------------------------------------------------

    class extendedSourcePDF {

       extendedSources::iSource& m_source;

       double m_gamma,m_prefactor,m_a,m_b;
       mutable double m_u;
       double m_gamma_pow_gamma;

       HypergeometricFunction2F1 m_hF;
       HypergeometricFunction2F1 m_dHF;

     public:

       extendedSourcePDF& operator=(const extendedSourcePDF& other){
          m_source=other.m_source;
	  m_gamma=other.m_gamma;
          m_prefactor=m_prefactor;
	  m_a=other.m_a;
	  m_b=other.m_b;
          m_u=other.m_u;
          m_hF.set(m_a,m_b,1.);
          m_dHF.set(m_a+1.,m_b+1.,2.);
	  return *this;
       };
       
       extendedSourcePDF(extendedSources::iSource& source,double g=2.0);

       double value(double u) const ;
       std::vector<double> gradient(double u, int npar) const;
       
       double operator() (double v) const ;

       double gpdf(double v) const ;
       double du_gpdf(double v) const ;

     private:
       class gradUFunctor {
             const extendedSources::iSource& m_source;
	     const extendedSourcePDF& m_pdf;
           public:
	     gradUFunctor(const extendedSources::iSource& src, const extendedSourcePDF& pdf):m_source(src),m_pdf(pdf){};  
             double operator()(double v) const { return m_source(v)*m_pdf.du_gpdf(v); };
	};     
       
        class gradSFunctor {
             const extendedSources::iSource& m_source;
	     const extendedSourcePDF& m_pdf;
	     int component;
           public:
	     gradSFunctor(const extendedSources::iSource& src,const extendedSourcePDF& pdf,int i):m_source(src),m_pdf(pdf),component(i){};  
             double operator()(double v) const { return m_source.grad(v,component)*m_pdf.gpdf(v); };
	};     
     };
 


    //--------------------------------------------------------------------------------------------------
    //  ExtendedSourcePseudoPsf
    //   
    //  Sell the PDF for an extended source to pointlike as an adjusted PSF
    //  This allows to keep most of the pointlike structures unchanged. Only minimal adjustments in the
    //  parameter handling are necessary 
    //--------------------------------------------------------------------------------------------------


    class ExtendedSourcePseudoPSF {
	public:
            /// @param gamma the power-law factor
            ExtendedSourcePseudoPSF(extendedSources::iSource& src, double g=2.0)
        	: m_source(src)
		, m_gamma(g)
		, m_pdf(src,g)
		, m_pointSourcePDF(g)
		, m_isPointSource(false)
	    {
		if (m_source.max()<0) m_isPointSource=true;
		
	    };

            extendedSources::iSource& source(){return m_source;};

            extendedSourcePDF& pdf()            {return m_pdf;};
            skymaps::PsfFunction&       pointSourcePdf() {return m_pointSourcePDF;};
            
            double gamma() const {return m_gamma;};
	    extendedSources::iSource& source() const { return m_source;};
	    
            void setGamma(double gamma){ 
	        m_gamma=gamma; 
		m_pdf=extendedSourcePDF(m_source,m_gamma);
		m_pointSourcePDF=skymaps::PsfFunction(m_gamma);
	    };
            
            void setSource(extendedSources::iSource& src){ 
	        m_source=src; 
		m_pdf=extendedSourcePDF(m_source,m_gamma);
	    };
	    
            double operator () (const astro::SkyDir & r, const astro::SkyDir & r_prime, double sigma) ;
            std::vector<double> gradient(double u,int ncomp=0) const ;

            ///@return the value as a function of the scaled deviation squared
            double operator()(double u) const;   

	    double integral(double umax) const { return m_pointSourcePDF.integral(umax);};
	    double integralSquare(double umax) const {return m_pointSourcePDF.integralSquare(umax);};

	private:
	    extendedSources::iSource& m_source;
            double m_gamma;
	    double m_umax;
            double m_norm;
	    extendedSourcePDF m_pdf;
            skymaps::PsfFunction m_pointSourcePDF;
            bool m_isPointSource;

    };

};

#endif

