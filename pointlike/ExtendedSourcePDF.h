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
	     virtual double get(int) const {  return NAN;};
	     virtual std::vector<double> get() const {  return std::vector<double>(0);};
	     virtual double split(int) const {  return -1;};
	     virtual bool parameterShiftsMax(int) const { return false;  };
	     virtual bool parameterShiftsMin(int) const { return false;  };
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
	     diskSource():iSource(-1.),m_rmax(-1.),m_rmin(-1.){max(-1);min(-1);};

	     diskSource(double sigma, double rmax, double rmin=0.):iSource(sigma),m_rmax(rmax),m_rmin(rmin){
	        if(rmax>=0 && rmin>=0){
	  	   max(0.5*rmax*rmax/sigma2());
	           min(0.5*rmin*rmin/sigma2());
//	          std::cout<<sigma<<" "<<min()<<" "<<max()<<std::endl;
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

	    std::vector<double> get() const { 
		std::vector <double> result(2);
		result[0]=m_rmax;
		result[1]=m_rmin;
		return result;
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
	     gaussSource(double sigma, double rad):iSource(sigma),m_radius(rad),m_intLimit(30.) {  //int limit 4.6052 equals to 99% boundary..,.
		m_scaledRadius=m_radius*m_radius/(sigma2());
		max(m_intLimit*m_scaledRadius);
		m_intLimitNorm=1.-exp(-m_intLimit); 
	     };
	     void set(const std::vector<double>& param){ 
		if (param.size()>0) m_radius=param[0]; 
		m_scaledRadius=m_radius*m_radius/sigma2();
		max(m_intLimit*m_scaledRadius); 
	     };
	     
	     double get(int) const {return m_radius;};
	     
	     std::vector<double> get() const { 
		std::vector <double> result(1);
		result[0]=m_radius;
		return result;
	     };

	     double operator() (double v) const { return exp(-v/m_scaledRadius)/(2*M_PI*m_scaledRadius*m_intLimitNorm);};
	     
	     double grad(double v,int icomp) const { 
	         if (icomp==0) return (v/m_scaledRadius - 1) * exp(-v/m_scaledRadius)  /
		                              (2*M_PI*m_scaledRadius*m_scaledRadius*m_intLimitNorm) *2.*m_radius/sigma2();
                 return NAN;

	     };
      	     // double split(int) const {return m_scaledRadius; };
	};     

    //--------------------------------------------------------------------------------------------------
    //  radiusGaussSource
    //   
    //  an extended source with a radial profile that has 68% of events within 1 sigma 
    //  introduced to fit the GaussianSource produced by gtobssim, which has no true gaussian shape
    //--------------------------------------------------------------------------------------------------

	class radiusGaussSource: public iSource {
	   mutable double m_radius;
	   double m_intLimit;
	   double m_intLimitNorm;
           mutable double m_scaledRadius;
	   public:  
	     radiusGaussSource(double sigma, double rad):iSource(sigma),m_radius(rad),m_intLimit(30.) {  //int limit 4.6052 equals to 99% boundary..,.
		m_scaledRadius=m_radius*m_radius/(sigma2());
		max(m_intLimit*m_scaledRadius);
		m_intLimitNorm=1.-exp(-m_intLimit); 
	     };
	     void set(const std::vector<double>& param){ 
		if (param.size()>0) m_radius=param[0]; 
		m_scaledRadius=m_radius*m_radius/sigma2();
		max(m_intLimit*m_scaledRadius); 
	     };
	     
	     double get(int) const {return m_radius;};
	     
	     std::vector<double> get() const { 
		std::vector <double> result(1);
		result[0]=m_radius;
		return result;
	     };

	     double operator() (double v) const 
	        { return exp(-v/m_scaledRadius)/(2*M_PI*sqrt(M_PI*v*m_scaledRadius)*m_intLimitNorm);};
	     
	     double grad(double v,int icomp) const { 
	         if (icomp==0) return (2. * v/ m_scaledRadius -1. ) * exp(-v/m_scaledRadius)  /
		                      (4.*m_intLimitNorm*M_PI*m_scaledRadius*sqrt(M_PI*m_scaledRadius*v)) *
			              2.*m_radius/sigma2();
                 return NAN;
	     };

     	     double split(int) const {return 0.5*m_scaledRadius; };

	};     


    //--------------------------------------------------------------------------------------------------
    //  generator
    //   
    //  a factory method to generate the different source types by name
    //--------------------------------------------------------------------------------------------------

	static iSource& generator(const std::string & name, double sigma){
	     if (name=="point") return *(new diskSource(sigma,-1.));
	     if (name=="calib") return *(new diskSource(sigma,-1.));
	     if (name=="pseudopoint") return *(new diskSource(sigma,1.e-4));
	     if (name=="disk") return *(new diskSource(sigma,1.e-2));
	     if (name=="gauss") return *(new gaussSource(sigma,1.e-2));
	     if (name=="radius_gauss") return *(new radiusGaussSource(sigma,1.e-2));
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

            void operator = (const ExtendedSourcePseudoPSF& s){
	       m_source=s.m_source;
               m_gamma=s.m_gamma;
	       m_umax=s.m_umax;
               m_norm=s.m_norm;
	       m_pdf=s.m_pdf;
               m_pointSourcePDF=skymaps::PsfFunction(m_gamma);
               m_isPointSource=s.m_isPointSource;
	    };

            extendedSources::iSource& source(){return m_source;};

            const extendedSourcePDF& pdf() const {return m_pdf;};
            const skymaps::PsfFunction& pointSourcePdf() const {return m_pointSourcePDF;};
            
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


    class ExtendedSourcePseudoPSF2 {
        static extendedSources::diskSource defaultSource;
    
	public:
            /// @param gamma the power-law factor
            ExtendedSourcePseudoPSF2(extendedSources::iSource& src,  double g1=2.0, double s1=-1,
	                             extendedSources::iSource& src2=defaultSource, 
				     double g2=-1., double s2=-1.,double f2=0.)
        	: m_source(src)
		, m_source2(src2)
		, m_gamma(g1)
		, m_gamma2(g2)
		, m_sigma(s1)
		, m_sigma2(s2)
		, m_frac2(f2)
		, m_psf(src,g1)
		, m_psf2(src2,g2)
		, m_pointSourcePDF(g1)
		, m_pointSourcePDF2(g2)
		, m_isPointSource(false)
	    {
		if (m_source.max()<0) m_isPointSource=true;
		if(m_sigma<0) m_sigma=m_source.sigma();
		m_s2ratio=(s1*s1/(s2*s2));
//		std::cout<<"PSF s2ratio="<<m_s2ratio<<" gamma="<<m_gamma<<" sigma="<<m_sigma<<" gamma2="<<m_gamma2<<" sigma2="<<m_sigma2<<std::endl;

	    };

            extendedSources::iSource& source(){return m_source;};

            const extendedSourcePDF& pdf() const {return m_psf.pdf();};
            const extendedSourcePDF& pdf2() const {return m_psf2.pdf();};
	    
            skymaps::PsfFunction&       pointSourcePdf() {return m_pointSourcePDF;};
            skymaps::PsfFunction&       pointSourcePdf2() {return m_pointSourcePDF2;};
            
            double gamma() const {return m_gamma;};
            double gamma2() const {return m_gamma2;};
            double sigma() const {return m_sigma;};
            double sigma2() const {return m_sigma2;};
	    
	    extendedSources::iSource& source() const { return m_source;};
	    extendedSources::iSource& source2() const { return m_source2;};
	    
            void setGamma(double gamma,double gamma2=-1){ 
	        m_gamma=gamma; 
	        m_gamma2=gamma2; 
		m_psf = ExtendedSourcePseudoPSF(m_source,m_gamma);
		m_psf2= ExtendedSourcePseudoPSF(m_source2,m_gamma2);
		m_pointSourcePDF=skymaps::PsfFunction(m_gamma);
		m_pointSourcePDF2=skymaps::PsfFunction(m_gamma2);
	    };
            
            void setSource(extendedSources::iSource& src,extendedSources::iSource& src2=defaultSource){ 
	        m_source=src; 
	        m_source2=src2;
		m_psf=ExtendedSourcePseudoPSF(m_source,m_gamma);
		m_psf2=ExtendedSourcePseudoPSF(m_source2,m_gamma2);
	    };
	    
            double operator () (const astro::SkyDir & r, const astro::SkyDir & r_prime, double sigma) ;
            
	    std::vector<double> gradient(double u,int ncomp=0) const ;

            ///@return the value as a function of the scaled deviation squared
            double operator()(double u) const;  

	    double integral(double umax) const { 
	        if(m_gamma2>0) return  (1-m_frac2)*m_pointSourcePDF.integral(umax)
		                      +  m_frac2  *m_pointSourcePDF2.integral(m_s2ratio*umax);
		else return m_pointSourcePDF.integral(umax);
	     };
	    double integralSquare(double umax) const {
	        if(m_gamma2>0) return  (1-m_frac2)*m_pointSourcePDF.integralSquare(umax)
		                      +  m_frac2  *m_pointSourcePDF2.integralSquare(m_s2ratio*umax);
		else return m_pointSourcePDF.integral(umax);
	     };

	private:
	    extendedSources::iSource& m_source;
	    extendedSources::iSource& m_source2;
	    ExtendedSourcePseudoPSF m_psf;
	    ExtendedSourcePseudoPSF m_psf2;
            double m_gamma,m_gamma2;
	    double m_sigma,m_sigma2;
	    double m_s2ratio;
	    double m_frac2;
	    double m_umax;
            double m_norm;
            skymaps::PsfFunction m_pointSourcePDF;
            skymaps::PsfFunction m_pointSourcePDF2;
            bool m_isPointSource;

    };



};

#endif

