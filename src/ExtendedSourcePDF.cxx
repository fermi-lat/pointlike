#include "pointlike/ExtendedSourcePDF.h"

#include "st_facilities/GaussianQuadrature.h"

#include <stdexcept>
#include <map>

    //--------------------------------------------------------------------------------------------------
    //  ExtendedSourcePDF
    //   
    //  Semi-analytic calculation of the PDF for an extended source given the GLAST PSF parametrization 
    //  v corresponds to the radial coordinate of the source extension in scaled space v= 0.5 * r**2 / sigma**2
    //  u is the deviation from the source center u= 0.5 * r**2 / sigma**2
    //  the integral over v is performed numerically with dgaus8
    //  the integral over azimuth angles was done analytically
    //--------------------------------------------------------------------------------------------------

#undef TEST_GRADIENT

namespace pointlike{

      static extendedSources::diskSource ExtendedSourcePseudoPSF2::defaultSource;

     static int extSource_instance_id;
     
     extendedSourcePDF::extendedSourcePDF(extendedSources::iSource& source, double g, bool usecache):
          m_source(source),
          m_gamma(g),
	  m_prefactor(2*M_PI*(g-1.)/g),
	  m_a(0.5*g),
	  m_b(0.5+0.5*g),
	  m_u(-1.),
          m_gamma_pow_gamma(pow(g,g)),
          m_hF(m_a,m_b,1.),
          m_dHF(m_a+1.,m_b+1.,2.),
	  m_useCaching(usecache),
	  m_cacheIpolTolerance(2.e-2),
	  cached_srcParam(0),
	  cacheID(extSource_instance_id++)
         {
             if(m_useCaching) {
		   cached_srcParam=m_source.get();
             };   

	 };   


    double extendedSourcePDF::valueFromIntegral(double u) const { 
          double error=0, result=0.;
	  int    ierr=-1;
	  m_u = u;
	  try{
	     result= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF>(*this,m_source.min(),m_source.max(),error,ierr);
          } catch (st_facilities::GaussianQuadrature::dgaus8Exception ex) { 
	     error=1e-5;
	     result= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF>(*this,m_source.min(),m_source.max(),error,ierr);
	  };  
	  if(ierr!=1) throw std::runtime_error("WARNING. Error in dgaus8 integration.");
	  return result;
       };
       
     std::vector<double> extendedSourcePDF::gradient(double u, int npar) const {
          std::vector<double> result(npar+1);
          double error=0.; int    ierr=-1;
	  m_u = u;
	  static int catch_counter;
	  
	  gradUFunctor gu(m_source,*this);
	  try{
 	     result[0]= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradUFunctor>(gu,m_source.min(),m_source.max(),error,ierr);
	  } catch (st_facilities::GaussianQuadrature::dgaus8Exception ex) { 
//	     std::cerr<<"gu catch counter:" <<(++catch_counter)<<std::endl;
	     error=1e-5;
	     result[0]= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradUFunctor>(gu,m_source.min(),m_source.max(),error,ierr);
	  };  
	  if(ierr!=1) throw std::runtime_error("WARNING. Error in dgaus8 integration.");
              	
#ifdef TEST_GRADIENT	       
	  static const double h=1.e-7;
	  m_u+=h;
    	  double test2 = st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF>(*this,m_source.min(),m_source.max(),error,ierr);
     	  m_u-=2*h;
	  double test1 = st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF>(*this,m_source.min(),m_source.max(),error,ierr);
          m_u+=h;
	  double ngrad=(test2-test1)/(2*h);          
          std::cout<<"GRADIENT TEST: gradu analytic="<<std::scientific<<std::setprecision(4)<<result[0]<<" numeric="<<ngrad<<" ratio="<<(ngrad/result[0])<<std::endl;
	  std::vector<double> oldp(npar);		
	  for (int i = 0; i<npar; i++) oldp[i]=m_source.get(i);
#endif
	  	
	  for (int i = 1; i<=npar; i++) {
	      gradSFunctor gs(m_source,*this,i-1);

	      try{
	         if (m_source.split(i-1)>m_source.min()){
 		   result[i]= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradSFunctor>(gs,m_source.min(),m_source.split(i-1),error,ierr);
 		   result[i]+=st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradSFunctor>(gs,m_source.split(i-1),m_source.max(),error,ierr);
 		 } else {
		   result[i]= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradSFunctor>(gs,m_source.min(),m_source.max(),error,ierr);
	         };
	      } catch (st_facilities::GaussianQuadrature::dgaus8Exception ex) { 
//	         std::cerr<<"gs catch counter:" <<(++catch_counter)<<std::endl;
		 error=1e-5;
	         if (m_source.split(i-1)>m_source.min()){
 		   result[i]= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradSFunctor>(gs,m_source.min(),m_source.split(i-1),error,ierr);
 		   result[i]+=st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradSFunctor>(gs,m_source.split(i-1),m_source.max(),error,ierr);
 		 } else {
		   result[i]= st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF::gradSFunctor>(gs,m_source.min(),m_source.max(),error,ierr);
	         };
	      };  
       	      if(m_source.parameterShiftsMax(i-1)) result[i]+= operator()(m_source.max())*m_source.get(i-1)/m_source.sigma2();
	      if(m_source.parameterShiftsMin(i-1)) result[i]-= operator()(m_source.min())*m_source.get(i-1)/m_source.sigma2();
	      if(ierr!=1) throw std::runtime_error("WARNING. Error in dgaus8 integration.");

#ifdef TEST_GRADIENT	       
       	      double test0 = operator()(m_source.max())*m_source.get(i-1)/m_source.sigma2();
 	      std::cout<<"GRADIENT TEST: parameter value="<<std::scientific<<std::setprecision(4)<<oldp[i-1]<<" u="<<m_u<<" gamma="<<m_gamma<<" f="<<test0<<std::endl;
              std::vector<double> testp(oldp);	
              testp[i-1]+=h;
              m_source.set(testp);
       	      double test2 = st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF>(*this,m_source.min(),m_source.max(),error,ierr);
              testp[i-1]-=2*h;
              m_source.set(testp);
       	      double test1 = st_facilities::GaussianQuadrature::dgaus8<const extendedSourcePDF>(*this,m_source.min(),m_source.max(),error,ierr);
	      m_source.set(oldp);
	      double ngrad=(test2-test1)/(2*h);
	      std::cout<<"GRADIENT TEST: parameter="<<(i-1)<<" analytic="<<result[i]<<" numeric="<<ngrad<<" ratio="<<(result[i]/ngrad)<<std::endl;

#endif	      
	  };   
	
	  return result;
     };
     
 
    double extendedSourcePDF::operator() (double v) const {
          return m_source(v)*gpdf(v);
    };

    double extendedSourcePDF::gpdf(double v) const {
	  double uvg=m_gamma+m_u+v;
          return m_prefactor*pow(uvg/m_gamma,-m_gamma)*m_hF(4*m_u*v/(uvg*uvg));
       };

//     double extendedSourcePDF::uvgpdf(double x) const {
//	  return m_prefactor*m_gamma_pow_gamma*pow(x,m_gamma-2.)*hypergeometric2F1(m_a,m_b,1.,4*m_u*x*(1.-x*(m_u+m_gamma)));
//       };


    double extendedSourcePDF::du_gpdf(double v) const {
	  double uvg = m_gamma+m_u+v;
	  double x   = 4*m_u*v/(uvg*uvg) ; 
	  double p0  = pow(uvg/m_gamma,-m_gamma) / (uvg*uvg*uvg) ;
	  double s1  = -m_gamma * uvg * uvg * m_hF(x);
	  double s2  = 4 * m_a * m_b * v * ( m_gamma-m_u+v ) *m_dHF(x);
	  return m_prefactor*p0*(s1+s2);
       };

    double extendedSourcePDF::valueFromCache(double u) const {
        static int cache_miss;
	static int cache_hit;
    
        cacheValidityCheck();
	std::map<double,double>::iterator psfit1 = psfCache.lower_bound(u);
	if (psfit1==psfCache.begin() || psfit1==psfCache.end()){
	   double psf=valueFromIntegral(u);
	   psfCache[u]=psf;
	   cache_miss++;
	   return psf;
	};
	std::map<double,double>::iterator psfit0 = psfit1; psfit0--;
	double val0=(*psfit0).second,val1=(*psfit1).second;
	if(fabs(val1-val0)>m_cacheIpolTolerance*val0){
	   double psf=valueFromIntegral(u);
	   psfCache[u]=psf;
	   cache_miss++;
	   return psf;
	};
	double u0=(*psfit0).first,u1=(*psfit1).first;
	double psf=val0+(val1-val0)/(u1-u0)*(u-u0);    
        cache_hit++;
/*	if(cache_hit%100000==0)
	   std::cout<<"id="<<cacheID<<" u="<<u<<" u0="<<u0<<" u1="<<u1<<" v1="<<val1<<" v0="<<val0
	            <<" psf="<<psf<<" true="<<valueFromIntegral(u)
	            <<" hit="<<cache_hit<<" miss="<<cache_miss<<" size="<<psfCache.size()
		    <<std::endl;
*/	return psf;
    };

    bool extendedSourcePDF::cacheValidityCheck() const{
       if(cached_srcParam!=m_source.get()){
//          std::cout<<"Source parameters changed. Cache "<<cacheID<<" invalid"<<std::endl;
          psfCache.clear();
	  cached_srcParam=m_source.get();
	  return false;
       };	
       return true;    
    };


    //--------------------------------------------------------------------------------------------------
    //  ExtendedSourcePseudoPsf
    //   
    //  Sell the PDF for an extended source to pointlike as an adjusted PSF
    //  This allows to keep most of the pointlike structures unchanged. Only minimal adjustments in the
    //  parameter handling are necessary 
    //--------------------------------------------------------------------------------------------------
    
    double ExtendedSourcePseudoPSF::operator () (const astro::SkyDir & r, const astro::SkyDir & r_prime, double sigma) 
    {
	double u( 0.5*(r_prime() - r()).mag2()/sigma/sigma);
	return operator()(u);
    };

    double ExtendedSourcePseudoPSF::operator()(double u) const {   
        if (m_isPointSource) return m_pointSourcePDF(u);
	return m_pdf.value(u); 
    };   

    std::vector<double> ExtendedSourcePseudoPSF::gradient(double u,int ncomp) const {   
        std::vector<double> grad(ncomp+1);
	if (m_isPointSource) grad[0]= - m_pointSourcePDF(u)/(1+u/m_gamma);
	else grad= m_pdf.gradient(u,ncomp);
        return grad;
    };   

    double ExtendedSourcePseudoPSF2::operator () (const astro::SkyDir & r, const astro::SkyDir & r_prime, double sigma) 
    {
	double u( 0.5*(r_prime() - r()).mag2()/sigma/sigma);
	return operator()(u);
    };

    double ExtendedSourcePseudoPSF2::operator()(double u) const {   
        double r1,r2=0;
        if (m_isPointSource) {
	  r1= m_pointSourcePDF(u); 
	  if (m_gamma2>0) r2= m_s2ratio*m_pointSourcePDF2(m_s2ratio*u); 
	} else {  
	  r1= pdf().value(u); 
	  if (m_gamma2>0) r2= m_s2ratio*pdf2().value(m_s2ratio*u); 
	};  
//	std::cout<<r1<<" "<<r2<<" "<<m_frac2<<" "<<u<<" "<<(m_s2ratio*u)<<std::endl;
        return (1-m_frac2)*r1+m_frac2*r2;
	  
	  
    };   
    

    std::vector<double> ExtendedSourcePseudoPSF2::gradient(double u,int ncomp) const {   
        std::vector<double> grad(ncomp+1);

        if (m_isPointSource) {
           double grad1,grad2=0;
	   grad1= - m_pointSourcePDF(u)/(1+u/m_gamma);
	   if (m_gamma2>0) grad2= - m_s2ratio*m_pointSourcePDF(m_s2ratio*u)/(1+m_s2ratio*u/m_gamma);
	   grad[0] = (1-m_frac2)*grad1+m_frac2*grad2;
	} else {
	   std::vector<double> grad2(ncomp+1);
           grad = m_psf.pdf().gradient(u,ncomp);
	   if (m_gamma2>0) {
	      grad2= m_psf2.pdf().gradient(m_s2ratio*u,ncomp);
	      for(int i=0; i<ncomp;i++) grad[i]=(1-m_frac2)*grad[i]+m_s2ratio*m_frac2*grad2[i];
	   };
	};    
        return grad;
    };   

};
