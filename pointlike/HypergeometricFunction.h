#ifndef Hypergeometric_function_h
#define Hypergeometric_function_h


#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <cmath>

#if defined(WIN32) && !defined(NAN)
static const unsigned int nan[2]={0xffffffff, 0x7fffffff};
#define NAN (*(const double *) nan)
#endif



// lazy man's evaluation of the hypergeometric function
// only valid between 0<= x < 1 :
// for x<0.7 evaluate series, for x >=0.7 use the reflection formula 15.3.7 from Abramowitz and Stegun before
// evaluating the series to accelerate convergence.
 

namespace {

#ifdef WIN32 // quiet this for Seattle folk
    int isinf(double){return 0;}
    double gamma(double){return 0;}
 #include <float.h>
#endif
//  --------------- helper functions to evaluate the hypergeometric function

   static const double GSL_DBL_EPSILON   =  2.2204460492503131e-16;
   static const double GSL_SQRT_DBL_MIN  =  1.4916681462400413e-154;
   static const double GSL_SQRT_DBL_MAX  =  1.34e154;
   static const double GSL_LOGSQRT_DBL_MAX = log(GSL_SQRT_DBL_MAX);

   static const double LogRootTwoPi =  0.9189385332046727418;
   
   
   static const double lanczos_7_c[9] = {
	 0.99999999999980993227684700473478,
	 676.520368121885098567009190444019,
	-1259.13921672240287047156078755283,
	 771.3234287776530788486528258894,
	-176.61502916214059906584551354,
	 12.507343278686904814458936853,
	-0.13857109526572011689554707,
	 9.984369578019570859563e-6,
	 1.50563273514931155834e-7
      };


//  --------------- lngamma function from GSL, lanczos algorithm

}

class GammaFunction {
  public:  
    static double lngamma(double x) {
    

      x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

      double Ag = lanczos_7_c[0];
      for(int k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

      /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
      double term1 = (x+0.5)*log((x+7.5)/M_E);
      double term2 = LogRootTwoPi + log(Ag);
      return term1 + (term2 - 7.0);
    }
    
    static double gamma(double x){
       
       if(x==0) return -log(0.);
       else if(x>0) return exp(lngamma(x));
       else return M_PI/sin(M_PI*x)*exp(lngamma(1-x));
    };   
       

//  --------------- ratio of gamma functions (with possibly negative arguments)

    static double gammaRatio(double w, double z){
       if (w>0 && z>0) return exp(lngamma(w)-lngamma(z));

       double sw = sin(-M_PI*w);
       if(fabs(sw)<GSL_DBL_EPSILON && w<=0) throw std::invalid_argument("Trying to evaluate gamma function at negative n");
       double sz = sin(-M_PI*z);
       if(fabs(sz)<GSL_DBL_EPSILON && w<=0) throw std::invalid_argument("Trying to evaluate gamma function at negative n");

       if(w>0 && z<0) return exp(lngamma(w)+lngamma(1-z))*sz/M_PI;
       if(w<0 && z>0) return exp(lngamma(z)-lngamma(1-w))*M_PI/sw;
       if(w<0 && z<0) return exp(lngamma(1-z)-lngamma(1-w))*sz/sw;
#ifndef WIN32
       return NAN;
#else
       return 1e30; // can't figure out how to set it, should not matter
#endif
    };
};


/// Hypergeometric series to max 2,1 order. To get 2,0 or 1,1 just set c=0 or b=0 

    class HypergeometricSeries {
      public:
	enum seriesType {UNDEFINED,h2F1,h1F1,h2F0,h1F0} ;
      private:        
	double m_a,m_b,m_c;
        double m_err;
	seriesType m_stype;
	int m_cmax;
	
        std::vector<double> seriesCoeff;

       void precalc(int cmax){
           switch(m_stype) {
 	      case h1F1: 
	         for(int k=0; k< cmax; k+=1) seriesCoeff[k] = (m_a+k) / ((m_c+k) * (k+1.0)); break; /* Gauss series */
	      case h2F0:
	         for(int k=0; k< cmax; k+=1) seriesCoeff[k] = (m_a+k)*(m_b+k) / (k+1.0); break; /* Gauss series */
	      case h1F0: 
	         for(int k=0; k< cmax; k+=1) seriesCoeff[k] = (m_a+k) / (k+1.0); break; /* Gauss series */
	      case h2F1:
   	         for(int k=0; k< cmax; k+=1) seriesCoeff[k] = (m_a+k)*(m_b+k) / ((m_c+k) * (k+1.0)); break;  /* Gauss series */
	      default:
	         throw std::runtime_error("Error in HypergeometricSeries: Series type not supported."); 
	   };
	   	  
        };

     public:
        HypergeometricSeries(double a,double b, double c, double err=GSL_DBL_EPSILON, seriesType stype=h2F1, int cmax=200)
	   : m_a(a), m_b(b), m_c(c), m_err(err), m_stype(stype), m_cmax(cmax), seriesCoeff(cmax){
           precalc(cmax);
        }

        void set(double a,double b, double c){ 
            m_a=a; m_b=b; m_c=c;
            precalc(m_cmax);
        };

        double operator () (double x) const {
	  double sum_pos = 1.0;
	  double sum_neg = 0.0;
	  double del = 1.0;
          double delabs = 1.0;

          std::vector<double>::const_iterator sCiter = seriesCoeff.begin();

          while(sCiter!=seriesCoeff.end()){
	     del *= (*sCiter++)*x;  /* Gauss series */
	     if(m_stype==h2F0) {
	         if(fabs(del)>delabs) return sum_pos+sum_neg;  //special case nonconverging series F20: return when coefficients
             };                                                // increase in size (same treatment as in GSL code)
	     if(del>0) { sum_pos+=del; delabs=del; }
             else      { sum_neg+=del; delabs=-del;};
	     if (delabs/(sum_pos+sum_neg)<m_err) return sum_pos+sum_neg;
          };
//          throw std::runtime_error("ERROR in hypergeometric series: Too many iterations.");
          std::cerr<<"WARNING in hypergeometric series: Iteration limit reached without convergence."<<std::endl;
          return sum_pos+sum_neg;
        };
    }; 


/////////////////////////////////////////////////////////////////////////////////////////////////

    class HypergeometricFunction2F1{
    private:
       double m_a,m_b,m_c,m_err;
       HypergeometricSeries  m_series;
       HypergeometricSeries  m_series_reflected1;
       HypergeometricSeries  m_series_reflected2;

    public:
       HypergeometricFunction2F1(double a,double b, double c, double err=GSL_DBL_EPSILON): m_a(a), m_b(b), m_c(c), m_err(err),
             m_series(a,b,c,err),
             m_series_reflected1(a,b,a+b-c+1,err),
             m_series_reflected2(c-a,c-b,c-a-b+1.,err)
             {};

       double operator() (double x) const {
         if(x < 0 || x>=1.) {
	    std::cerr<<"WARNING in hypergeometric series: Invalid parameters. "<<std::endl;
	    std::cerr<<m_a<<" "<<m_b<<" "<<m_c<<" "<<x<<std::endl;
	    throw std::invalid_argument("Invalid parameters");
            return 0;
         }
         if (x<0.70) return m_series(x);

         double t1 = GammaFunction::gammaRatio(m_c,m_c-m_a)*GammaFunction::gammaRatio(m_c-m_a-m_b,m_c-m_b);
         double t2 = GammaFunction::gammaRatio(m_c,m_a)*GammaFunction::gammaRatio(m_a+m_b-m_c,m_b)*pow(1-x,m_c-m_a-m_b);
	 double h1 = m_series_reflected1(1.-x);
	 double h2 = m_series_reflected2(1.-x);
	 return t1*h1+t2*h2;
       };

       void set(double a,double b, double c){
             m_series.set(a,b,c);
             m_series_reflected1.set(a,b,a+b-c+1);
             m_series_reflected2.set(c-a,c-b,c-a-b+1.);
	     m_a=a;  m_b=b;  m_c=c;
       };
    };


////////////////////////////////////////////////////////////////////////////////////
    
    class HypergeometricFunction1F1{
    private:
       double m_a,m_b,m_err,m_lf,m_a0;
       HypergeometricSeries  m_series1F1_a0;
       HypergeometricSeries  m_series1F1_a0p1;
       HypergeometricSeries  m_series2F0;
       HypergeometricSeries  m_series2F0_a0;
       HypergeometricSeries  m_series2F0_a0p1;
       

    public:
       HypergeometricFunction1F1(double a,double b, double err=GSL_DBL_EPSILON,double logf=0)
           : m_a(a), m_b(b), m_err(err), m_lf(logf),m_a0(a-int(a)),
             m_series1F1_a0(m_a0,0,b,err,HypergeometricSeries::h1F1,500),
             m_series1F1_a0p1(m_a0+1,0,b,err,HypergeometricSeries::h1F1,500),
             m_series2F0(b-a,1-a,0,err,HypergeometricSeries::h2F0,500),
             m_series2F0_a0(b-m_a0,1-m_a0,0,err,HypergeometricSeries::h2F0,500),
             m_series2F0_a0p1(b-m_a0-1,-m_a0,0,err,HypergeometricSeries::h2F0,500)
             {};

       double operator() (double x) const {
         if(x < 0 ) {
	    std::cerr<<"WARNING in hypergeometric series: Invalid parameters. "<<std::endl;
	    std::cerr<<m_a<<" "<<m_b<<" "<<" "<<x<<std::endl;
	    throw std::invalid_argument("Invalid parameters");
            return 0;
         }

         if(asymp_eval_ok(x,m_a)) return asymp(x);
	   double result;
           int nrescales=0;
	   if(fabs(m_a-m_a0)<GSL_DBL_EPSILON) result=exp(m_lf)*m_series1F1_a0(x);
	   else if(fabs(m_a-m_a0-1)<GSL_DBL_EPSILON) result=exp(m_lf)*m_series1F1_a0p1(x);
	   else {
	       double mm1,m0,mp1;
	       double logf;
	       int nrescales=0; 
	       
	       if(asymp_eval_ok(x,m_a0+1)){
	          logf=m_lf+x;
                  m0=exp(GammaFunction::lngamma(m_b)-GammaFunction::lngamma(m_a0+1)+log(x)*(m_a0+1-m_b))*m_series2F0_a0p1(1/x);
	          mm1=exp(GammaFunction::lngamma(m_b)-GammaFunction::lngamma(m_a0)+log(x)*(m_a0-m_b))*m_series2F0_a0(1/x);
//  	          std::cout<<m_a<<" "<<x<<" "<<m_lf<<" "<<m0<<" "<<mm1<<" "
//		           <<m_series2F0_a0(1/x)<<" "
//			   <<(m_lf+GammaFunction::lngamma(m_b)-GammaFunction::lngamma(m_a0)+log(x)*(m_a0-m_b)+x)<<std::endl;
	       } else {
                  logf=m_lf;
		  mm1=m_series1F1_a0(x);
	          m0=m_series1F1_a0p1(x);
               };
#ifndef WIN32
	       if(isinf(m0)||isnan(m0)) std::cout<<"NAN/INF in m0 of 1F1 series "<<m_a<<" "<<m_lf<<" "<<x<<" "<<m_a0<<std::endl;
	       if(isinf(mm1)||isnan(mm1)) std::cout<<"NAN/INF in mm1 of 1F1 series "<<m_a<<" "<<m_lf<<" "<<x<<" "<<m_a0<<std::endl;
#endif

	       for(double a=m_a0+1;a<m_a-GSL_DBL_EPSILON;a+=1.){
	          mp1=((m_b-a)*mm1+(2*a-m_b+x)*m0)/a;
		  mm1=m0;
		  m0=mp1;
		  if(m0>GSL_SQRT_DBL_MAX){m0/=GSL_SQRT_DBL_MAX; mm1/=GSL_SQRT_DBL_MAX; nrescales++;};
	       }; 
               result = exp(logf+nrescales*GSL_LOGSQRT_DBL_MAX)*m0;
	   }; 
#ifndef WIN32 
           if(isinf(result)||isnan(result)) std::cout<<"NAN/INF in 1F1 series "<<m_a<<" "<<m_lf<<" "<<x<<" "<<m_a0<<" "<<nrescales<<std::endl;
#endif
	   return result;
       };

       void set(double a,double b){
             m_a = a; m_b = b; m_a0=m_a-int(m_a);
             m_series1F1_a0.set(m_a0,0,b);
             m_series1F1_a0p1.set(m_a0+1,0,b);
             m_series2F0.set(b-a,1-a,0);
             m_series2F0_a0.set(b-m_a0,1-m_a0,0);
             m_series2F0_a0p1.set(b-m_a0-1,-m_a0,0);
       };
       
       private:
                   
       bool asymp_eval_ok(double x,double a) const {
         return (x>100 && std::max(fabs(m_b-a),1.)*std::max(fabs(1.-a),1.)<0.7*fabs(x) );
       };
       
       double asymp(double x) const {
//	 std::cout<<"asymp. "<<std::endl;
	 double exp_factor = exp(m_lf+GammaFunction::lngamma(m_b)-GammaFunction::lngamma(m_a)+log(x)*(m_a-m_b)+x) ;
	 double series     = m_series2F0(1/x);
	 double result     = exp_factor*series;
#ifndef WIN32
	 if(isinf(result)) std::cout<<"NAN in 1F1 asymp "<<m_a<<" "<<m_b<<" "<<x<<std::endl;
#endif
	 return result;
       };
    };
    
/////////////////////////////////////////////////////////////////////////////////

    
    class HypergeometricFunctionUbgt1{
    private:
       double m_a,m_b,m_err,m_lf;
       double m_f1,m_f2,m_sinb;
       HypergeometricFunction1F1  m_hyperg1F1_a;
       HypergeometricFunction1F1  m_hyperg1F1_b;
       double m_a0;
       int m_nIter;
       HypergeometricFunctionUbgt1* m_ua0;

    public:
       HypergeometricFunctionUbgt1(double a,double b, double err=GSL_DBL_EPSILON,double logf=0)
           : m_a(a), m_b(b), m_err(err), m_lf(logf),
             m_hyperg1F1_a(a,b,err),
             m_hyperg1F1_b(1+a-b,2-b,err)
             {
	       m_f1=-GammaFunction::lngamma(b)-GammaFunction::lngamma(1+a-b);
	       m_f2=-GammaFunction::lngamma(a)-GammaFunction::lngamma(2-b);
	       m_sinb=M_PI/sin(M_PI*b);
               m_nIter=int(m_a)-1;
	       m_a0=m_a-m_nIter;
	       m_ua0=0;
	       if(m_nIter>0) m_ua0 = new HypergeometricFunctionUbgt1(m_a0,m_b,m_err);
               if(b<1)throw std::runtime_error("Only b>1 supported in this class.");
	     };

         
       ~HypergeometricFunctionUbgt1(){ if(m_ua0) delete m_ua0; };

       
       double operator() (double x) const {
         if(x < 0 ) {
	    std::cerr<<"WARNING in hypergeometric series: Invalid parameters. "<<std::endl;
	    std::cerr<<m_a<<" "<<m_b<<" "<<x<<std::endl;
	    throw std::invalid_argument("Invalid parameters");
            return 0;
         }
   	 
         double result;
         if (asympEvalOK(m_a,m_b,x)) {
	    result=asymp(x); 
//	    std::cout<<"asymp OK"<<std::endl; 
//            if(isnan(result)) std::cout<<"NAN in asymp"<<m_a<<" "<<m_b<<" "<<x<<std::endl;
	    return result; 
	 }
         else if(seriesEvalOK(m_a,m_b,x)) { 
	     result=series(x); 
//	     std::cout<<"series OK"<<std::endl; 
//             if(isnan(result)) std::cout<<"NAN in series"<<m_a<<" "<<m_b<<" "<<x<<std::endl;
	     return result; 
	 }
	 else {
	   // this code seems to work, however I am fundamentally lacking understanding....
	 
//           std::cout<<"continuous fraction & iteration."<<std::endl; 
	     
	   const double  rescale_factor=1e150;
	   const double  log_rf=log(rescale_factor);
	    
	   double result;
	   double r1   = hyperg_U_CF1(m_a,m_b,x);
//	   std::cout<<"r1="<<r1<<std::endl;
	   double ua   =  GSL_SQRT_DBL_MIN;
	   double uap1 = r1/m_a * ua;
           double uam1,u0;

           int nrescales=0;
	   for(double ap=m_a; ap>m_a0+0.1; ap -= 1.0) {
             uam1 = -((m_b-2.0*ap-x)*ua + ap*(1.0+ap-m_b)*uap1);
             uap1 = ua;
             ua   = uam1;
	     if(ua> rescale_factor){ uap1/= rescale_factor; ua/= rescale_factor; nrescales++;};
	   }
           
          if (seriesEvalOK(m_a0,m_b,x) || asympEvalOK(m_a0,m_b,x)){
	      if(m_ua0) u0= m_ua0->operator()(x);
	      else throw std::runtime_error("This line should never be reached.");
	      
	      result= GSL_SQRT_DBL_MIN*u0/ua*exp(m_lf-nrescales*log_rf);
//	      std::cout<<" ua="<<ua<<" u0="<<u0<<" "<<nrescales<<" res="<<result<<std::endl;
//              if(isnan(result)) std::cout<<"NAN in iteration"<<m_a<<" "<<m_b<<" "<<x<<std::endl;
	      return result;
           } else {
	      throw std::runtime_error("HypergeometricU not supported for a,b,x chosen");
 	   };
         };
      };   

       void set(double a,double b){
             m_a=a; m_b=b;
             m_hyperg1F1_a.set(a,b);
             m_hyperg1F1_b.set(1+a-b,2-b);
             m_f1=-GammaFunction::lngamma(b)-GammaFunction::lngamma(1+a-b);
	     m_f2=-GammaFunction::lngamma(a)-GammaFunction::lngamma(2-b);
             m_sinb=M_PI/sin(M_PI*b);
             m_nIter=int(m_a)-1;
             m_a0=m_a-m_nIter;
	     if(m_ua0) delete m_ua0;
  	     m_ua0=0;
             if(m_nIter>0) m_ua0 = new HypergeometricFunctionUbgt1(m_a0,m_b,m_err);
       };
       
    private:
       bool seriesEvalOK(double a,double b,double x) const
	    { return ((fabs(a) < 5 && b < 5 && x < 2.0) || (fabs(a) <  10 && b < 10 && x < 1.0));};
       bool asympEvalOK(double a,double b,double x) const
	    { return (std::max(fabs(a),1.0)*std::max(fabs(1.0+a-b),1.0) < 0.99*fabs(x));};
       double series(double x) const {
            double a1,a2;
	    if(m_b>0 && 1.+m_a-m_b>0)       a1 = exp(m_lf+m_f1)*m_hyperg1F1_a(x);
	    else if(m_b==0|| 1+m_a-m_b==0)  a1 = 0;
	    else                            a1 = exp(m_lf)/gamma(m_b)/gamma(1+m_a-m_b)*m_hyperg1F1_a(x);

	    if(m_a>0 && 2.-m_b>0)           a2 = exp(m_lf+m_f2)*m_hyperg1F1_b(x);
	    else if(m_a==0|| 2.-m_b==0)     a2 = 0;
	    else                            a2 = exp(m_lf)/gamma(m_a)/gamma(2.-m_b)*m_hyperg1F1_b(x);

//	    std::cout<<"                 ->a1="<<a1<< " a2="<<a2<<" f1="<<m_f1<<" f2="<<m_f2<<" a="<<m_a<<" b="<<m_b<<std::endl;
	    return m_sinb * (a1 - pow(x,1.-m_b)*a2);
       };
       double asymp(double x) const{  
            return exp(m_lf-m_a*log(x))*d9chu(m_a,m_b,x); 
       };	    

       double hyperg_U_CF1(const double a, const double b, const double x,const double N=0) const{

	    const int maxiter = 1000000;
	    int n = 1;
	    double Anm2 = 1.0;
	    double Bnm2 = 0.0;
	    double Anm1 = 0.0;
	    double Bnm1 = 1.0;
	    double a1 = -(a + N);
	    double b1 =  (b - 2.0*a - x - 2.0*(N+1));
	    double An = b1*Anm1 + a1*Anm2;
	    double Bn = b1*Bnm1 + a1*Bnm2;
	    double an, bn;
	    double fn = An/Bn;
            double del;
	
	    while(n < maxiter) {
	      double old_fn;
	      n++;
	      Anm2 = Anm1;
	      Bnm2 = Bnm1;
	      Anm1 = An;
	      Bnm1 = Bn;
	      an = -(a + N + n - b)*(a + N + n - 1.0);
	      bn =  (b - 2.0*a - x - 2.0*(N+n));
	      An = bn*Anm1 + an*Anm2;
	      Bn = bn*Bnm1 + an*Bnm2;
	      
	      if(fabs(An)>GSL_SQRT_DBL_MAX || fabs(Bn)>GSL_SQRT_DBL_MAX){
	         An/=GSL_SQRT_DBL_MAX ; Bn/=GSL_SQRT_DBL_MAX;
	         Anm1/=GSL_SQRT_DBL_MAX ; Bnm1/=GSL_SQRT_DBL_MAX;
	         Anm2/=GSL_SQRT_DBL_MAX ; Bnm2/=GSL_SQRT_DBL_MAX;
	      };
	      
	      old_fn = fn;
	      fn = An/Bn;
	      del = old_fn/fn;
//              std::cout<<fabs(del - 1.0)<<" "<<An<<" "<<Bn<<std::endl;

	      if(fabs(del - 1.0) < 10.0*GSL_DBL_EPSILON) break;
	    }

	    if(n == maxiter) {
	       std::cerr<<"WARNING: continued fraction accuracy problem. max iterations reached: del="<<fabs(del - 1.0)<<std::endl;
               //throw std::runtime_error("HypergeometricU evaluation with continued fraction failed.");
	    } ;
	    return fn;
	  }

    double d9chu(const double a, const double b, const double x) const {
      const double EPS   = 8.0 * GSL_DBL_EPSILON; 
      const int maxiter = 500;
      double aa[4], bb[4];
      int i;

      double bp = 1.0 + a - b;
      double ab = a*bp;
      double ct2 = 2.0 * (x - ab);
      double sab = a + bp;

      double ct3 = sab + 1.0 + ab;
      double anbn = ct3 + sab + 3.0;
      double ct1 = 1.0 + 2.0*x/anbn;

      bb[0] = 1.0;
      aa[0] = 1.0;

      bb[1] = 1.0 + 2.0*x/ct3;
      aa[1] = 1.0 + ct2/ct3;

      bb[2] = 1.0 + 6.0*ct1*x/ct3;
      aa[2] = 1.0 + 6.0*ab/anbn + 3.0*ct1*ct2/ct3;

      for(i=4; i<maxiter; i++) {
	int j;
	double c2;
	double d1z;
	double g1, g2, g3;
	double x2i1 = 2*i - 3;
	ct1   = x2i1/(x2i1-2.0);
	anbn += x2i1 + sab;
	ct2   = (x2i1 - 1.0)/anbn;
	c2    = x2i1*ct2 - 1.0;
	d1z   = 2.0*x2i1*x/anbn;

	ct3 = sab*ct2;
	g1  = d1z + ct1*(c2+ct3);
	g2  = d1z - c2;
	g3  = ct1*(1.0 - ct3 - 2.0*ct2);

	bb[3] = g1*bb[2] + g2*bb[1] + g3*bb[0];
	aa[3] = g1*aa[2] + g2*aa[1] + g3*aa[0];

	if(fabs(aa[3]*bb[0]-aa[0]*bb[3]) < EPS*fabs(bb[3]*bb[0])) break;

	for(j=0; j<3; j++) {
	  aa[j] = aa[j+1];
	  bb[j] = bb[j+1];
	}
      }

      double result = aa[3]/bb[3];

      if(i == maxiter) {
	 throw std::runtime_error("Error in HypergeometricU rational approximation.");
      }
      
      return result;    
    };

};


class HypergeometricFunctionU {
  private:
     double m_a,m_b,m_err,m_lf;
     HypergeometricFunctionUbgt1 m_u;
     bool m_isReflected;

  public:
     HypergeometricFunctionU(double a, double b, double err=GSL_DBL_EPSILON,double logf=0)
        : m_a(a), m_b(b), m_err(err), m_lf(logf)
        , m_u( b>=1 ? a : 1+a-b , b>=1 ? b : 2.-b ,err,logf) 
	, m_isReflected(b<1)
     {};
     
     double operator()(double x) const{
        if(m_isReflected) return pow(x,1.-m_b)*m_u(x);
	return m_u(x);
     };
     
     void set(double a,double b){
       m_isReflected=(b<1);
       m_u.set(b>=1 ? a : 1+a-b , b>=1 ? b : 2.-b);
     };
     	
};

#endif // Hypergeometric_function_h
