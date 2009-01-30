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

namespace {
   static const double GSL_DBL_EPSILON   =  2.2204460492503131e-16;
}

namespace pointlike {


class GammaFunction {
  public:  
    static double lngamma(double x);
    static double gamma(double x);
//  --------------- ratio of gamma functions (with possibly negative arguments)
    static double gammaRatio(double w, double z);
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
 
        void precalc(int cmax);

      public:
        HypergeometricSeries(double a,double b, double c, double err=GSL_DBL_EPSILON, seriesType stype=h2F1, int cmax=200);
        void set(double a,double b, double c);
	double operator () (double x) const;
    }; 


/////////////////////////////////////////////////////////////////////////////////////////////////

class HypergeometricFunction2F1 {
    private:
       double m_a,m_b,m_c,m_err;
       HypergeometricSeries  m_series;
       HypergeometricSeries  m_series_reflected1;
       HypergeometricSeries  m_series_reflected2;

    public:
       HypergeometricFunction2F1(double a,double b, double c, double err=GSL_DBL_EPSILON);
       double operator() (double x) const;
       void set(double a,double b, double c);
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

       HypergeometricFunction1F1(double a,double b, double err=GSL_DBL_EPSILON,double logf=0);
       double operator() (double x) const;
       void set(double a,double b);
       
       private:
                   
       bool asymp_eval_ok(double x,double a) const {
         return (x>100 && std::max(fabs(m_b-a),1.)*std::max(fabs(1.-a),1.)<0.7*fabs(x) );
       };

       double asymp(double x) const ;
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
       HypergeometricFunctionUbgt1(double a,double b, double err=GSL_DBL_EPSILON,double logf=0);
       ~HypergeometricFunctionUbgt1();
       double operator() (double x) const ;
       void set(double a,double b);
       
    private:
       bool seriesEvalOK(double a,double b,double x) const
	    { return ((fabs(a) < 5 && b < 5 && x < 2.0) || (fabs(a) <  10 && b < 10 && x < 1.0));};
       bool asympEvalOK(double a,double b,double x) const
	    { return (std::max(fabs(a),1.0)*std::max(fabs(1.0+a-b),1.0) < 0.99*fabs(x));};
       double series(double x) const ;
       double asymp(double x) const;	    
       double hyperg_U_CF1(const double a, const double b, const double x,const double N=0) const;
       double d9chu(const double a, const double b, const double x) const;

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

};
#endif // Hypergeometric_function_h
