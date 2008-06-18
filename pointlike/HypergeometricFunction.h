#ifndef Hypergeometric_function_h
#define Hypergeometric_function_h


#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

// lazy man's evaluation of the hypergeometric function
// only valid between 0<= x < 1 :
// for x<0.7 evaluate series, for x >=0.7 use the reflection formula 15.3.7 from Abramowitz and Stegun before
// evaluating the series to accelerate convergence.


namespace {

//  --------------- helper functions to evaluate the hypergeometric function

   static const double GSL_DBL_EPSILON   =  2.2204460492503131e-16;
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

    inline double lngamma(double x) {

      x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

      double Ag = lanczos_7_c[0];
      for(int k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

      /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
      double term1 = (x+0.5)*log((x+7.5)/M_E);
      double term2 = LogRootTwoPi + log(Ag);
      return term1 + (term2 - 7.0);
    }

//  --------------- ratio of gamma functions (with possibly negative arguments)

    inline double gammaRatio(double w, double z){
       if (w>0 && z>0) return exp(lngamma(w)-lngamma(z));

       double sw = sin(-M_PI*w);
       if(fabs(sw)<GSL_DBL_EPSILON && w<=0) throw std::invalid_argument("Trying to evaluate gamma function at negative n");
       double sz = sin(-M_PI*z);
       if(fabs(sz)<GSL_DBL_EPSILON && w<=0) throw std::invalid_argument("Trying to evaluate gamma function at negative n");

       if(w>0 && z<0) return exp(lngamma(w)+lngamma(1-z))*sz/M_PI;
       if(w<0 && z>0) return exp(lngamma(z)-lngamma(1-w))*M_PI/sw;
       if(w<0 && z<0) return exp(lngamma(1-z)-lngamma(1-w))*sz/sw;
       return NAN;
    };
};


    class HypergeometricSeries {
        double m_a,m_b,m_c;
        double m_err;
        std::vector<double> seriesCoeff;

       void precalc(){
          for(int k=0; k<200.; k+=1)
	     seriesCoeff[k] = (m_a+k)*(m_b+k) / ((m_c+k) * (k+1.0));  /* Gauss series */
        };

     public:
        HypergeometricSeries(double a,double b, double c, double err=GSL_DBL_EPSILON): m_a(a), m_b(b), m_c(c), m_err(err), seriesCoeff(200){
           precalc();
        }

        void set(double a,double b, double c){ 
            m_a=a; m_b=b; m_c=c;
            precalc();
        };

        double operator () (double x) const {
	  double sum_pos = 1.0;
	  double sum_neg = 0.0;
	  double del = 1.0;
          double delabs;

          std::vector<double>::const_iterator sCiter = seriesCoeff.begin();

          while(sCiter!=seriesCoeff.end()){
	     del *= (*sCiter++)*x;  /* Gauss series */
             if(del>0) { sum_pos+=del; delabs=del; }
             else      { sum_neg+=del; delabs=-del;};
	     if (delabs/(sum_pos+sum_neg)<m_err) return sum_pos+sum_neg;
          };
          throw std::runtime_error("ERROR in hypergeometric series: Too many iterations.");
          return sum_pos+sum_neg;
        };
    }; 


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

         double t1 = gammaRatio(m_c,m_c-m_a)*gammaRatio(m_c-m_a-m_b,m_c-m_b);
         double t2 = gammaRatio(m_c,m_a)*gammaRatio(m_a+m_b-m_c,m_b)*pow(1-x,m_c-m_a-m_b);
	 double h1 = m_series_reflected1(1.-x);
	 double h2 = m_series_reflected2(1.-x);
	 return t1*h1+t2*h2;
       };

       void set(double a,double b, double c){
             m_series.set(a,b,c);
             m_series_reflected1.set(a,b,a+b-c+1);
             m_series_reflected2.set(c-a,c-b,c-a-b+1.);
       };
    };
    


#endif // Hypergeometric_function_h
