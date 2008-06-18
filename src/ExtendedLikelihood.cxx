/** @file ExtendedLikelihood.cxx
    @brief Implementation of class ExtendedLikelihood

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/ExtendedLikelihood.cxx,v 1.24 2008/02/19 21:03:54 burnett Exp $
*/

#include "pointlike/ExtendedLikelihood.h"

#include "skymaps/DiffuseFunction.h"
#include "astro/SkyDir.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <TMinuit.h>

#undef TEST_GRADIENT

using astro::SkyDir;
using skymaps::SkySpectrum;
using skymaps::DiffuseFunction;

using namespace pointlike;

//#define DEBUG_PRINT
double ExtendedLikelihood::s_defaultUmax =200;
double  ExtendedLikelihood::s_tolerance(0.05); // default:

namespace {

    static ExtendedLikelihood* gExtendedPointer = NULL;
    void extended_likelihood_wrapper(Int_t &npar, Double_t* derivative, 
				     Double_t & loglike, Double_t par[], Int_t iflag){ 
      double alpha = par[0];
      if (alpha < 0) alpha = gExtendedPointer->alpha();
      loglike = gExtendedPointer->operator()(alpha);
      std::cout << "\tAlpha: " << alpha << " " << loglike << std::endl;

    }

    bool debug(false);
    double tolerance(0.05); // integral tolerance
    inline double sqr(double x){return x*x;}
    std::ostream * psf_data = &std::cout;

    // for the binning in 1/q when q is negative, pixels far from source
    double binsize(0.025); // bin size for 1/q. set zero to disable
    std::map<double,int> qmap; // used for binning 1/q


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class Convert  {  // Object to run over all pairs of HealPixel and photon number pairs 

    public:           // 
        Convert( const SkyDir& dir, ExtendedSourcePseudoPSF& f,
		 const astro::SkyFunction& background,
            double sigma, double umax,
            std::vector< std::pair<double, int> >& vec2,
            std::vector<int>& vec4,
             bool subset
            )
            : m_dir(dir), m_f(f), m_sigma(sigma)
            , m_back(background )
            , m_vec2(vec2)
            , m_vec4(vec4)
            , m_umax(umax)
            , m_F( f.integral(umax) ) // for normalization of the PSF
            , m_sum(0), m_count(0)
            , m_sumu(0) // average u, useful for calibration
            , m_sumb(0) // average b, useful for calibration
            , m_pixels(0)
            , m_subset(subset)
            , m_maxu_found(0)
        {
	  // SF: NO MORE NEEDED
//             double angle(sqrt(2.*umax)*sigma);
//             if(angle<0.1) { 
// 	        m_umax= 0.01/(2*sigma*sigma);
// //                std::cout<<"umax adjusted to "<<m_umax<<std::endl;
// 	    };	

//             if( ExtendedLikelihood::diffuse()!=0){
//                ExtendedLikelihood::diffuse()->setEnergyRange(emin, emax);
//                m_back_norm = ExtendedLikelihood::diffuse()->average(dir, angle, ExtendedLikelihood::tolerance());
//                if( m_back_norm==0){
//                    std::cerr << "Warning: normalization zero" << std::endl;
//                    m_back_norm=0.1; // kluge, like below
//                }
//             }
            m_first = m_vec4.size()==0;
	    
            if(debug){
                //psf_data = new std::ofstream("d:/users/burnett/temp/psf.txt");
                (*psf_data) << "u        f(u)      count    q" << std:: endl;
            }
        }
        ~Convert(){
             //if(debug){ psf_data->close();}
        }
        void operator()(const std::pair<SkyDir, int> x){
            double diff = x.first.difference(m_dir);
            double  u = sqr(diff/m_sigma)/2.;
            if(u>m_maxu_found) m_maxu_found = u;
            double signal(m_f(u)/m_F)
                 , bkg(m_back(x.first()))
                 , q( bkg/(signal*m_umax-bkg)); 

            if(isnan(signal) || isinf(signal))
	       std::cerr<<"WARNING: nan/inf encountered in signal ."<<std::endl;

            m_sum   += x.second*signal;
            m_count += x.second;
            m_sumu  += x.second*u;
            m_sumb  += bkg;
            m_pixels+= 1;
            if(debug){
                (*psf_data) << std::left<< std::setw(12) 
                    << u << std::setw(12) 
                    << signal << std::setw(5)
                    <<  x.second 
                    << std::setw(10)<<  bkg << std::endl;
            }
            if( binsize==0. || q>0){
                // not binning: save all in the array
                m_vec2.push_back(std::make_pair(q, x.second) );
            }else{
                // bin (in 1/q) the pixels with large u, on periphery
                // discretize these values of q in bins of 1/q
                double qbin = 1./binsize/floor(1./binsize/q+0.5);
                qmap[qbin]+=x.second;
            }
        }
        
        void consolidate() {
            if( binsize==0) return;
            // recombine the qmap with m_vec2.
            for(std::map<double,int>::iterator it = qmap.begin();it!=qmap.end();++it) {
                m_vec2.push_back(std::make_pair(it->first,it->second));
            }
            qmap.clear();
        }
        double average_f()const {return m_count>0? m_sum/m_count : -1.;}
        double average_u()const {return m_count>0? m_sumu/m_count : -1;}
        double average_b()const {return m_count>0? m_sumb/m_pixels: -1;}
        double count()const{return m_count;}

    private:
        SkyDir m_dir;

        ExtendedSourcePseudoPSF& m_f;
        const astro::SkyFunction& m_back;
        double m_sigma;
        std::vector<std::pair<double, int> >& m_vec2;
        std::vector<int>& m_vec4;
        
	std::vector<double> sourceTemplate;
       
        double m_umax, m_F;
        double m_sum, m_count, m_sumu, m_sumb, m_pixels;
        double m_back_norm;
        bool m_first;   //first call?
        bool m_subset;
        double m_maxu_found;
    };

     //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // The log likelihood function
    class LogLike  {
    public:
        LogLike(double a):m_a(a){}

        double operator()(double prev, const std::pair<float, int>& x)
        {
            return prev - x.second * log(m_a/x.first+1.);
        }
        double m_a;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // for accumulating the derivatives of the likelihood
    class Derivatives  {
    public:

        Derivatives(double x): m_x(x){}

        std::pair<double, double> operator()(std::pair<double, double> prev, const std::pair<float, int>& v)
        {
            double t(m_x+v.first)
                , d1(v.second/t)
                , d2(d1/t) ; // first, second derivatives (except sign)

            return std::make_pair(prev.first- d1, prev.second + d2);
        }
        double m_x;
    };
} // anon namespace
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// private nested class to manage normalized background

class ExtendedLikelihood::NormalizedBackground : public astro::SkyFunction {
  
public:
  NormalizedBackground(const skymaps::SkySpectrum* diffuse, const SkyDir& dir, double angle, double emin, double emax)
    : m_diffuse(diffuse)
      , m_emin(emin)
      , m_emax(emax){
    if( m_diffuse!=0){
      m_diffuse->setEnergyRange(emin, emax);
      m_back_norm = m_diffuse->average(dir, angle, ExtendedLikelihood::tolerance());
      if( m_back_norm==0){
	std::cerr << "Warning: normalization zero" << std::endl;
	m_back_norm=0.1; // kluge, like below
      }
    }
  }
    //! the normalized (in 0<u<umax)  background in the direction dir
  double operator()(const SkyDir& dir)const{
    double val(1.);
    if( m_diffuse!=0){
      m_diffuse->setEnergyRange(m_emin, m_emax);
      val=(*m_diffuse)(dir)/m_back_norm;
    }
    // a return value in the given direction
    if( val==0){ val=0.1;} // prevent zero, which is bad
    return val;
  }
  
private:
  const skymaps::SkySpectrum* m_diffuse;
    double m_back_norm;
    double m_emin, m_emax;
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ExtendedLikelihood::ExtendedLikelihood(const skymaps::Band& band,
				       const astro::SkyDir& dir,
				       const std::string& name, 
				       const std::vector<double>& src_param,
				       double umax, 
				       const skymaps::SkySpectrum* diffuse)
        : m_band(band)
        , m_averageF(0)
  	, m_psf(extendedSources::generator(name,band.sigma()),band.gamma())
	, m_src_param(src_param)
        , m_sigma(band.sigma())
        , m_alpha(-1)
        , m_curv(-1)
        , m_background(-1)  // default: no estimate
        , m_umax(umax)
        , m_emin(band.emin()), m_emax(band.emax())
        , m_diffuse(diffuse)
	, m_gradient(src_param.size()+2)
        , m_vloglike(0)
        , m_vjacobian(src_param.size()+2)
        , m_back(0)
{ 

    m_psf.source().set(m_src_param);
    setDir(dir);

    m_fint = m_psf.integral(m_umax);// integral out to umax

    // the integral of square, used by the quick estimator
    m_fint2 = m_psf.integralSquare(m_umax);
}

// legacy constructor

ExtendedLikelihood::ExtendedLikelihood(const skymaps::Band& band,
				       const astro::SkyDir& dir,
				       double umax, 
				       const skymaps::SkySpectrum* diffuse)
        : m_band(band)
        , m_averageF(0)
  	, m_psf(extendedSources::generator("point",band.sigma()),band.gamma())
        , m_sigma(band.sigma())
        , m_alpha(-1)
        , m_curv(-1)
        , m_background(-1)  // default: no estimate
        , m_umax(umax)
        , m_emin(band.emin()), m_emax(band.emax())
        , m_diffuse(diffuse)
	, m_gradient(2)
        , m_vloglike(0)
        , m_vjacobian(2)
        , m_back(0)
{ 

    setDir(dir);
    m_fint = m_psf.integral(m_umax);// integral out to umax
    // the integral of square, used by the quick estimator
    m_fint2 = m_psf.integralSquare(m_umax);
}

ExtendedLikelihood::~ExtendedLikelihood()
{
    delete m_back;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void ExtendedLikelihood::setDir(const astro::SkyDir& dir, bool subset)
{
    m_dir = dir;
    reload(subset);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void ExtendedLikelihood::reload(bool subset)
{
    // create set of ( 1/(f(u)/Fbar-1), weight) pairs in  m_vec2
    m_vec.clear();
    m_vec2.clear();
    delete m_back;
    double angle(sqrt(2.*m_umax)*sigma());
    if(angle<0.1) { 
      m_umax= 0.01/(2*sigma()*sigma());
      //                std::cout<<"umax adjusted to "<<m_umax<<std::endl;
    };	

    m_back = new NormalizedBackground(m_diffuse, m_dir, angle, m_emin, m_emax);
    m_band.query_disk(m_dir, angle, m_vec);
    Convert conv(m_dir, m_psf, *m_back, sigma(), m_umax, m_vec2, m_vec4, subset);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);
    result.consolidate();

    m_photon_count = static_cast<int>(result.count());

    if( m_photon_count==0) return; //throw std::invalid_argument("ExtendedLikelihood: no data after transform");

    m_averageF = result.average_f();
    m_avu = result.average_u();
    m_avb = result.average_b();

    // initialize to estimate.
    if( m_alpha<0) m_alpha = estimate();

}

void ExtendedLikelihood::setDir(const astro::SkyDir& dir, const std::vector<double>& src_param, bool subset){
    m_psf.source().set(src_param);
    setDir(dir,subset);								
};				

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double ExtendedLikelihood::operator()( double a) const
{
    if( a<0) a=m_alpha;
    return std::accumulate(m_vec2.begin(), m_vec2.end(), poissonLikelihood(a), LogLike(a));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ExtendedLikelihood::estimate() const
{
    if( m_photon_count==0) return -1;
    double est((m_umax*m_averageF-m_fint)/(m_umax*m_fint2-m_fint));
    if (est<0 ) est = 1e-4;
    if (est>1 ) est = 1-1e-4;
    return est;
}

/*
std::pair<double,double> ExtendedLikelihood::maximize(bool newStyle){

  if (newStyle){
     std::cout << "Using the new style" << std::endl;
    return maximizeMinuit();
  } else{ 
     std::cout << "Using the old style" << std::endl;
    return maximize();
  }
}
*/

std::pair<double,double> ExtendedLikelihood::maximizeMinuit()
{
  static int itermax(5000);
  if( m_photon_count==0) return std::make_pair(0.,0.);
  double x(m_alpha); // starting point

  int npar = 1;
  TMinuit* gMinuit = new TMinuit(npar);
//   std::cout << "Setting likelihood for : " << this << std::endl;
  
    // Set the pointer for access in the minuit function
  gExtendedPointer = this;
  gMinuit->SetFCN(extended_likelihood_wrapper);

  double par[npar];
  double stepSize[npar];
  double minVal[npar];
  double maxVal[npar];
  std::string parName[npar];
  par[0] = x;
  stepSize[0] = 0.001;
  minVal[0] = 0.;
  maxVal[0] = 1.0;
  parName[0] = std::string("alpha");

  for (int i = 0; i < npar; ++i)
    gMinuit->DefineParameter(i, parName[i].c_str(), par[i], stepSize[i], minVal[i], maxVal[i]);

  Int_t ierflag = 0;
  Double_t arglist[10];
  arglist[0] = 0.5; 
  int nargs = 1;
  gMinuit->mnexcm("SET ERR", arglist, nargs, ierflag);
  arglist[0] = itermax; 
  arglist[1] = 0.01; 
  nargs = 2;
  gMinuit->mnexcm("MIGRAD", arglist, nargs, ierflag);

  m_alpha = x;
  double errX;
  gMinuit->GetParameter(0, x, errX);
  delete gMinuit;
  return std::make_pair(x, errX);

}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::pair<double,double> ExtendedLikelihood::maximize()
{
    static double tol(1e-3), tolchange(0.05);
    static int itermax(5);
    if( m_photon_count==0) return std::make_pair(0.,0.);
    double x(m_alpha); // starting point
    int iteration(0);
    std::pair<double, double> dw(100,100);

    double delta(0);
    for(; fabs(dw.first) > tol && iteration<itermax; ++iteration, x-= delta) {

        // accumulate first and second derivatives
      dw = accumulate(m_vec2.begin(), m_vec2.end(), poissonDerivatives(x), Derivatives(x));
      if( debug) std::cout 
	<< std::setw(5)  << iteration 
	<< std::setw(10) << std::setprecision(3) << x 
	<< std::setw(15) << dw.first 
	<< std::setw(15) << dw.second 
	<< std::setw(10) << this->operator()(x) 
	<< std::setw(15) << TS() << std::endl;
      
      
      if( fabs(dw.second) < 1e-6 ) {
	dw.second = 1.0;
	break;
      }
      delta =dw.first/dw.second;
      //        std::cout<<std::setprecision(13)<<delta<<std::endl;

        // prevent from going negative or greater than 1
      if( delta > x ) delta = x/2.; 
      if( x-delta>1 ) delta = (x-1.)/2.;
      
      // quit if change small compared with resolution
      double sigdel(1./sqrt(dw.second));
      if( fabs(delta/sigdel) < tolchange ){x-=delta; break;}
    }
//    std::cout<<std::setprecision(13)<<x<<std::endl;
//    std::cout<<iteration<<std::endl;
    if(isnan(x)) {
      std::cout<<std::setprecision(5)<<"WARNING: isnan encountered in maximize: x="
	       <<x<<" delta="<<delta<<std::endl;
      x=m_alpha;
    }
    m_alpha = x;
    m_sigma_alpha = 1./sqrt(dw.second);
    return std::make_pair(x, m_sigma_alpha);
    
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ExtendedLikelihood::poissonLikelihood(double a)const
{
    if( m_background<0) return 0; 
    double expect(signal(a)+background());
    if( expect<=0) return 0; //?
    return expect - m_photon_count*log(expect);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::pair<double,double> ExtendedLikelihood::poissonDerivatives(double a)
{
    double d1(0), d2(0);

    if( m_background>=0 ) {
        double t(a + background()/m_photon_count);
        d1 = m_photon_count*(1-1/t);  // first log likelihood derivative
        d2 = m_photon_count/t/t;      // second
    }
    // none specified 
    return std::make_pair(d1, d2);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

const std::vector<double>& ExtendedLikelihood::gradient(const Hep3Vector& ex,
							const Hep3Vector& ey) const {

#ifdef TEST_GRADIENT
    static const double h=1.e-7;

    const SkyDir old_dir = m_dir;
    setDir(SkyDir(old_dir()+h*ex));
    double l2x = operator()(m_alpha);    
    setDir(SkyDir(old_dir()-h*ex));
    double l1x = operator()(m_alpha);    
    setDir(SkyDir(old_dir()+h*ey));
    double l2y = operator()(m_alpha);    
    setDir(SkyDir(old_dir()-h*ey));
    double l1y = operator()(m_alpha);    
    setDir(old_dir());
 
    double gradx=(l2x-l1x)/(2*h);
    double grady=(l2y-l1y)/(2*h);
    
    std::cout<<"TEST GRADIENT: sky dir="<<old_dir()<<" numeric_x="<<gradx<<" numeric_y="<<grady<<std::endl;
#endif

    double sig2( sqr(m_sigma) );
 
    for(unsigned int k=0;k<m_src_param.size()+2;k++) m_gradient[k]=0.;

    if(m_vec2.size()==0) return m_gradient;
    PixelList::const_iterator it = m_vec.begin();
    std::vector<std::pair<double,int> >::const_iterator it2 = m_vec2.begin();
//     std::vector<double>::const_iterator it3 = m_vec3.begin();
    
    m_vloglike.clear();
    for(unsigned int k=0;k<m_src_param.size()+2;k++) m_vjacobian[k].clear();

    for( ; it< m_vec.end(); ++it){
        const std::pair<SkyDir,int>& h = *it;
        SkyDir d( h.first );
        int nphoton( h.second);
        double q = (*it2).first;

        Hep3Vector delta( m_dir()-d() ); 
        double u = 0.5*delta.mag2()/sig2 ;
	double x = delta.dot(ex);
	double y = delta.dot(ey);
//	std::cout<<std::setprecision(4)<<"delta="<<delta<<" x="<<x<<" y="<<y<<" frac="<<(sqrt(x*x+y*y)/delta.mag())<<std::endl; 
	
	
// #if 0 // original version
//         if(u>m_umax) continue;
// #else // Marshal mod
//         if((u>m_umax&&m_vec4.size()==0)||(find(m_vec4.begin(),m_vec4.end(),h.first.index())==m_vec4.end()&&m_vec4.size()!=0)) continue;
// #endif

//         if(fabs(u-(*it3))>1e-8 || nphoton!=(*it2).second) {
        if(nphoton!=(*it2).second) {
            ++it; continue;
            if (it==m_vec.end()) { 
               std::cerr<<"WARNING m_vec.size="<<m_vec.size()<<" m_vec2.size="<<m_vec2.size()<<" nph="<<nphoton<<"/"<<(*it2).second<<std::endl;
               throw std::runtime_error("gradient: length of pixel and q vectors don't agree.");
            };
        };

 //       double fhat = m_psf(u)*m_umax/m_fint;
 //       double qq    = 1./(fhat-1) ;
 //       double Aold    = (1+qq)/(m_alpha+qq) ;
  
        double psf=m_psf(u);
        double A= psf>0 ? -m_alpha * (1+q)/(m_alpha+q) / psf : 0;
        std::vector<double> grad = m_psf.gradient(u,m_src_param.size());

#ifdef USE_MINPACK 
        double sa= 2.*M_PI*sqr(m_sigma);
        double likelihood= psf/(sa* m_fint) * (m_alpha+q) /(1+q) ; // how else can one reconstruct the likelihood from q ?
        std::cout<<q<<" "<<m_alpha<<" "<<sa<<" "<<likelihood<<std::endl;
        double sqrtlike = (likelihood >0) ? sqrt( - nphoton * log(likelihood) ) : 0 ;
        std::cout<<q<<" "<<m_alpha<<" "<<sa<<" "<<likelihood<<" "<<sqrtlike<<std::endl;
        m_vloglike.push_back(sqrtlike); 
        if (sqrtlike==0) sqrtlike = 1.e305;
        m_vjacobian[0].push_back(0.5/sqrtlike*nphoton*A*grad[0]*x/sig2);
        m_vjacobian[1].push_back(0.5/sqrtlike*nphoton*A*grad[0]*y/sig2);
	for(unsigned int k=2;k<m_gradient.size();k++) m_vjacobian[k].push_back(0.5/sqrtlike*nphoton*A*grad[k-1]);
#endif

	m_gradient[0]+= nphoton*A*grad[0]*x/sig2;
	m_gradient[1]+= nphoton*A*grad[0]*y/sig2;
	for(unsigned int k=2;k<m_gradient.size();k++) m_gradient[k] += nphoton*A*grad[k-1];
//        std::cout<<std::setprecision(4)<<std::scientific<<"             : "<<" A="<<A<<" nphoton="<<nphoton<<" psf="<<m_psf(u)<<std::endl;

       ++it2;
    }

#ifdef TEST_GRADIENT
    std::cout<<"TEST GRADIENT: sky dir="<<old_dir()<<" ana_x="<<m_gradient[0]<<"("
             <<(m_gradient[0]/gradx)<<") ana_y="<<m_gradient[1]<<"("<<(m_gradient[1]/grady)<<")"<<std::endl;
      for(unsigned int k=0;k<m_gradient.size();k++) std::cout<<"gradient("<<k<<")="<< m_gradient[k]<<std::endl;
 //   m_gradient[0]=gradx;
 //    m_gradient[1]=grady;
#endif


   return m_gradient; 
};
 

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hep3Vector ExtendedLikelihood::ps_gradient() const
{
    Hep3Vector grad;
    m_curv=m_w=0;

    double sig2( sqr(sigma()) );
    double gamma( m_psf.gamma() );
    //double w(0); // -log likelihood (check)
    int count(0);

    Hep3Vector perp(m_dir().orthogonal());
    PixelList::const_iterator it = m_vec.begin();

    for( ; it< m_vec.end(); ++it){
        const std::pair<SkyDir,int>& h = *it;
        SkyDir d( h.first );
        int nphoton( h.second);
        Hep3Vector delta( m_dir() - d() ); 
        double u( 0.5*delta.mag2()/sig2);
        double y = perp*delta; // pick an arbitrary direction for the curvature

        double fhat( m_psf(u)*m_umax/m_fint)
            ,  q( 1./(fhat-1) )
            ,  A( (1+q)/(m_alpha+q) )
            ,  B( 1./(1.+u/gamma) );
        grad   += nphoton * delta*A*B;
        m_curv += nphoton * A*B*(1-y*y/sig2*B*(q*(1-m_alpha)/(m_alpha+q) + 1/gamma));
        m_w    -= nphoton * log(m_alpha/q+1);
        count  += nphoton;

    }
    grad -= m_dir()*(m_dir()*grad); //subtract component in direction of dir
    m_curv *=  m_alpha/sig2;   // save curvature for access
    m_photon_count = count;
    return grad*m_alpha/sig2;

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ExtendedLikelihood::ps_curvature() const
{
    if( m_curv==-1.) ps_gradient(); // in not initialized. (negative curvature is bad)
    return m_curv;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ExtendedLikelihood::TS(double alpha)const {
    return 2.*((*this)(0) - (*this)(alpha));
    //return -2*(*this)(alpha);
}
double ExtendedLikelihood::solidAngle()const{
    return 2*M_PI* sqr(sigma())*m_umax;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ExtendedLikelihood::operator()(const astro::SkyDir& dir)const
{
    double diff =dir.difference(m_dir); 
    double  u = sqr(diff/sigma())/2.;
    // just to see what is there
    // astro::SkyDir r(x.first()); double ra(r.ra()), dec(r.dec());
    double jacobian (2.*M_PI * sqr(sigma()) );
    return signal()*m_psf(u)/ m_fint/jacobian;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double ExtendedLikelihood::display(const astro::SkyDir& dir, int mode) const
{
    double diff( dir.difference(m_dir) ); 
    double u   ( sqr(diff/sigma())/2.  );
    double jacobian( 2.*M_PI* sqr(sigma()) );
    if( mode ==0){
        // default is density, as above
        return signal()*m_psf(u)/ m_fint/jacobian;
    }
    SkyDir pdir(m_band.dir(m_band.index(dir)));
    u = sqr( pdir.difference(m_dir)/sigma())/2.; // recalculate u for pixel
    if( u> m_umax ) return 0;

    // now look for it in the data
    PixelList::const_iterator it;
    int index( m_band.index(pdir));
    for( it=m_vec.begin(); it!=m_vec.end(); ++it ) {
        if( m_band.index(it->first) == index) break;
    }

    int counts (it==m_vec.end()? 0:  it->second);

    if( mode==1 ) return counts;   // observed in the pixel
    double back( (*m_back)(pdir) );           // normalized background at the pixel: average should be 1.0
    if( mode==2 ) return back; 
    double area( m_band.pixelArea() );

    double prediction ( area/jacobian * photons() * ( m_alpha*m_psf(u)/m_fint + (1.-m_alpha)*back/m_umax ) );
    if( mode==3 ) return prediction; // prediction of fit
    if( mode==4 ) return counts-prediction;
    return 0; 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


const skymaps::SkySpectrum* ExtendedLikelihood::diffuse() const
{
    return m_diffuse;
}

void ExtendedLikelihood::setDiffuse(skymaps::SkySpectrum* diff)
{
    m_diffuse = diff;
    setDir(m_dir); // reset data
}

double ExtendedLikelihood::tolerance()
{
    return s_tolerance;
}

double ExtendedLikelihood::defaultUmax()
{
    return s_defaultUmax;
}
void ExtendedLikelihood::setDefaultUmax(double umax)
{
    s_defaultUmax = umax;
}

void ExtendedLikelihood::setTolerance(double tol)
{
    s_tolerance = tol;
}

void ExtendedLikelihood::recalc(bool subset) 
{
    // create set of ( 1/(f(u)/Fbar-1), weight) pairs in  m_vec2
    m_vec2.clear();
    Convert conv(m_dir, m_psf, *m_back, sigma(), m_umax, m_vec2, m_vec4, subset);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);
    result.consolidate();

    if(m_photon_count!= static_cast<int>(result.count())) {
        m_photon_count=m_photon_count;
    } 

    if( m_photon_count==0) return; //throw std::invalid_argument("ExtendedLikelihood: no data after transform");

    m_averageF = result.average_f();
    m_avu = result.average_u();
    m_avb = result.average_b();
}
