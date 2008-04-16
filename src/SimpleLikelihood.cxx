/** @file SimpleLikelihood.cxx
@brief Implementation of class SimpleLikelihood

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SimpleLikelihood.cxx,v 1.29 2008/04/14 18:23:32 burnett Exp $
*/

#include "pointlike/SimpleLikelihood.h"
#include "skymaps/DiffuseFunction.h"
#include "astro/SkyDir.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <iomanip>
#include <fstream>

using astro::SkyDir;
using healpix::HealPixel;
using skymaps::SkySpectrum;
using skymaps::DiffuseFunction;
using skymaps::PsfFunction;

using namespace pointlike;

//#define DEBUG_PRINT
#define QBIN
double SimpleLikelihood::s_defaultUmax =50;

skymaps::SkySpectrum* SimpleLikelihood::s_diffuse(0);
double  SimpleLikelihood::s_tolerance(0.05); // default:

int SimpleLikelihood::s_display_mode(0); 
void SimpleLikelihood::setDisplayMode(int newmode){
    s_display_mode=newmode;}

namespace {

    bool debug(false);
    double tolerance(0.05); // integral tolerance
    inline double sqr(double x){return x*x;}
    std::ostream * psf_data = &std::cout;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class Convert  {
    public:
        Convert( const SkyDir& dir, PsfFunction& f, 
            const astro::SkyFunction& background,
            double sigma, double umax,
            std::vector< std::pair<double, int> >& vec2,
            std::vector<int>& vec4,
             bool subset
            )
            : m_dir(dir), m_f(f), m_sigma(sigma)
            , m_vec2(vec2)
            , m_vec4(vec4)
            , m_umax(umax)
            , m_F( f.integral(umax) ) // for normalization of the PSF
            , m_sum(0), m_count(0)
            , m_sumu(0) // average u, useful for calibration
            , m_pixels(0)
            , m_subset(subset)
            , m_back(background )
        {
            m_first = m_vec4.size()==0;
            if(debug){
                //psf_data = new std::ofstream("d:/users/burnett/temp/psf.txt");
                (*psf_data) << "u        f(u)      count    q" << std:: endl;
            }
        }
        ~Convert(){
             //if(debug){ psf_data->close();}
        }

        void operator()(const std::pair<HealPixel, int>& x){
            double diff =x.first().difference(m_dir); 
            double  u = sqr(diff/m_sigma)/2.;
            std::vector<int>::iterator it = find(m_vec4.begin(),m_vec4.end(),x.first.index());
            //return if 
            //1)normal mode and outside of cone 
            //2)first time in selection mode and the pixel is outside of cone
            //3)subsequent time in selection mode and the pixel is not in the original data set
            if((!m_subset&&u>m_umax)||(m_subset&&m_first&&u>m_umax)||(!m_first&&it==m_vec4.end()&&m_subset)) return;
            // just to see what is there
            // astro::SkyDir r(x.first()); double ra(r.ra()), dec(r.dec());
            double signal(m_f(u)/m_F)
                 , bkg(m_back(x.first()))
                 , q( bkg/(signal*m_umax-bkg)); 
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
#ifndef QBIN
            // todo: combine elements with vanishing t
            m_vec2.push_back(std::make_pair(q, x.second) );
#else
            if(q>0) {
                m_vec2.push_back(std::make_pair(q, x.second) );
            }
            else {
                int bin = (int)(-10./q);
                std::map<int,int>::iterator it = m_qbin.find(bin);
                if(it==m_qbin.end()) {
                    m_qbin.insert(std::map<int,int>::value_type(bin,x.second));
                } else {
                    ++(it->second);
                }
            }
#endif

            // save list of healpix id's
            if(m_first) m_vec4.push_back(x.first.index());
        }
        
        void consolidate() {
            for(std::map<int,int>::iterator it = m_qbin.begin();it!=m_qbin.end();++it) {
                double invq = 0.1*(it->first)+0.05;
                m_vec2.push_back(std::make_pair(-1/invq,it->second));
            }
            m_qbin.clear();
        }

        double average_f()const {return m_count>0? m_sum/m_count : -1.;}
        double average_u()const {return m_count>0? m_sumu/m_count : -1;}
        double average_b()const {return m_count>0? m_sumb/m_pixels: -1;}
        double count()const{return m_count;}

    private:
        SkyDir m_dir;
        PsfFunction& m_f;
        const astro::SkyFunction& m_back;
        double m_sigma;
        std::vector<std::pair<double, int> >& m_vec2;
        std::vector<int>& m_vec4;
        std::map<int,int> m_qbin;
        double m_umax, m_F;
        double m_sum, m_count, m_sumu, m_sumb, m_pixels;
        double m_back_norm;
        bool m_first;   //first call?
        bool m_subset;
    };

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // The log likelihood function
    class LogLike  {
    public:
        LogLike(double a):m_a(a){}

        double operator()(double prev, const std::pair<float, int>& x)
        {
            // std::cout <<  x.first  << " "<< x.second << std::endl;
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
class SimpleLikelihood::NormalizedBackground : public astro::SkyFunction {
public:
    NormalizedBackground(const SkyDir& dir, double angle, double emin, double emax){
        if( SimpleLikelihood::diffuse()!=0){
            SimpleLikelihood::diffuse()->setEnergyRange(emin, emax);
            m_back_norm = SimpleLikelihood::diffuse()->average(dir, angle, SimpleLikelihood::tolerance());
            if( m_back_norm==0){
                std::cerr << "Warning: normalization zero" << std::endl;
                m_back_norm=0.1; // kluge, like below
            }
        }
    }
    //! the normalized (in 0<u<umax)  background in the direction dir
    double operator()(const SkyDir& dir)const{
        double val(1.);
        if( SimpleLikelihood::diffuse()!=0){
            val=(*SimpleLikelihood::diffuse())(dir)/m_back_norm;
        }
        // a return value in the given direction
        if( val==0){ val=0.1;} // prevent zero, which is bad
        return val;
    }

private:
    double m_back_norm;
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SimpleLikelihood::SimpleLikelihood(const std::vector<std::pair<healpix::HealPixel, int> >& vec,
        const astro::SkyDir& dir, 
        double gamma, double sigma, double background, double umax, double emin,double emax)
        : m_vec(vec)
        , m_averageF(0)
        , m_psf(gamma)
        , m_sigma(sigma)
        , m_alpha(-1)
        , m_curv(-1)
        , m_background(background)  // default: no estimate
        , m_umax(umax)
        , m_emin(emin), m_emax(emax)
        , m_back(0)
{ 

    m_fint = m_psf.integral(m_umax);// integral out to umax

    // the integral of square, used by the quick estimator
    m_fint2 = m_psf.integralSquare(m_umax);

    setDir(dir);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SimpleLikelihood::~SimpleLikelihood()
{
    delete m_back;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SimpleLikelihood::setDir(const astro::SkyDir& dir, bool subset)
{
    m_dir = dir;
    // create set of ( 1/(f(u)/Fbar-1), weight) pairs in  m_vec2
    m_vec2.clear();
    delete m_back;
    m_back = new NormalizedBackground(dir, sqrt(2.*m_umax)*m_sigma, m_emin, m_emax);
    Convert conv(m_dir, m_psf, *m_back, m_sigma, m_umax, m_vec2, m_vec4,  subset);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);
    result.consolidate();

    m_photon_count = static_cast<int>(result.count());

    if( m_photon_count==0) return; //throw std::invalid_argument("SimpleLikelihood: no data after transform");

    m_averageF = result.average_f();
    m_avu = result.average_u();
    m_avb = result.average_b();

    // initialize to estimate.
    if( m_alpha<0) m_alpha = estimate();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double SimpleLikelihood::operator()( double a) const
{
    if( a<0) a=m_alpha;
    return std::accumulate(m_vec2.begin(), m_vec2.end(), poissonLikelihood(a), LogLike(a));
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::estimate() const
{
    if( m_photon_count==0) return -1;
    double est((m_umax*m_averageF-m_fint)/(m_umax*m_fint2-m_fint));
    if (est<0 ) est = 1e-4;
    if (est>1 ) est = 1-1e-4;
    return est;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::pair<double,double> SimpleLikelihood::maximize()
{
    static double tol(1e-3), tolchange(0.05);
    static int itermax(5);
    if( m_photon_count==0) return std::make_pair(0.,0.);
    double x(m_alpha); // starting point
    int iteration(0);
    std::pair<double, double> dw(100,100);

    for(double delta(0); fabs(dw.first) > tol && iteration<itermax; ++iteration, x-= delta) {

        // accumulate first and second derivatives
        dw = accumulate(m_vec2.begin(), m_vec2.end(), poissonDerivatives(x), Derivatives(x));
        if( debug) std::cout 
            << std::setw(5)  << iteration 
            << std::setw(10) << std::setprecision(3) << x 
            << std::setw(15) << dw.first 
            << std::setw(15) << dw.second << std::endl;


        if( fabs(dw.second) < 1e-6 ) {
            dw.second = 1.0;
            break;
        }
        delta =dw.first/dw.second;

        // prevent from going negative or greater than 1
        if( delta > x ) delta = x/2.; 
        if( x-delta>1 ) delta = (x-1.)/2.;

        // quit if change small compared with resolution
        double sigdel(1./sqrt(dw.second));
        if( fabs(delta/sigdel) < tolchange ){x-=delta; break;}
    }
    m_alpha = x;
    m_sigma_alpha = 1./sqrt(dw.second);
    return std::make_pair(x, m_sigma_alpha);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::poissonLikelihood(double a)const
{
    if( m_background<0) return 0; 
    double expect(signal(a)+background());
    if( expect<=0) return 0; //?
    return expect - m_photon_count*log(expect);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::pair<double,double> SimpleLikelihood::poissonDerivatives(double a)
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
Hep3Vector SimpleLikelihood::gradient() const
{
    Hep3Vector grad;
    m_curv=m_w=0;

    double sig2( sqr(m_sigma) );
    double gamma( m_psf.gamma() );
    //double w(0); // -log likelihood (check)
    int count(0);

    Hep3Vector perp(m_dir().orthogonal());

    std::vector<std::pair<healpix::HealPixel,int> >::const_iterator it = m_vec.begin();

    for( ; it< m_vec.end(); ++it){
        const std::pair<healpix::HealPixel,int>& h = *it;

        SkyDir d( h.first );
        int nphoton( h.second);
        Hep3Vector delta( m_dir() - d() ); 
        double u( 0.5*delta.mag2()/sig2);
        if(u>m_umax) continue;
        if((u>m_umax&&m_vec4.size()==0)||(find(m_vec4.begin(),m_vec4.end(),h.first.index())==m_vec4.end()&&m_vec4.size()!=0)) continue;
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
double SimpleLikelihood::curvature() const
{
    if( m_curv==-1.) gradient(); // in not initialized. (negative curvature is bad)
    return m_curv;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::TS(double alpha)const { 
    return 2.*((*this)(0) - (*this)(alpha));
    //return -2*(*this)(alpha);
}
double SimpleLikelihood::solidAngle()const{
    return 2*M_PI* sqr(m_sigma)*m_umax;
}

double SimpleLikelihood::feval(double k) {

    if(!m_vec.size()) return -1.0;
    double F = m_psf.integral(k*m_umax);
    double acc = 0;
    for(std::vector<std::pair<HealPixel,int> >::const_iterator ite=m_vec.begin();ite!=m_vec.end();++ite) {
        double diff =ite->first().difference(m_dir); 
        double u = sqr(diff/m_sigma)/2.;
        if(u>m_umax) continue;
        // just to see what is there
        // astro::SkyDir r(x.first()); double ra(r.ra()), dec(r.dec());
        double f = m_psf(k*u);
        acc-=ite->second*log(m_alpha*f*k/F+(1-m_alpha)/(m_umax));
    }
    return acc;
}

double SimpleLikelihood::geval(double k) {
   if(!m_vec.size()) return -1.0;
    m_vec2.clear();
    PsfFunction ps(k*m_psf.gamma());
    Convert conv(m_dir, ps, *m_back, m_sigma, m_umax, m_vec2, m_vec4, true);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);
    result.consolidate();
    double ts = -TS(m_alpha);
    m_vec2.clear();
    return ts;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::operator()(const astro::SkyDir& dir)const
{
    double diff( dir.difference(m_dir) ); 
    double u   ( sqr(diff/m_sigma)/2.  );

    if( s_display_mode>0 ){
        // residual
        if( u> m_umax ) return 0;
        // get the pixel
        int level = m_vec[0].first.level();
        healpix::HealPixel pixel(m_dir, level);
        std::vector<std::pair<healpix::HealPixel,int> >::const_iterator it;
#if 0 // fails to compile?
        it=  std::find(m_vec.begin(), m_vec.end(), pixel);
#else
        for( it=m_vec.begin(); it!=m_vec.end(); ++it ) {
            if( it->first == pixel) break;
        }
#endif
        if( it==m_vec.end() ) return 0;

        int counts (it->second);
        if( s_display_mode==1 ) return counts;   // observed in the pixel
        double back( (*m_back)(dir) );           // normalized background prediction at the pixel
        if( s_display_mode==2 ) return back; 
        

    }
    // density version
    // just to see what is there
    // astro::SkyDir r(x.first()); double ra(r.ra()), dec(r.dec());

    return signal()*m_psf(u)/ m_fint/(2*M_PI*sqr(m_sigma));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


skymaps::SkySpectrum* SimpleLikelihood::diffuse()
{
    return s_diffuse;
}

void SimpleLikelihood::setDiffuse(skymaps::SkySpectrum* diff)
{
    s_diffuse = diff;
}

double SimpleLikelihood::tolerance()
{
    return s_tolerance;
}

double SimpleLikelihood::defaultUmax()
{
    return s_defaultUmax;
}
void SimpleLikelihood::setDefaultUmax(double umax)
{
    s_defaultUmax = umax;
}

void SimpleLikelihood::setTolerance(double tol)
{
    s_tolerance = tol;
}

void SimpleLikelihood::recalc(bool subset) 
{
    // create set of ( 1/(f(u)/Fbar-1), weight) pairs in  m_vec2
    m_vec2.clear();
    Convert conv(m_dir, m_psf, *m_back, m_sigma, m_umax, m_vec2,  m_vec4, subset);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);
    result.consolidate();

    if(m_photon_count!= static_cast<int>(result.count())) {
        m_photon_count=m_photon_count;
    } 

    if( m_photon_count==0) return; //throw std::invalid_argument("SimpleLikelihood: no data after transform");

    m_averageF = result.average_f();
    m_avu = result.average_u();
    m_avb = result.average_b();
}