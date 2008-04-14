/** @file SimpleLikelihood.cxx
@brief Implementation of class SimpleLikelihood

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SimpleLikelihood.cxx,v 1.25 2008/03/31 09:01:37 burnett Exp $
*/

#include "pointlike/SimpleLikelihood.h"
#include "skymaps/DiffuseFunction.h"
#include "astro/SkyDir.h"
#include <algorithm>
#include <cmath>
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
double SimpleLikelihood::s_defaultUmax =50;
#define UMAXCUT

skymaps::SkySpectrum* SimpleLikelihood::s_diffuse(0);
double  SimpleLikelihood::s_tolerance(0.05); // default:

namespace {

    bool debug(false);
    double tolerance(0.05); // integral tolerance
    inline double sqr(double x){return x*x;}
    std::ostream * psf_data = &std::cout;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class Convert  {
    public:
        Convert( const SkyDir& dir, PsfFunction& f, 
            double sigma, double umax,
            std::vector< std::pair<double, int> >& vec2,
            std::vector<int>& vec4,
            double emin, double emax, bool subset
            )
            : m_dir(dir), m_f(f), m_sigma(sigma)
            , m_vec2(vec2)
            , m_vec4(vec4)
            , m_umax(umax)
            , m_F( f.integral(umax) ) // for normalization of the PSF
            , m_sum(0), m_count(0)
            , m_sumu(0) // average u, useful for calibration
            , m_pixels(0)
            , m_back_norm(1)
            , m_subset(subset)
        {
            if( SimpleLikelihood::diffuse()!=0){
                SimpleLikelihood::diffuse()->setEnergyRange(emin, emax);
                double angle(sqrt(2.*umax)*sigma);
                m_back_norm = SimpleLikelihood::diffuse()->average(dir, angle, SimpleLikelihood::tolerance());
                if( m_back_norm==0){
                    std::cerr << "Warning: normalization zero" << std::endl;
                    m_back_norm=0.1; // kluge, like below
                }
            }
            m_first = m_vec4.size()==0;
            if(debug){
                //psf_data = new std::ofstream("d:/users/burnett/temp/psf.txt");
                (*psf_data) << "u        f(u)      count    q" << std:: endl;
            }
        }
        ~Convert(){
            //if(debug){ psf_data->close();}
        }
        //! the normalized (in 0<u<umax)  background in the direction dir
        double b(const SkyDir& dir){
            double val(1.);
            if( SimpleLikelihood::diffuse()!=0){
                val=(*SimpleLikelihood::diffuse())(dir)/m_back_norm;
            }
            // a return value in the given direction
            if( val==0){ val=0.1;} // prevent zero, which is bad
            return val;
        }

        void operator()(const std::pair<HealPixel, int>& x){
            double diff =x.first().difference(m_dir); 
            double  u = sqr(diff/m_sigma)/2.;
#ifndef  UMAXCUT // old
            if( u>m_umax ) return;
#else // new
            std::vector<int>::iterator it = find(m_vec4.begin(),m_vec4.end(),x.first.index());
            //return if 
            //1)normal mode and outside of cone 
            //2)first time in selection mode and the pixel is outside of cone
            //3)subsequent time in selection mode and the pixel is not in the original data set
            if((!m_subset&&u>m_umax)||(m_subset&&m_first&&u>m_umax)||(!m_first&&it==m_vec4.end()&&m_subset)) return;
#endif
            // just to see what is there
            // astro::SkyDir r(x.first()); double ra(r.ra()), dec(r.dec());
            double signal(m_f(u)/m_F)
                , bkg(b(x.first()))
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

            // todo: combine elements with vanishing t
            m_vec2.push_back(std::make_pair(q, x.second) );
            if(m_first) m_vec4.push_back(x.first.index());
        }
        double average_f()const {return m_count>0? m_sum/m_count : -1.;}
        double average_u()const {return m_count>0? m_sumu/m_count : -1;}
        double average_b()const {return m_count>0? m_sumb/m_pixels: -1;}
        double count()const{return m_count;}

    private:
        SkyDir m_dir;
        PsfFunction& m_f;
        double m_sigma;
        std::vector<std::pair<double, int> >& m_vec2;
        std::vector<int>& m_vec4;
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
{ 

    m_fint = m_psf.integral(m_umax);// integral out to umax

    // the integral of square, used by the quick estimator
    m_fint2 = m_psf.integralSquare(m_umax);

    setDir(dir);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SimpleLikelihood::setDir(const astro::SkyDir& dir, bool subset)
{
    m_dir = dir;
    // create set of ( 1/(f(u)/Fbar-1), weight) pairs in  m_vec2
    m_vec2.clear();

    Convert conv(m_dir, m_psf, m_sigma, m_umax, m_vec2, m_vec4, m_emin, m_emax, subset);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);

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
#ifndef UMAXCUT // original version
        if(u>m_umax) continue;
#else // Marshal mod
        if((u>m_umax&&m_vec4.size()==0)||(find(m_vec4.begin(),m_vec4.end(),h.first.index())==m_vec4.end()&&m_vec4.size()!=0)) continue;
#endif
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
    Convert conv(m_dir, ps, m_sigma, m_umax, m_vec2, m_vec4, m_emin, m_emax,true);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);
    m_vec2.clear();
    return -TS(m_alpha);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::operator()(const astro::SkyDir& dir)const
{
    double diff =dir.difference(m_dir); 
    double  u = sqr(diff/m_sigma)/2.;
    // just to see what is there
    // astro::SkyDir r(x.first()); double ra(r.ra()), dec(r.dec());

    return signal()*m_psf(u)/ m_fint/(2*M_PI*sqr(m_sigma));
}


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
    Convert conv(m_dir, m_psf, m_sigma, m_umax, m_vec2, m_vec4, m_emin, m_emax, subset);
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);

    if(m_photon_count!= static_cast<int>(result.count())) {
        m_photon_count=m_photon_count;
    } 

    if( m_photon_count==0) return; //throw std::invalid_argument("SimpleLikelihood: no data after transform");

    m_averageF = result.average_f();
    m_avu = result.average_u();
    m_avb = result.average_b();
}