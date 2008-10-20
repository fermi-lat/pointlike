/** @file SimpleLikelihood.cxx
@brief Implementation of class SimpleLikelihood

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SimpleLikelihood.cxx,v 1.48 2008/10/20 02:58:32 burnett Exp $
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
using skymaps::SkySpectrum;
using skymaps::DiffuseFunction;
using skymaps::PsfFunction;

using namespace pointlike;

//#define DEBUG_PRINT
double SimpleLikelihood::s_defaultUmax =50;

double  SimpleLikelihood::s_tolerance(0.05); // default:

bool SimpleLikelihood::s_enable_extended(false);
void SimpleLikelihood::enable_extended_likelihood(bool q)
{
    s_enable_extended=q;
}
bool SimpleLikelihood::extended_likelihood(){return s_enable_extended;}



namespace {

    bool debug(false);
    double tolerance(0.05); // integral tolerance
    inline double sqr(double x){return x*x;}
    std::ostream * psf_data = &std::cout;

    // for the binning in 1/q when q is negative, pixels far from source
    double binsize(0.005); // 0.025); // bin size for 1/q. set zero to disable
    std::map<double,int> qmap; // used for binning 1/q


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    class Convert  {
    public:
        Convert( const SkyDir& dir, const PsfFunction& f, 
            const astro::SkyFunction& background,
            double sigma, double umax,
            std::vector< std::pair<double, int> >& vec2 
            )
            : m_dir(dir), m_f(f), m_sigma(sigma)
            , m_vec2(vec2)
            //, m_vec4(vec4)
            , m_umax(umax)
            , m_F( f.integral(umax) ) // for normalization of the PSF
            , m_sum(0), m_count(0)
            , m_sumu(0) // average u, useful for calibration
            , m_pixels(0)
            , m_back(background )
            , m_maxu_found(0)
        {
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
        const PsfFunction& m_f;
        const astro::SkyFunction& m_back;
        double m_sigma;
        std::vector<std::pair<double, int> >& m_vec2;
        //std::vector<int>& m_vec4;
        double m_umax, m_F;
        double m_sum, m_count, m_sumu, m_sumb, m_pixels;
        double m_back_norm;
        bool m_first;   //first call?
        //bool m_subset;
        double m_maxu_found;
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
// private nested class to  normalized background
class SimpleLikelihood::NormalizedBackground : public astro::SkyFunction {
public:
    NormalizedBackground(const astro::SkyFunction* diffuse, const SkyDir& dir, double angle)
    : m_diffuse(diffuse)
    {
        if( m_diffuse!=0){
            m_back_norm = m_diffuse->average(dir, angle, SimpleLikelihood::tolerance());
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
            val=(*m_diffuse)(dir)/m_back_norm;
        }
        // a return value in the given direction
        if( val==0){ val=0.1;} // prevent zero, which is bad
        return val;
    }

    //! The expected background, or total normalization (integral of the background function)
    double total()const{return m_back_norm;}

private:
    const astro::SkyFunction* m_diffuse;
    double m_back_norm;
};


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SimpleLikelihood::SimpleLikelihood(const skymaps::Band& band,
        const astro::SkyDir& dir, 
        double umax, 
        const astro::SkyFunction* background)
        : m_averageF(0)
        , m_psf(band.gamma())
        , m_alpha(-1)
        , m_curv(-1)
        , m_umax(umax)
        , m_sigma(band.sigma())
        , m_back(0)
        , m_diffuse(background)
{ 
    m_bands.push_back(&band);

    m_fint =  m_psf.integral(umax) ; // for normalization of the PSF

    // the integral of square, used by the quick estimator
    m_fint2 = m_psf.integralSquare(m_umax);

    m_dir = dir;
    setDir(dir);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SimpleLikelihood::addBand(const skymaps::Band& moredata)
{
    m_bands.push_back(&moredata);
    //reload();
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
    reload(subset);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void SimpleLikelihood::reload(bool subset)
{
    using skymaps::Band;

    if( ! subset) {
        // filll m_vec with weighted pixels unless operating on the current set
        m_vec.clear();
        delete m_back;
        m_back = new NormalizedBackground(m_diffuse, m_dir, sqrt(2.*m_umax)*sigma());
        // select pixels within u_max, sum bands
        for( std::vector<const Band*>::const_iterator bandit(m_bands.begin()); bandit!=m_bands.end();++bandit){ 
            (*bandit)->query_disk(m_dir, sigma()*sqrt(2.*m_umax), m_vec);
        }
    }    

    // create set of ( 1/(f(u)/Fbar-1), weight) pairs in  m_vec2
    m_vec2.clear();
    Convert conv(m_dir, m_psf, *m_back, sigma(), m_umax, m_vec2); 
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
    double v(0); 
    if( extended_likelihood() )
        v = poissonLikelihood(a); 

    return std::accumulate(m_vec2.begin(), m_vec2.end(),  v, LogLike(a));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double SimpleLikelihood::logLikelihood( double count) const
{
    double denom( photons() );
    if( background()>0 ) denom=count+background();
    return operator()(count/denom);
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
std::pair<double, double> SimpleLikelihood::derivatives(double x)
{
    std::pair<double,double> v( std::make_pair(0.,0.) );
    if( extended_likelihood() ){
        v = poissonDerivatives(x);
    }
    return accumulate(m_vec2.begin(), m_vec2.end(), v,   Derivatives(x));
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::pair<double,double> SimpleLikelihood::maximize()
{
    static double tol(1e-3), tolchange(0.05);
    static int itermax(5);
    if( m_vec.size() ==0 ) reload(); // first time maybe
    if( m_photon_count==0) return std::make_pair(0.,0.);
    double x(m_alpha); // starting point
    int iteration(0);
    std::pair<double, double> dw(100,100);

    for(double delta(0); fabs(dw.first) > tol && iteration<itermax; ++iteration, x-= delta) {

        // accumulate first and second derivatives
        dw = derivatives(x);
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

double SimpleLikelihood::background()const{return m_back!=0? m_back->total()*solidAngle(): 0;}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::poissonLikelihood(double a)const
{
    double expect(signal(a)+background());
    if( expect<=0) return 0; //?
    return expect - m_photon_count*log(expect);

}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
std::pair<double,double> SimpleLikelihood::poissonDerivatives(double a)
{
    double d1(0), d2(0);

    double t(a + background()/m_photon_count);
    d1 = m_photon_count*(1-1/t);  // first log likelihood derivative
    d2 = m_photon_count/t/t;      // second
    return std::make_pair(d1, d2);
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hep3Vector SimpleLikelihood::gradient() const
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
    return 2*M_PI* sqr(sigma())*m_umax;
}

double SimpleLikelihood::feval(double k) {
    if(!m_vec.size()) return -1.0;
    double F = m_psf.integral(k*m_umax);
    double acc = 0;
    for(PixelList::const_iterator ite=m_vec.begin();ite!=m_vec.end();++ite) {
        double diff =ite->first.difference(m_dir); 
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
    //maximize();
    double F = ps.integral(m_umax);
    double acc = 0;
    for(PixelList::const_iterator ite=m_vec.begin();ite!=m_vec.end();++ite) {
        double diff =ite->first.difference(m_dir); 
        double u = sqr(diff/m_sigma)/2.;
        if(u>m_umax) continue;
        // just to see what is there
        // astro::SkyDir r(x.first()); double ra(r.ra()), dec(r.dec());
        double f = ps(u);
        acc-=ite->second*log(m_alpha*f/F+(1-m_alpha)/(m_umax));
    }
    
    m_vec2.clear();
    return acc;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::operator()(const astro::SkyDir& dir)const
{
    double diff( dir.difference(m_dir) ); 
    double u   ( sqr(diff/sigma())/2.  );
    double jacobian( 2.*M_PI* sqr(sigma()) );
    return signal()*m_psf(u)/ m_fint/jacobian;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
double SimpleLikelihood::display(const astro::SkyDir& dir, int mode) const
{
    SimpleLikelihood* self = const_cast<SimpleLikelihood*>(this);
    if( mode==5 ) return self->TSmap(dir);

    double diff( dir.difference(m_dir) ); 
    double u   ( sqr(diff/sigma())/2.  );
    double jacobian( 2.*M_PI* sqr(sigma()) );
    if( mode ==0){
        // default is density, as above
        return signal()*m_psf(u)/ m_fint/jacobian;
    }
    SkyDir pdir(band().dir(band().index(dir)));
    u = sqr( pdir.difference(m_dir)/sigma())/2.; // recalculate u for pixel
    if( u> m_umax ) return 0;

    // now look for it in the data
    PixelList::const_iterator it;
    int index( band().index(pdir));
    for( it=m_vec.begin(); it!=m_vec.end(); ++it ) {
        if( band().index(it->first) == index) break;
    }

    int counts (it==m_vec.end()? 0:  it->second);

    if( mode==1 ) return counts;   // observed in the pixel
    double back( (*m_back)(pdir));           // normalized background at the pixel: average should be 1.0
    if( mode==2 ) return back; 
    double area( band().pixelArea() );

    double prediction ( area/jacobian * photons() * ( m_alpha*m_psf(u)/m_fint + (1.-m_alpha)*back/m_umax ) );
    if( mode==3 ) return prediction; // prediction of fit
    if( mode==4 ) return counts-prediction;
    return 0; 
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


const astro::SkyFunction* SimpleLikelihood::diffuse() const
{
    return m_diffuse;
}

void SimpleLikelihood::setDiffuse(astro::SkyFunction* diff)
{
    m_diffuse = diff;
    setDir(m_dir); // reset data
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
    Convert conv(m_dir, m_psf, *m_back, sigma(), m_umax, m_vec2); //,  m_vec4, subset);
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

double SimpleLikelihood::TSmap(astro::SkyDir sdir)const
{
    SimpleLikelihood* self = const_cast<SimpleLikelihood*>(this);
#if 0
    SkyDir old_dir(m_dir);
    m_dir = sdir;
    self->reload( true); // recalculate
    double ts( TS() ); 
    m_dir=old_dir;
    self->reload( true); // should restore
    return ts;

#else
    std::vector<std::pair<double, int> > vec2;  //stores <log-like,nphotons>
    double oldbinsize(binsize);
    binsize =0; // disable the binning, makes jumps
    // load vec2 from current list of pixels and this direction
    Convert conv(sdir, m_psf, *m_back, sigma(), m_umax, vec2 );
    Convert result=std::for_each(m_vec.begin(), m_vec.end(), conv);
    result.consolidate();

    // now compute TS, using current alpha, but with vec2
    double ret =2.*(
        std::accumulate(vec2.begin(), vec2.end(),  poissonLikelihood(0), LogLike(0))
        -std::accumulate(vec2.begin(), vec2.end(),  poissonLikelihood(m_alpha), LogLike(m_alpha))
        );
    binsize = oldbinsize; // restore binning 
    return ret;

#endif
}

