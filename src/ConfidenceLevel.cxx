/** @file ConfidenceLevel.cxx
@brief implementation of ConfidenceLevel

$Header$

*/

#include "pointlike/ConfidenceLevel.h"
#include "pointlike/SimpleLikelihood.h"
#include "skymaps/PsfFunction.h"

#include "skymaps/Band.h"
#include "astro/SkyDir.h"
#include <cmath>

using namespace pointlike;
using astro::SkyDir;
using skymaps::Band;
using skymaps::PsfFunction;


namespace {

    int ran_seed(-1); // seed to use wiht ran1
    // following copied from numerical recipes
    double gammln(const double xx)
    {
        int j;
        double x,y,tmp,ser;
        static const double cof[6]={76.18009172947146,-86.50532032941677,
            24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
            -0.5395239384953e-5};

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<6;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
    }

    double ran1(int &idum)
    {
        const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
        const int NDIV=(1+(IM-1)/NTAB);
        const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
        static int iy=0;
        //static Vec_INT iv(NTAB);
        //static std::vector<int>iv(NTAB);
        static int iv[NTAB];
        int j,k;
        double temp;

        if (idum <= 0 || !iy) {
            if (-idum < 1) idum=1;
            else idum = -idum;
            for (j=NTAB+7;j>=0;j--) {
                k=idum/IQ;
                idum=IA*(idum-k*IQ)-IR*k;
                if (idum < 0) idum += IM;
                if (j < NTAB) iv[j] = idum;
            }
            iy=iv[0];
        }
        k=idum/IQ;
        idum=IA*(idum-k*IQ)-IR*k;
        if (idum < 0) idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = idum;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
    }

    double poidev(const double xm, int &idum)
    {
        const double PI=3.141592653589793238;
        static double sq,alxm,g,oldm=(-1.0);
        double em,t,y;

        if (xm < 12.0) {
            if (xm != oldm) {
                oldm=xm;
                g=exp(-xm);
            }
            em = -1;
            t=1.0;
            do {
                ++em;
                t *= ran1(idum);
            } while (t > g);
        } else {
            if (xm != oldm) {
                oldm=xm;
                sq=sqrt(2.0*xm);
                alxm=log(xm);
                g=xm*alxm-gammln(xm+1.0);
            }
            do {
                do {
                    y=tan(PI*ran1(idum));
                    em=sq*y+xm;
                } while (em < 0.0);
                em=floor(em);
                t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
            } while (ran1(idum) > t);
        }
        return em;
    }

}

ConfidenceLevel::ConfidenceLevel(const SimpleLikelihood& slike, 
                                 const SkyDir& dir,
                                 int realizations)
: m_slike(slike)
, m_realizations(realizations)
{
    setup(dir);
}

void ConfidenceLevel::setup(const SkyDir& dir)
{
    int iband(0); // for now

    const Band& band(m_slike.band());
    m_events = m_slike.photons(); 
    double radius(band.sigma()*sqrt(2.*m_slike.umax()));

    // get a list of pixels and counts, including the empty ones
    band.query_disk( dir, radius, m_pixels, true); 
    unsigned int npix(m_pixels.size());

    // determine expected normalized background for each
    m_back.reserve(npix);
    m_signal.reserve(npix);

    // psf function used (get parameters from band)
    double sigma(band.sigma()), gamma(band.gamma());
    PsfFunction psf(gamma);

    // get the predictions for each pixel
    double backsum(0), sigsum(0);
    for( std::vector<std::pair<SkyDir, int> >::const_iterator ib(m_pixels.begin()); ib!=m_pixels.end(); ++ib){
        const std::pair<SkyDir,int>& pixel(*ib);
        double backval(m_slike.background_function()(pixel.first)),
            sigval(psf( dir, pixel.first, sigma));
        m_back.push_back(backval);
        m_signal.push_back(sigval);
        backsum += backval;
        sigsum += sigval;
    }

    // finally normalize
    for( unsigned int i(0); i<npix; ++i)  {
        m_signal[i]/=sigsum;
        m_back[i] /=backsum;
    }
}

double ConfidenceLevel::likelihood(double alpha, bool save) const
{
    double w(0);
    if( alpha<0 ) alpha=0;
    else if(alpha>1.0) alpha=1;
    for( unsigned int i(0); i<m_signal.size(); ++i ){
        int n( m_pixels[i].second );
        if( n>0 || save ) {
            double value(alpha * m_signal[i] + (1-alpha) * m_back[i] );
            if(save) m_predict.push_back(value);
            if(n>0)  w -= n * log( value );
        }
    }
    return w;
}

// predicate for search
class not_greater{
public:
    not_greater(double x):m_x(x){};
    double operator()(double y)const{ 
#if 1
        bool ret(y>m_x);
        if(ret){
            return true;
        }else{
            return false;
        }
#else
        return y>m_x; 
#endif
    }
    double m_x;
};
double ConfidenceLevel::operator()() const
{

   double alpha(m_slike.alpha()); // the fit fraction from the fit
   double sigma_alpha(m_slike.sigma_alpha());

   // check that fit is ok: (error should increase -loglike by 0.5)
   double w( likelihood(alpha, true) ); // recalculate and save predicted values
#if 0 // keep for testing for now
   double wa( likelihood(alpha-sigma_alpha)-w ),
    wb( likelihood(alpha+sigma_alpha)-w ),
    w2a( likelihood(alpha-2*sigma_alpha)-w ),
    w2b( likelihood(alpha+2*sigma_alpha)-w );

   // compare with the original

    double y0 ( m_slike() ), 
        ya( m_slike(alpha-sigma_alpha)-y0 ), 
        yb( m_slike(alpha+sigma_alpha)-y0 ),
        y2a( m_slike(alpha-2*sigma_alpha)-y0 ), 
        y2b( m_slike(alpha+2*sigma_alpha)-y0 );
#endif

    // generate experiments
    std::vector<double> wlist(m_realizations);
    for(std::vector<double>::iterator i(wlist.begin()); i< wlist.end(); ++i){
       *i= realize(); 
    }
    sort(wlist.begin(), wlist.end());
#if 1
    int k(0);
    for(std::vector<double>::iterator i(wlist.begin()); i< wlist.end(); ++i, ++k){
        if( w< *i) break;
    }
    double c(double(k)/m_realizations);
    
#else // why does this not invoke the local not_greater?
    std::vector<double>::iterator wloc( 
        std::find_if(wlist.begin(), wlist.end(), 
            not_greater(w)
        ));
    double c( 1.0- (wloc-wlist.begin())/m_realizations);

#endif
    return c;

}

double ConfidenceLevel::realize()const
{
    double w(0);
    for( unsigned int i(0); i<m_signal.size(); ++i ){
        double mu(m_predict[i]*m_events );
        int n( poidev(mu, ran_seed ));
        if( n>0 ) {
            w -= n * log( m_predict[i] );
        }
    }
    return w;
}

