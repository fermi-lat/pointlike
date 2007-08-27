/** @file UbPointSourceLikelihood.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/UbPointSourceLikelihood.cxx,v 1.5 2007/07/14 22:09:20 burnett Exp $

*/

#include "pointlike/UbPointSourceLikelihood.h"

#include "map_tools/PhotonMap.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>


using namespace astro;
using namespace pointlike;
namespace {

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // from fits 2/26/2006 
    double fit_gamma[]={0,0,0,0,0,0,
        2.27,2.22,2.25,2.25,2.29,2.14,2.02,1.87};
#if 0 // old
    double fit_sigma[]={0,0,0,0,0,0,
        0.335,0.319,0.332,0.352,0.397,0.446,0.526,0.657};
#else // from empirical study of a single source 
    // set level 8 -13 using handoff (7/14/07) tweak level 12 from .52 to .62
    double fit_sigma[]={0,0,0,0,0,0,
        0.335,0.319, //0.422, 0.9, 0.9, 1.0, 1.1, 1.0 
    0.42899897,  0.45402442,  0.4742285,   0.60760651,  0.62,  0.94671
    };

#endif
    // the scale_factor used: 2.5 degree at level 6, approximately the 68% containment
    int base_level(6);
    double scale_factor(int level){return 2.5*pow(2.0, base_level-level)*M_PI/180.;}

    double s_TScut(2.);  // only combine energy bands
}
// set these from preliminary data above
std::vector<double> UbPointSourceLikelihood::gamma_level(fit_gamma, fit_gamma+sizeof(fit_gamma)/sizeof(double)); 
std::vector<double> UbPointSourceLikelihood::sigma_level(fit_sigma, fit_sigma+sizeof(fit_gamma)/sizeof(double)); 

UbPointSourceLikelihood::UbPointSourceLikelihood(const std::vector<astro::Photon> m_data_vec,
                                             std::string name,
                                             const astro::SkyDir& dir, 
                                             double radius,
                                             int minlevel, int maxlevel
                                             )
                                             : m_name(name)
                                             , m_dir(dir)
                                             , m_verbose(false)
                                             , m_out(&std::cout)
{
    int begin(minlevel), end(maxlevel+1); // levels: extract from data?

    for( int level=begin; level<end; ++level){
        
        std::vector<astro::Photon> cut = extract(m_data_vec,level);
        // get PSF parameters from fits
        double gamma( gamma_level[level] );
        // and create the simple likelihood object
        (*this)[level] = new UbSimpleLikelihood(cut, dir, gamma);
        if( false ) { // make table of parameters
            out() << std::setw(6) << level 
                << " " << std::setw(10) << std::left << gamma 
                << std::setw(10) 
                << std::right << std::endl;
        }

    }
}
UbPointSourceLikelihood::~UbPointSourceLikelihood()
{
    for( iterator it = begin(); it!=end(); ++it){
        delete it->second;
    }
}

std::vector<astro::Photon> UbPointSourceLikelihood::extract(const std::vector<astro::Photon> data,int level) {
    std::vector<astro::Photon> result;
    double lowsig = 0.058*pow(2.,1.*(6-level));
    double hisig = lowsig/2;
    for(std::vector<astro::Photon>::const_iterator it = data.begin();it!=data.end();++it) {
        double p0,p1;
        if(it->eventClass()==0) {
            p0=0.058;
            p1=0.000377;
        } else {
            p0=0.096;
            p1=0.0013;
        }
        double sigma = sqrt((p0*pow(it->energy()/100,-0.8))*(p0*pow(it->energy()/100,-0.8))+p1*p1);
        if(sigma<=lowsig&&sigma>hisig) {
            result.push_back((*it));
        }
    }
    return result;
}

double UbPointSourceLikelihood::maximize(int skip)
{
    m_TS = 0;
    iterator it = begin();
    for( int i = 0; i< skip; ++it, ++i);

    for( ; it!=end(); ++it){
        it->second->maximize();
        m_TS+= it->second->TS();
    }
    return m_TS;
}

void UbPointSourceLikelihood::setBackgroundDensity(const std::vector<double>& density)
{
    std::vector<double>::const_iterator id = density.begin();
    for( iterator it = begin(); it!=end(); ++it, ++id){
        double bk(*id);
        it->second->setBackgroundDensity(bk);
    }
}

void UbPointSourceLikelihood::setDir(const astro::SkyDir& dir){
    for( iterator it = begin(); it!=end(); ++it){
        it->second->setDir(dir);
    }
    m_dir = dir;
}

Hep3Vector UbPointSourceLikelihood::gradient(int skip) const{
    Hep3Vector t(0);
    const_iterator it = begin();
    for( int i = 0; i< skip; ++it, ++i);
    for( ; it!=end(); ++it){
        if( it->second->TS()< s_TScut) continue;
        Hep3Vector grad(it->second->gradient());
        double curv(it->second->curvature());
        if( curv > 0 ) t+= grad;
    }
    return t;
}

double UbPointSourceLikelihood::curvature(int skip) const{
    double t(0);
    const_iterator it = begin();
    for( int i = 0; i< skip; ++it, ++i);
    for( ; it!=end(); ++it){
        if( it->second->TS()< s_TScut) continue;
        double curv(it->second->curvature());
        if( curv>0 )  t+= curv;
    }
    return t;
}

void UbPointSourceLikelihood::printSpectrum()
{

    using std::setw; using std::left; using std::setprecision; 
    if( verbose() ){
        out() << "\nSpectrum of source " << m_name << " at ra, dec=" 
            << setprecision(6) << m_dir.ra() << ", "<< m_dir.dec() << std::endl;

        out() << "level events  backgnd  sig fraction  TS " << std::endl;
        //  level events  sig fraction    TS
        //    6  592  0.56 +/-  0.036     193.9
    }
    m_TS =0;
    for( const_iterator it = begin(); it!=end(); ++it){

        UbSimpleLikelihood& levellike = *it->second;
        int level = it->first;

        if( verbose() ){
            double bkg(levellike.background());
            out()  << std::setw(5) << std::fixed << level 
                   << setw(6) << levellike.photons()
                    << setw(10);
            if(bkg>=0) {
                out() << setprecision(1) << levellike.background();
            }else{
                out() << "     -    ";
            }
            if( levellike.photons()==0)  out() << std::endl; 
        }

        if( levellike.photons()==0) {
            continue;
        }
        std::pair<double,double> a(levellike.maximize());
        double ts(levellike.TS()); m_TS+=ts;

        if( verbose() ){
            out() << setprecision(2) << setw(6)<< a.first<<" +/- "
                << std::setw(4)<< a.second 
                << setw(6)<< setprecision(0)<< ts << std::endl;
        }
    }
    if( verbose() ){
        out() << std::setw(32) << m_TS << std::endl;
    }
}
double UbPointSourceLikelihood::localize(int skip1, int skip2)
{
    double t(100);
    for( int skip=skip1; skip<=skip2; ++skip){
        t = localize(skip);
        if (t<1) return t;
    }
    return t;
}

double UbPointSourceLikelihood::localize(int skip)
{
    using std::setw; using std::left; using std::setprecision; using std::right;
    using std::fixed;
    int wd(10), iter(0), maxiter(10);
 
        if( verbose()){
            out() 
                << "      Searching for best position, skipping first "<< skip <<" layers \n"
                << setw(wd) << left<< "Gradient   " 
                << setw(wd) << left<< "delta  "   
                << setw(wd) << left<< "ra"
                << setw(wd) << left<< "dec "
                << setw(wd) << left<< "error "
                << setw(wd) << left<< "Ts "
                <<std::endl;
        }
        double oldTs( maximize(skip)); // initial (partial) TS

        for( ; iter<maxiter; ++iter){
            Hep3Vector grad( gradient(skip) );
            double     curv( curvature(skip) );

            // check that resolution is ok: if curvature gets small or negative we are lost
            if( curv < 1.e4){
                if( verbose()) out() << "  >>>aborting, lost" << std::endl;
                return 98.;
                break;
            }
            double     sig( 1/sqrt(curv))
                ,      gradmag( grad.mag() )
                ;
            //,      oldTs( TS() );
            Hep3Vector delta = grad/curv;

            if( verbose() ){
                out() << fixed << setprecision(0)
                    <<  setw(wd-2)  << right<< gradmag << "  "
                    <<  setprecision(4)
                    <<  setw(wd) << left<< delta.mag()*180/M_PI
                    <<  setw(wd) << left<< m_dir.ra()
                    <<  setw(wd) << left<< m_dir.dec() 
                    <<  setw(wd) << left<< sig*180/M_PI
                    <<  setw(wd) << left <<setprecision(0)<< oldTs
                    <<  right <<setprecision(3) <<std::endl;
            }
            // check for too large step, limit to 3 sigma
            if( delta.mag() > 3.* sig) delta = 3.*sig* delta.unit(); 

            // here decide to back off if likelihood does not increase
            Hep3Vector olddir(m_dir()); int count(5);
            while( count-->0){
                m_dir = olddir -delta;
                setDir(m_dir);
                double newTs(maximize(skip));
                if( newTs > oldTs ){
                    oldTs=newTs;
                    break;
                }
                delta*=0.5;
                if( verbose() ){ out()<< setw(54) <<setprecision(0) << newTs << " --back off"
                    <<setprecision(3) << std::endl; }

            }

            if( gradmag < 0.1 || delta.mag()< 0.1*sig) break;
        }// iter loop
        if( iter==maxiter){
            if( verbose() ) out() << "   >>>did not converge" << std::endl;
            return 99.;
        }
        if( iter<5 ){
            if(verbose() ) out() << "    *** good fit *** " << std::endl;
        }
        return errorCircle(skip);

}

double UbPointSourceLikelihood::localize(int skip1, int skip2, int itermax, double TSmin )
{
    double sig(99);

    double currentTS(TS());
    if(verbose()) printSpectrum();

    for( int iter(0); iter<itermax; ++iter){
        if( TS()>TSmin) {
            sig = localize(skip1, skip2); // may skip low levels
            if( sig<1) { // good location?
                maximize();
                printSpectrum();
            }
        }
        if( TS() < currentTS+0.1 ) break; // done if didn't improve
        currentTS = TS();
    }
    return sig;
}

