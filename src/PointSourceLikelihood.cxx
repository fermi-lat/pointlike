/** @file PointSourceLikelihood.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/PointSourceLikelihood.cxx,v 1.13 2007/11/04 22:11:32 burnett Exp $

*/

#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/DiffuseFunction.h"

#include "map_tools/PhotonMap.h"
#include "embed_python/Module.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>


using namespace astro;
using namespace pointlike;
namespace {

   //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // from fits 2/26/2006 
    double fit_gamma[]={0,0,0,0,0,2.25,
        2.27,2.22,2.25,2.25,2.29,2.14,2.02,1.87};
#if 0 // old
    double fit_sigma[]={0,0,0,0,0,0,
        0.335,0.319,0.332,0.352,0.397,0.446,0.526,0.657};
#else // from empirical study of a single source 
    // set level 8 -13 using handoff (7/14/07) tweak level 12 from .52 to .62
    double fit_sigma[]={0,0,0,0,0,0.343,
        0.335,0.319, //0.422, 0.9, 0.9, 1.0, 1.1, 1.0 
//    0.42899897,  0.45402442,  0.4742285,   0.60760651,  0.62,  0.94671
    0.431, 0.449, 0.499, 0.566, 0.698, 0.818 // from Marshall
    };

#endif
    // the scale_factor used: 2.5 degree at level 6, approximately the 68% containment
    int base_level(6);
    double scale_factor(int level){return 2.5*pow(2.0, base_level-level)*M_PI/180.;}

    double s_TScut(2.);  // only combine energy bands

    // temporary kluge: energy to assume per level, or band
    // todo: setup more central place for this, worry about event class
    double energy_level(int level){ return 120*pow(2.35, level-6);} 

}

//---    Static (class) variables: defaults here, that can be set by a static function.
// set these from preliminary data above
std::vector<double> PointSourceLikelihood::gamma_level(fit_gamma, fit_gamma+sizeof(fit_gamma)/sizeof(double)); 
std::vector<double> PointSourceLikelihood::sigma_level(fit_sigma, fit_sigma+sizeof(fit_gamma)/sizeof(double));


double PointSourceLikelihood::s_radius(7.0);
int    PointSourceLikelihood::s_minlevel(6);
int    PointSourceLikelihood::s_maxlevel(13);
double PointSourceLikelihood::s_minalpha(0.15);
int    PointSourceLikelihood::s_skip1(1);
int    PointSourceLikelihood::s_skip2(2);
int    PointSourceLikelihood::s_itermax(2);
double PointSourceLikelihood::s_TSmin(5.0);
int    PointSourceLikelihood::s_verbose(0);

void PointSourceLikelihood::setParameters(embed_python::Module& par)
{
    
    par.getValue("pslradius",   s_radius,   s_radius);
    par.getValue("minlevel", s_minlevel, s_minlevel);
    par.getValue("maxlevel", s_maxlevel, s_maxlevel);
    par.getValue("minalpha", s_minalpha, s_minalpha);

    par.getValue("skip1",    s_skip1, s_skip1);
    par.getValue("skip2",    s_skip2, s_skip2);
    par.getValue("itermax",  s_itermax, s_itermax);
    par.getValue("TSmin",    s_TSmin, s_TSmin);
    par.getValue("minlevel", s_minlevel, s_minlevel);
    par.getValue("verbose",  s_verbose, s_verbose);

    par.getValue("tolerance",  SimpleLikelihood::s_tolerance, SimpleLikelihood::s_tolerance);

    // needed by SimpleLikelihood
    par.getValue("umax", SimpleLikelihood::s_defaultUmax, SimpleLikelihood::s_defaultUmax);
    std::string diffusefile;
    par.getValue("diffusefile", diffusefile, "");
    if( ! diffusefile.empty() ) {
        SimpleLikelihood::s_diffuse = new DiffuseFunction(diffusefile);
        std::cout << "Using diffuse definition "<< diffusefile << std::endl; 
    }
}

PointSourceLikelihood::PointSourceLikelihood(
    const map_tools::PhotonMap& data,    
    std::string name,
    const astro::SkyDir& dir)
    : m_energies( data.energyBins() ) // load energies from the data object
    , m_minlevel (data.minLevel() )   // and the minimum level
    , m_nlevels( data.levels() )      // and the number of levels
    , m_name(name)
    , m_dir(dir)
    , m_out(&std::cout)
{
    m_verbose = s_verbose!=0;
    if( s_minlevel< m_minlevel || s_maxlevel>m_minlevel+m_nlevels ){
        throw std::invalid_argument("PointSourceLikelihood: invalid levels for data");
    }
    m_energies.push_back(1e5); // put guard at 100GeV
    setup( data, s_radius, s_minlevel, s_maxlevel);
}


void PointSourceLikelihood::setup(const map_tools::PhotonMap& data,double radius, int minlevel, int maxlevel)
{
    for( int level=minlevel; level<maxlevel+1; ++level){

        // create and fill the vector of data for this level 
        data.extract(  m_dir, radius, m_data_vec[level], -1, level);

        // get PSF parameters from fits
        double gamma( gamma_level[level] ),
            sigma ( scale_factor(level)* sigma_level[level]);

        double emin( m_energies[level-m_minlevel]), emax( m_energies[level-m_minlevel+1]);

        // and create the simple likelihood object
        SimpleLikelihood* sl = new SimpleLikelihood(m_data_vec[level], m_dir, 
            gamma, sigma,
            -1, // background (not used now)
            SimpleLikelihood::s_defaultUmax, 
            emin, emax
            );
        (*this)[level] = sl;

        bool debug_print(false);
        if( debug_print ) { // make table of parameters
            out() << std::setw(6) << level 
                << " " << std::setw(10) << std::left << gamma 
                << std::setw(10)<< std::setprecision(5) << sigma  
                << std::right << std::endl;
        }

    }
}

PointSourceLikelihood::~PointSourceLikelihood()
{
    for( iterator it = begin(); it!=end(); ++it){
        delete it->second;
    }
}

double PointSourceLikelihood::maximize(int skip)
{
    m_TS = 0;
    iterator it = begin();
    for( int i = 0; i< skip; ++it, ++i);

    for( ; it!=end(); ++it){
        SimpleLikelihood& like = *(it->second);
        std::pair<double,double> a(like.maximize());
        if( a.first > s_minalpha ) {
            m_TS+= like.TS();
        }
    }
    return m_TS;
}
#if 0
void PointSourceLikelihood::setBackgroundDensity(const std::vector<double>& density)
{
    std::vector<double>::const_iterator id = density.begin();
    for( iterator it = begin(); it!=end(); ++it, ++id){
        double bk(*id);
        it->second->setBackgroundDensity(bk);
    }
}
#endif
void PointSourceLikelihood::setDir(const astro::SkyDir& dir){
    for( iterator it = begin(); it!=end(); ++it){
        it->second->setDir(dir);
    }
    m_dir = dir;
}

Hep3Vector PointSourceLikelihood::gradient(int skip) const{
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

double PointSourceLikelihood::curvature(int skip) const{
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

void PointSourceLikelihood::printSpectrum()
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

        SimpleLikelihood& levellike = *it->second;
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
        double ts(levellike.TS()); 
        if( a.first > s_minalpha ) {
           m_TS+=ts;
        }

        if( verbose() ){
            out() << setprecision(2) << setw(6)<< a.first<<" +/- "
                << std::setw(4)<< a.second 
                << setw(6)<< setprecision(0)<< ts
                << setw(8) << setprecision(2) << levellike.average_b()
#if 1 // temporary log likelihood
                << setw(10) << setprecision(2) << levellike() 
#endif
                << std::endl;
        }
    }
    if( verbose() ){
        if( s_minalpha>0){
            out() << "\tTS sum  (alpha>"<<s_minalpha<<")  ";
        }else{
            out() << "\tTS sum                            ";
        }
        out() <<  m_TS << std::endl;
    }
}
double PointSourceLikelihood::localize(int skip1, int skip2)
{
    double t(100);
    for( int skip=skip1; skip<=skip2; ++skip){
        t = localize(skip);
        if (t<1) return t;
    }
    return t;
}

double PointSourceLikelihood::localize(int skip)
{
    using std::setw; using std::left; using std::setprecision; using std::right;
    using std::fixed;
    int wd(10), iter(0), maxiter(5);
 
        if( verbose()){
            out() 
                << "      Searching for best position, start at level "<< skip+s_minlevel<<"\n"
                << setw(wd) << left<< "Gradient   " 
                << setw(wd) << left<< "delta  "   
                << setw(wd) << left<< "ra"
                << setw(wd) << left<< "dec "
                << setw(wd) << left<< "error "
                << setw(wd) << left<< "Ts "
                <<std::endl;
        }
        double oldTs( maximize(skip)); // initial (partial) TS
        SkyDir last_dir(dir()); // save current direction

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
                if( verbose() ){ out()<< setw(52) <<setprecision(0) << newTs << " --back off"
                    <<setprecision(3) << std::endl; }

            }

            if( gradmag < 0.1 || delta.mag()< 0.1*sig) break;
        }// iter loop
        if( iter==maxiter){
            if( verbose() ) out() << "   >>>did not converge" << std::endl;
            setDir(last_dir()); // restore position
            return 99.;
        }
        if( iter<5 ){
            if(verbose() ) out() << "    *** good fit *** " << std::endl;
        }
        return errorCircle(skip);

}

double PointSourceLikelihood::localize()
{
    int skip1(s_skip1), skip2(s_skip2), itermax(s_itermax);
    double TSmin(s_TSmin);

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

double PointSourceLikelihood::value(const astro::SkyDir& dir, double energy) const
{
    int level(m_minlevel);
    for(; energy> m_energies[level-m_minlevel] && level<= s_maxlevel; ++level);
    std::map<int, SimpleLikelihood*>::const_iterator it = find(level-1);
    if( it==end() ){
        throw std::invalid_argument("PointSourceLikelihood::value--no fit for the requested energy");
    }
    return it->second->operator()(dir);

}
    ///@brief integral for the energy limits, in the given direction
double PointSourceLikelihood::integral(const astro::SkyDir& dir, double emin, double emax)const
{
    // implement by just finding the right bin
    return value(dir, sqrt(emin*emax) );

void PointSourceLikelihood::recalc(int level) {
    // get PSF parameters from fits
        double gamma( gamma_level[level] ),
            sigma ( scale_factor(level)* sigma_level[level]);

        // and create the simple likelihood object
        SimpleLikelihood* sl = new SimpleLikelihood(m_data_vec[level], m_dir, 
            gamma, sigma,
            -1, // background level?
            SimpleLikelihood::s_defaultUmax, energy_level(level));
        (*this)[level] = sl;
}
