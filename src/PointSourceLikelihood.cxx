/** @file PointSourceLikelihood.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/PointSourceLikelihood.cxx,v 1.32 2008/04/19 23:52:04 burnett Exp $

*/

#include "pointlike/PointSourceLikelihood.h"

#include "skymaps/DiffuseFunction.h"
#include "skymaps/CompositeSkySpectrum.h"
#include "skymaps/PhotonMap.h"

#include "embed_python/Module.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>


using namespace astro;
using namespace pointlike;
using skymaps::CompositeSkySpectrum;
using skymaps::DiffuseFunction;
namespace {
    // the scale_factor used: 2.5 degree at level 6, approximately the 68% containment
    int base_level(6);
    double scale_factor(int level){return 2.5*pow(2.0, base_level-level)*M_PI/180.;}

    double s_TScut(2.);  // only combine energy bands
    // default PSF
    double gamma_list[] ={0,0,0,0,0,
        2.25,  2.27,  2.22,  2.31,  2.30,  2.31,  2.16,  2.19,  2.07};
    double sigma_list[] ={0,0,0,0,0,
        0.343, 0.4199,0.4249 ,0.4202 ,0.4028 ,0.4223 ,0.4438 ,0.5113 ,0.5596 };



}

//  ----- static (class) variables -----
skymaps::SkySpectrum* PointSourceLikelihood::s_diffuse(0);
std::vector<double> PointSourceLikelihood::s_gamma_level(gamma_list,gamma_list+14);
std::vector<double> PointSourceLikelihood::s_sigma_level(sigma_list,sigma_list+14);


int    PointSourceLikelihood::s_minlevel(6);
int    PointSourceLikelihood::s_maxlevel(13);
double PointSourceLikelihood::s_minalpha(0.05);
int    PointSourceLikelihood::s_skip1(1);
int    PointSourceLikelihood::s_skip2(3);
int    PointSourceLikelihood::s_itermax(1);
double PointSourceLikelihood::s_TSmin(5.0);
int    PointSourceLikelihood::s_verbose(0);
double PointSourceLikelihood::s_maxstep(0.25);  // if calculated step is larger then this (deg), abort localization

void PointSourceLikelihood::setParameters(const embed_python::Module& par)
{
    static std::string prefix("PointSourceLikelihood.");

    par.getValue(prefix+"minlevel", s_minlevel, s_minlevel);
    par.getValue(prefix+"maxlevel", s_maxlevel, s_maxlevel);
    par.getValue(prefix+"minalpha", s_minalpha, s_minalpha);

    par.getValue(prefix+"skip1",    s_skip1, s_skip1);
    par.getValue(prefix+"skip2",    s_skip2, s_skip2);
    par.getValue(prefix+"itermax",  s_itermax, s_itermax);
    par.getValue(prefix+"TSmin",    s_TSmin, s_TSmin);
    par.getValue(prefix+"minlevel", s_minlevel, s_minlevel);
    par.getValue(prefix+"verbose",  s_verbose, s_verbose);
    par.getValue(prefix+"maxstep",  s_maxstep, s_maxstep); // override with global
    par.getValue("verbose",  s_verbose, s_verbose); // override with global

    // needed by SimpleLikelihood
    double umax(SimpleLikelihood::defaultUmax());
    par.getValue(prefix+"umax", umax, umax);
    SimpleLikelihood::setDefaultUmax(umax);

    double tolerance(SimpleLikelihood::tolerance());
    par.getValue("Diffuse.tolerance",  tolerance, tolerance);
    SimpleLikelihood::setTolerance(tolerance);

    // load parameters from the setup 
    s_gamma_level.clear(); s_sigma_level.clear();
    par.getList(prefix+"gamma_list", s_gamma_level);
    par.getList(prefix+"sigma_list", s_sigma_level);
    // require that all  levels were set
    int sigsize(s_sigma_level.size());
    if( s_gamma_level.size() !=14 || sigsize !=14){
        throw std::invalid_argument("PointSourceLikelihood::setParameters: gamma or sigma parameter not set properly");
    }

    std::string diffusefile;
    par.getValue("Diffuse.file", diffusefile);
    double exposure(1.0);
    par.getValue("Diffuse.exposure", exposure);
    int interpolate(0);
    par.getValue("interpolate", interpolate, interpolate);
    if( ! diffusefile.empty() ) {

        set_diffuse(new CompositeSkySpectrum(
            new DiffuseFunction(diffusefile, interpolate!=0), exposure) );

        std::cout << "Using diffuse definition "<< diffusefile 
            << " with exposure factor " << exposure << std::endl; 
    }
}


PointSourceLikelihood::PointSourceLikelihood(
    const skymaps::PhotonMap& data,    
    std::string name,
    const astro::SkyDir& dir)
    : m_energies( data.energyBins() ) // load energies from the data object
    , m_minlevel (data.minLevel() )   // and the minimum level
    , m_nlevels( data.levels() )      // and the number of levels
    , m_name(name)
    , m_dir(dir)
    , m_out(&std::cout)
    , m_background(0)
{
    if( s_gamma_level.size()==0){
        s_gamma_level.resize(14,2.2);
        s_sigma_level.resize(14,0.4); 
        std::cerr << "Warning, PointSourceLikelihood: PSF not set up, setting default gamma, sigma to 2.2, 0.4" 
            << std::endl;
    }
    m_verbose = s_verbose!=0;
    if( s_minlevel< m_minlevel || s_maxlevel>m_minlevel+m_nlevels ){
        throw std::invalid_argument("PointSourceLikelihood: invalid levels for data");
    }
    m_energies.push_back(2e5); // put guard at 200GeV

    
    if( s_diffuse !=0){
        m_background = new skymaps::CompositeSkySpectrum(s_diffuse);
    }else {
        // may not be valid?
        m_background = 0; //new skymaps::CompositeSkySpectrum();
    }

    setup( data, s_minlevel, s_maxlevel);

}


void PointSourceLikelihood::setup(const skymaps::PhotonMap& data, int minlevel, int maxlevel)
{
    bool debug_print(false);
    out() << "level   gamma   sigma   roi  pixels" << std::endl;
    for( int level=minlevel; level<maxlevel+1; ++level){


        // get PSF parameters from fits
        double gamma( gamma_level(level) ),
            sigma ( scale_factor(level)* sigma_level(level));

        // and then the radius of the roi for extraction

        double roi_radius( sigma*sqrt(2.*SimpleLikelihood::defaultUmax()) * 180/M_PI);

        // create and fill the vector of data for this level (note fudge) 
        data.extract_level(  m_dir, 1.2* roi_radius, 
            m_data_vec[level], level,  false);
        double emin( m_energies[level-m_minlevel]), emax( m_energies[level-m_minlevel+1]);

        // and create the simple likelihood object
        SimpleLikelihood* sl = new SimpleLikelihood(m_data_vec[level], m_dir, 
            gamma, sigma,
            -1, // background (not used now)
            SimpleLikelihood::defaultUmax(), 
            emin, emax
            ,m_background
            );
        (*this)[level] = sl;

        if( debug_print ) { // make table of parameters
            out() << std::setw(6) << level 
                << " " << std::setw(10) << std::left << gamma 
                << std::setw(10)<< std::setprecision(5) << sigma  
                << std::setw(10)<< std::setprecision(5) << roi_radius  
                << std::setw(10)<< std::setprecision(5) << m_data_vec[level].size()  
                << std::right << std::endl;
        }

    }
}

PointSourceLikelihood::~PointSourceLikelihood()
{
    for( iterator it = begin(); it!=end(); ++it){
        delete it->second;
    }
    delete m_background;
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
#if 0 // obsolete?
void PointSourceLikelihood::setBackgroundDensity(const std::vector<double>& density)
{
    std::vector<double>::const_iterator id = density.begin();
    for( iterator it = begin(); it!=end(); ++it, ++id){
        double bk(*id);
        it->second->setBackgroundDensity(bk);
    }
}
#endif
void PointSourceLikelihood::setDir(const astro::SkyDir& dir, bool subset){
    for( iterator it = begin(); it!=end(); ++it){
        it->second->setDir(dir,subset);
    }
    m_dir = dir;
}

const Hep3Vector& PointSourceLikelihood::gradient(int skip) const{
    m_gradient=Hep3Vector(0);  
    const_iterator it = begin();
    for( int i = 0; i< skip; ++it, ++i);
    for( ; it!=end(); ++it){
        if( it->second->TS()< s_TScut) continue;
        Hep3Vector grad(it->second->gradient());
        double curv(it->second->curvature());
        if( curv > 0 ) m_gradient+= grad;
    }
    return m_gradient;
}

double PointSourceLikelihood::curvature(int skip) const{
    double t(0);
    const_iterator it = begin();
    for( int i = 0; i< skip; ++it, ++i);
    for( ; it!=end(); ++it){
        if( it->second->TS()< s_TScut) continue;
#if 0 // Marshall?
        it->second->gradient();
#endif
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

        out() << "level events   signal_fract  TS " << std::endl;
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
                //out() << "     -    ";
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
            double avb(levellike.average_b());
            out() << setprecision(2) << setw(6)<< a.first<<" +/- "
                << std::setw(4)<< a.second 
                << setw(6)<< setprecision(0)<< ts;
#if 0 // debug output for average background check
            out() << setprecision(2) << std::scientific << " " <<levellike.average_b()<< std::fixed ;
#endif
            out() << std::endl;
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
        if (t<1) break;
    }
    return t;
}

double PointSourceLikelihood::localize(int skip)
{
    using std::setw; using std::left; using std::setprecision; using std::right;
    using std::fixed;
    int wd(10), iter(0), maxiter(20);
    double steplimit(10.0), // in units of sigma
        stepmin(0.1);     // quit if step this small
    double backoff_ratio(0.5); // scale step back if TS does not increase
    int backoff_count(2);

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
    SkyDir last_dir(dir()); // save current direction
    setDir(dir(), true);    // initialize
    double oldTs( maximize(skip)); // initial (partial) TS
    bool backingoff;  // keep track of backing


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
        double step(delta.mag());

        if( verbose() ){
            out() << fixed << setprecision(0)
                <<  setw(wd-2)  << right<< gradmag << "  "
                <<  setprecision(4)
                <<  setw(wd) << left<< step*180/M_PI
                <<  setw(wd) << left<< m_dir.ra()
                <<  setw(wd) << left<< m_dir.dec() 
                <<  setw(wd) << left<< sig*180/M_PI
                <<  setw(wd) << left <<setprecision(1)<< oldTs
                <<  right <<setprecision(3) <<std::endl;
        }
#if 0
        if( step*180/M_PI > s_maxstep) {
            if( verbose() ){ out() << " >>> aborting, attempted step " 
                << (step*180/M_PI) << " deg  greater than limit " 
                << s_maxstep << std::endl;
            }
            return 97.;
            break;
        }
#endif
        // check for too large step, limit to steplimt* sigma
        if( step > steplimit* sig ) {
            delta = steplimit*sig* delta.unit();
            if( verbose() ) out() << setw(52) << "reduced step to "
                <<setprecision(5)<< delta.mag() << std::endl;
        }

        // here decide to back off if likelihood does not increase
        Hep3Vector olddir(m_dir()); int count(backoff_count); 
        backingoff =true;
        while( count-->0){
            m_dir = olddir -delta;
            setDir(m_dir,true);
            double newTs(maximize(skip));
            if( newTs > oldTs-0.01 ){ // allow a little slop
                oldTs=newTs;
                backingoff=false;
                break;
            }
            delta*=backoff_ratio; 
            if( verbose() ){ out()<< setw(56) <<setprecision(1) << newTs << " --back off "
                <<setprecision(5)<< delta.mag() << std::endl; }

        }

        if( gradmag < 0.1 || delta.mag()< stepmin*sig) break;
    }// iter loop
    if( iter==maxiter && ! backingoff ){
        if( verbose() ) out() << "   >>>did not converge" << std::endl;
        setDir(last_dir()); // restore position
        return 99.;
    }
    if(verbose() ) out() << "    *** good fit *** " << std::endl;
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

double PointSourceLikelihood::display(const astro::SkyDir& dir, double energy, int mode) const
{
    int level(m_minlevel);
    for(; energy> m_energies[level-m_minlevel] && level<= s_maxlevel; ++level);
    std::map<int, SimpleLikelihood*>::const_iterator it = find(level-1);
    if( it==end() ){
        throw std::invalid_argument("PointSourceLikelihood::display--no fit for the requested energy");
    }
    return it->second->display(dir, mode);

}


///@brief integral for the energy limits, in the given direction
double PointSourceLikelihood::integral(const astro::SkyDir& dir, double emin, double emax)const
{
    // implement by just finding the right bin
    return value(dir, sqrt(emin*emax) );
}

void PointSourceLikelihood::recalc(int level) {
    // get PSF parameters from fits
    double gamma( gamma_level(level) ),
        sigma ( scale_factor(level)* sigma_level(level));
    find(level)->second->setgamma(gamma);
    find(level)->second->setsigma(sigma);
    find(level)->second->recalc();
}

double PointSourceLikelihood::sigma(int level)const
{
    std::map<int, SimpleLikelihood*>::const_iterator it = find(level);
    if( it==end() ){
        throw std::invalid_argument("PointSourceLikelihood::sigma--no fit for the requested level");
    }
    return it->second->sigma();
}

skymaps::SkySpectrum* PointSourceLikelihood::set_diffuse(skymaps::SkySpectrum* diffuse, double exposure)
{  
    // save current to return
    skymaps::SkySpectrum* ret =   s_diffuse;

    s_diffuse = diffuse;
    
    return ret;
}


void PointSourceLikelihood::addBackgroundPointSource(const PointSourceLikelihood* fit)
{
    if( fit==this){
        throw std::invalid_argument("PointSourceLikelihood::setBackgroundFit: cannot add self as background");
    }
    if( m_background==0){
        throw std::invalid_argument("PointSourceLikelihood::setBackgroundFit: no diffuse background");
    }

    if( s_verbose>0 ) {
        std::cout << "Adding source " << fit->name() << " to background" << std::endl;
    }
    m_background->add(fit);
    setDir(dir()); // recomputes background for each SimpleLikelihood object
}

void PointSourceLikelihood::clearBackgroundPointSource()
{
    if( m_background==0){
        throw std::invalid_argument("PointSourceLikelihood::setBackgroundFit: no diffuse to add to");
    }
    while( m_background->size()>1) {
        m_background->pop_back();
    }
    setDir(dir()); // recomputes background for each SimpleLikelihood object
}

const skymaps::SkySpectrum * PointSourceLikelihood::background()const
{
    return m_background;
}


/// @brief set radius for individual fits
void PointSourceLikelihood::setDefaultUmax(double umax)
{ 
    SimpleLikelihood::setDefaultUmax(umax); 
}


double PointSourceLikelihood::set_gamma_level(int level, double v)
{
    double t = s_gamma_level[level]; s_gamma_level[level]=v; 
    return t;
}

double PointSourceLikelihood::set_sigma_level(int level, double v)
{
    double t = s_sigma_level[level]; s_sigma_level[level]=v; 
    return t;
}



/// @brief get the starting, ending levels used
int PointSourceLikelihood::minlevel(){return s_minlevel;}

int PointSourceLikelihood::maxlevel(){return s_maxlevel;}

void PointSourceLikelihood::set_levels(int minlevel, int maxlevel)
{ s_minlevel= minlevel; s_maxlevel = maxlevel;
}

double PointSourceLikelihood::set_tolerance(double tol)
{
    double old(SimpleLikelihood::tolerance());
    SimpleLikelihood::setTolerance(tol);
    return old;
}

double PointSourceLikelihood::gamma_level(int i)
{
    return s_gamma_level.at(i);
}
double PointSourceLikelihood::sigma_level(int i)
{
    return s_sigma_level.at(i);
}

//=======================================================================
//         PSLdisplay implementation
PSLdisplay::PSLdisplay(const PointSourceLikelihood & psl, int mode)
: m_psl(psl)
, m_mode(mode)
{}

double PSLdisplay::value(const astro::SkyDir& dir, double e)const{
        return m_psl.display(dir, e, m_mode);
    }

    ///@brief integral for the energy limits, in the given direction -- not impleme
double PSLdisplay::integral(const astro::SkyDir& dir, double a, double b)const{
        return value(dir, sqrt(a*b));
    }

std::string PSLdisplay::name()const
{
    static std::string type[]={"density", "data", "background", "fit", "residual"};
    if( m_mode<0 || m_mode>4) return "illegal";
    return m_psl.name()+"--"+type[m_mode];

}


