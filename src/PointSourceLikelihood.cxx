/** @file PointSourceLikelihood.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/PointSourceLikelihood.cxx,v 1.35 2008/04/22 17:57:26 mar0 Exp $

*/

#include "pointlike/PointSourceLikelihood.h"

#include "skymaps/DiffuseFunction.h"
#include "skymaps/CompositeSkySpectrum.h"
#ifdef OLD
#include "skymaps/PhotonMap.h"
#else
#include "skymaps/BinnedPhotonData.h"
#endif

#include "embed_python/Module.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>

//#define LEVELS

using namespace astro;
using namespace pointlike;
using skymaps::CompositeSkySpectrum;
using skymaps::DiffuseFunction;


namespace {
    // the scale_factor used: 2.5 degree at level 6, approximately the 68% containment
    double s_TScut(2.);  // only combine energy bands

}

//  ----- static (class) variables -----
skymaps::SkySpectrum* PointSourceLikelihood::s_diffuse(0);

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
    const skymaps::BinnedPhotonData& data,
    std::string name,
    const astro::SkyDir& dir)
    : m_name(name)
    , m_dir(dir)
    , m_out(&std::cout)
    , m_background(0)
{
    
    if( s_diffuse !=0){
        m_background = new skymaps::CompositeSkySpectrum(s_diffuse);
    }else {
        // may not be valid?
        m_background = 0; //new skymaps::CompositeSkySpectrum();
    }
    setup( data, s_minlevel, s_maxlevel);

}


void PointSourceLikelihood::setup(
                                  const skymaps::BinnedPhotonData& data, int minlevel, int maxlevel
                                  )
{

    using skymaps::Band;
    for( skymaps::BinnedPhotonData::const_iterator bit = data.begin(); bit!=data.end(); ++bit){
        const Band& b = *bit;
        // get previous level from nside
        int nside(b.nside()), level(1);
        for( ; level<14;++level){ 
            nside/=2; if( (nside&1)!=0) break;
        }
        if( level<minlevel || level>maxlevel) continue;

        SimpleLikelihood* sl = new SimpleLikelihood(b, m_dir, 
            SimpleLikelihood::defaultUmax() 
            ,m_background
            );
        (*this)[level] = sl;

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

        out() << "level events   signal_fract  TS " << std::right << std::endl;
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
double PointSourceLikelihood::set_tolerance(double tol)
{
    double old(SimpleLikelihood::tolerance());
    SimpleLikelihood::setTolerance(tol);
    return old;
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


