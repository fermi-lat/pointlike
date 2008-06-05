/** @file PointSourceLikelihood.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/PointSourceLikelihood.cxx,v 1.39 2008/05/27 16:46:41 burnett Exp $

*/

#include "pointlike/PointSourceLikelihood.h"

#include "skymaps/DiffuseFunction.h"
#include "skymaps/CompositeSkySpectrum.h"
#include "skymaps/BinnedPhotonData.h"

#include "embed_python/Module.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>


using namespace pointlike;
using astro::SkyDir;
using skymaps::CompositeSkySpectrum;
using skymaps::BinnedPhotonData;
using skymaps::DiffuseFunction;
using skymaps::Band;



namespace {
    // the scale_factor used: 2.5 degree at level 6, approximately the 68% containment
    double s_TScut(2.);  // only combine energy bands

}

//  ----- static (class) variables -----
skymaps::SkySpectrum* PointSourceLikelihood::s_diffuse(0);

// manage energy range for selection of bands to fit 
double PointSourceLikelihood::s_emin(100.); 
double PointSourceLikelihood::s_emax(1e6);
void PointSourceLikelihood::set_energy_range(double emin, double emax){
    s_emin = emin; s_emax=emax;
}

double PointSourceLikelihood::s_minalpha(0.05);
int    PointSourceLikelihood::s_skip1(1);
int    PointSourceLikelihood::s_skip2(3);
int    PointSourceLikelihood::s_itermax(1);
double PointSourceLikelihood::s_TSmin(5.0);

int    PointSourceLikelihood::s_verbose(0);
void PointSourceLikelihood::set_verbose(bool verbosity){s_verbose=verbosity;}
bool PointSourceLikelihood::verbose(){return s_verbose;}

double PointSourceLikelihood::s_maxstep(0.25);  // if calculated step is larger then this (deg), abort localization

void PointSourceLikelihood::setParameters(const embed_python::Module& par)
{
    static std::string prefix("PointSourceLikelihood.");

    par.getValue(prefix+"emin",     s_emin, s_emin);
    par.getValue(prefix+"emax",     s_emax, s_emax);
    par.getValue(prefix+"minalpha", s_minalpha, s_minalpha);

    par.getValue(prefix+"skip1",    s_skip1, s_skip1);
    par.getValue(prefix+"skip2",    s_skip2, s_skip2);
    par.getValue(prefix+"itermax",  s_itermax, s_itermax);
    par.getValue(prefix+"TSmin",    s_TSmin, s_TSmin);
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
    const BinnedPhotonData& data,
    std::string name,
    const SkyDir& dir)
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
    setup( data);

}


void PointSourceLikelihood::setup( const skymaps::BinnedPhotonData& data )
{

    for( skymaps::BinnedPhotonData::const_iterator bit = data.begin(); bit!=data.end(); ++bit){
        const Band& b = *bit;

        double emin(floor(b.emin()+0.5) ), emax(floor(b.emax()+0.5));
        if( emin < s_emin && emax < s_emin ) continue;
        if( emax > s_emax ) break;

        SimpleLikelihood* sl = new SimpleLikelihood(b, m_dir, 
            SimpleLikelihood::defaultUmax() 
            ,m_background
            );
        this->push_back( sl );
    }
    if( this->empty()){
        throw std::invalid_argument("PointSourceLikelihood::setup: no bands to fit1");
    }
}

PointSourceLikelihood::~PointSourceLikelihood()
{
    for( iterator it = begin(); it!=end(); ++it){
        delete *it;
    }
    delete m_background;
}

double PointSourceLikelihood::maximize()
{
    m_TS = 0;
    iterator it = begin();
    for( ; it!=end(); ++it){
        SimpleLikelihood& like = **it;
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
        (*it)->setDir(dir,subset);
    }
    m_dir = dir;
}

const Hep3Vector& PointSourceLikelihood::gradient() const{
    m_gradient=Hep3Vector(0);  
    const_iterator it = begin();
    for( ; it!=end(); ++it){
        if( (*it)->TS()< s_TScut) continue;
        Hep3Vector grad((*it)->gradient());
        double curv((*it)->curvature());
        if( curv > 0 ) m_gradient+= grad;
    }
    return m_gradient;
}

double PointSourceLikelihood::curvature() const{
    double t(0);
    const_iterator it = begin();
    for( ; it!=end(); ++it){
        if( (*it)->TS()< s_TScut) continue;
#if 0 // Marshall?
        it->second->gradient();
#endif
        double curv((*it)->curvature());
        if( curv>0 )  t+= curv;
    }
    return t;
}

void PointSourceLikelihood::printSpectrum()
{

    using std::setw; using std::left; using std::setprecision; 
    out() << "\nSpectrum of source " << m_name << " at ra, dec=" 
        << setprecision(6) << m_dir.ra() << ", "<< m_dir.dec() << std::endl;

    out() << "  emin eclass events   signal_fract    TS " << std::right << std::endl;

    m_TS =0;
    for( const_iterator it = begin(); it!=end(); ++it){

        SimpleLikelihood& levellike = **it;
        const skymaps::Band& band ( levellike.band() );
 
        double bkg(levellike.background());
        out()  << std::fixed << std::right 
            << setw(7) << static_cast<int>( band.emin()+0.5 )
            << setw(5) << band.event_class()
            << setw(8) << levellike.photons()
            << setw(10);
        if(bkg>=0) {
            out() << setprecision(1) << levellike.background();
        }else{
            //out() << "     -    ";
        }
        if( levellike.photons()==0)  out() << std::endl; 

        if( levellike.photons()==0) {
            continue;
        }
        std::pair<double,double> a(levellike.maximize());
        double ts(levellike.TS()); 
        if( a.first > s_minalpha ) {
            m_TS+=ts;
        }

        double avb(levellike.average_b());
        out() << setprecision(2) << setw(6)<< a.first<<" +/- "
            << setw(4)<< a.second 
            << setw(6)<< setprecision(0)<< ts;
#if 0 // debug output for average background check
        out() << setprecision(2) << std::scientific << " " <<levellike.average_b()<< std::fixed ;
#endif
        out() << std::endl;
    }
    if( s_minalpha>0.1){
        out() << "\tTS sum  (alpha>"<<s_minalpha<<")  ";
    }else{
        out() << setw(30) << "sum  ";
    }
    out() << setw(14) << m_TS << std::endl;
}

std::vector<double> PointSourceLikelihood::energyList()const
{
    
    std::vector<double> energies;
    if( size()>0) {
        const_iterator it (begin());
        for( ; it!= end(); ++it){
            double emin((*it)->band().emin());
            if( energies.size()==0 || energies.back()!= emin){
                energies.push_back(emin);
            }
        }
        energies.push_back( back()->band().emax() );
    }
    return energies;

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
        stepmin(0.25);     // quit if step this small, in units of sigma
    double backoff_ratio(0.5); // scale step back if TS does not increase
    int backoff_count(2);

    if( verbose()){
        out() 
            << "      Searching for best position, start at band "<< skip <<"\n"
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
    double oldTs( maximize()); // initial (partial) TS
    bool backingoff;  // keep track of backing


    for( ; iter<maxiter; ++iter){
        Hep3Vector grad( gradient() );
        double     curv( curvature() );

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
            double newTs(maximize());
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
    return errorCircle();

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
            }
        }
        if( TS() < currentTS+0.1 ) break; // done if didn't improve
        currentTS = TS();
    }
    return sig;
}

double PointSourceLikelihood::value(const astro::SkyDir& dir, double energy) const
{
    double result(0);
    const_iterator it = begin();
    for( ; it!=end(); ++it){
        const Band& band ( (*it)->band() );
        if( energy >= band.emin() && energy < band.emax() ){
            result += (**it)(dir);
        }
    }
    return result;

}

double PointSourceLikelihood::display(const astro::SkyDir& dir, double energy, int mode) const
{
    const_iterator it = begin();
    for( ; it!=end(); ++it){
        const Band& band ((*it)->band());
        if( energy >= band.emin() && energy < band.emax() )break;
    }
    if( it==end() ){
        throw std::invalid_argument("PointSourceLikelihood::display--no fit for the requested energy");
    }
    return (*it)->display(dir, mode);
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


