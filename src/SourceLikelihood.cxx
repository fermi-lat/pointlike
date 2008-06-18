/** @file SourceLikelihood.cxx

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SourceLikelihood.cxx,v 1.1 2008/06/18 01:19:24 funk Exp $

*/
#define USE_GRADIENT

#include "pointlike/SourceLikelihood.h"
#include "skymaps/DiffuseFunction.h"
#include "skymaps/CompositeSkySpectrum.h"
#include "skymaps/BinnedPhotonData.h"

#include "embed_python/Module.h"
#include "astro/SkyDir.h"

#ifdef USE_MINPACK
#include "minpack/minpack.h"
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include <TMinuit.h>

//to be tested....
using namespace pointlike;

namespace {

  int gFitCounter = 0;
  pointlike::SourceLikelihood* gSourcePointer = NULL;
  CLHEP::Hep3Vector gFitStartDir(0.,0.,0.);
  CLHEP::Hep3Vector gFitDeltaX(0.,0.,0.);
  CLHEP::Hep3Vector gFitDeltaY(0.,0.,0.);
  double gTSvalue=0;
  
  void minuit_likelihood_wrapper(Int_t &npar, Double_t* derivative, 
				 Double_t & loglike, Double_t par[], Int_t iflag){
    if(!gSourcePointer) 
      throw std::invalid_argument("SourceLikelihood::gSourcePointer not set");
    
    double x = par[0];
    double y = par[1];
    //       int skip = gSourcePointer->minuitLevelsToSkip();
    // set source parameters     
    std::vector<double> src_par(npar-2,0);
    for (int i=0;i<npar-2;i++) src_par[i]=par[i+2];
    
    // calculate new direction vector from fit parameters (delta(theta),delta(phi))
    
    CLHEP::Hep3Vector newDir = gFitStartDir+ x* gFitDeltaX + y* gFitDeltaY;
    
    //      std::cout<<"new dir: "<<newDir<<" "<<newDir.mag()<<std::endl;
    newDir.setMag(1.);
    
    astro::SkyDir dirNew(newDir);
    //      std::cout<<"new dir: dec="<<dirNew.dec()<<" ra="<<dirNew.ra()<<std::endl;
    
    //      astro::SkyDir dirNew = astro::SkyDir(x, y, astro::SkyDir::EQUATORIAL);
    
    gSourcePointer->setDir(dirNew,src_par,false);
    
    //       double  TS   = gSourcePointer->maximize(skip);
    double  TS   = gSourcePointer->maximize();
    
// # if 0
//     double TS=0;
//     for( pointlike::SourceLikelihood::iterator it = gSourcePointer->begin(); 
// 	 it!= gSourcePointer->end(); ++it){
//       pointlike::ExtendedLikelihood& like = *(it->second);
//       TS+= like.TS();
//     };
//     if(TS<1e-10 || fabs(gTSvalue-TS)/TS > 0.00 ) {
//       std::cout<<"maximize."<<" "<< fabs(gTSvalue-TS)/TS <<" "<<TS<<" "<<gTSvalue<<std::endl;
//       //         TS   = gSourcePointer->maximize(skip);
//       TS   = gSourcePointer->maximize();
//       gTSvalue = TS;
//     };
// #endif
    loglike = -0.5*TS; //m_loglike;
    
#ifdef USE_GRADIENT
    //       const std::vector<double>& gradient = gSourcePointer->gradient(skip);
    const std::vector<double>& gradient = gSourcePointer->gradient();
    for (int i=0;i<npar;i++) derivative[i]=gradient[i];
#endif
    
    std::cout << "**** Iteration "<<gFitCounter<<" **** Testing parameters:"
	      << std::fixed<< std::setprecision(4) 
	      << " delta(X)=" << (par[0]*180./M_PI)
	      <<" deg\tdelta(Y)="<<(par[1]*180./M_PI)<<" deg"<< std::scientific;
    for (int i=2;i<npar;i++)  std::cout <<"\tp"<<i<<"="<<par[i];
    //       std::cout << "\tTS: " << TS <<" "<<skip<<std::endl;
    std::cout << "\tTS: " << TS <<" "<<std::endl;
    ++gFitCounter;
    
  }
/// minpack likelihood wrapper
 
  int minpack_likelihood_wrapper (void *p, int m, int npar, const double *par, 
				  double *fvec,int iflag ) {
    
    if(!gSourcePointer) 
      throw std::invalid_argument("SourceLikelihood::gSourcePointer not set");
    
    double x = par[0]/180.*M_PI;
    double y = par[1]/180.*M_PI;
    //       int skip = gSourcePointer->minuitLevelsToSkip();
    // set source parameters     
    std::vector<double> src_par(npar-2,0);
    for (int i=0;i<npar-2;i++) src_par[i]=par[i+2]/180.*M_PI;
    
    // calculate new direction vector from fit parameters (delta(theta),delta(phi))
    
    CLHEP::Hep3Vector newDir = gFitStartDir+ x* gFitDeltaX + y* gFitDeltaY;
    newDir.setMag(1.);
    astro::SkyDir dirNew(newDir);
    gSourcePointer->setDir(dirNew,src_par,false);
    
    //       double TS = gSourcePointer->maximize(skip);
    double TS = gSourcePointer->maximize();
    
    int im =0;
    double llike=0;
//     for( pointlike::SourceLikelihood::iterator it = gSourcePointer->begin(); 
// 	 it!= gSourcePointer->end(); ++it){
//       pointlike::ExtendedLikelihood& like = *(it->second);
//       const std::vector<double>& gradient = like.gradient(gFitDeltaX,gFitDeltaY);
//       const std::vector<double>& sqrtLike = like.residual();
//       for(int k=0;k<sqrtLike.size();k++) {
// 	fvec[im++]=sqrtLike[k];
// 	llike+=sqrtLike[k]*sqrtLike[k];
// 	if(im==m) throw std::runtime_error("minpack vector size too small.");
//       };
//     };
    for (int k=im;k<m;k++) fvec[k]=0.;
    
    std::cout << "**** Iteration "<<gFitCounter<<" **** Testing parameters:"
	      << std::fixed<< std::setprecision(4) 
	      << " delta(X)=" << par[0]<<" deg\tdelta(Y)="<<par[1] <<" deg"
	      << std::scientific;

    for (int i=2;i<npar;i++)  std::cout <<"\tp"<<i<<"="<<par[i];
    std::cout << "\tTS: " << TS <<" im="<<im<<" llike="<<(2*llike)<<std::endl;
    ++gFitCounter;
    
    return iflag;
    
  }
  
}

namespace {
  
  // the scale_factor used: 2.5 degree at level 6, approximately the 68% containment
  
  // SF: Apparently no more used ...
  //     int base_level(6);
  //     double scale_factor(int level){return 2.5*pow(2.0, base_level-level)*M_PI/180.;}
  
  double s_TScut(2.);  // only combine energy bands
}
// SF: also no more used ...
//     // default PSF
//     double gamma_list[] ={0,0,0,0,0,
//         2.25,  2.27,  2.22,  2.31,  2.30,  2.31,  2.16,  2.19,  2.07};
//     double sigma_list[] ={0,0,0,0,0,
//         0.343, 0.4199,0.4249 ,0.4202 ,0.4028 ,0.4223 ,0.4438 ,0.5113 ,0.5596 };

skymaps::SkySpectrum* SourceLikelihood::s_diffuse(0);

// manage energy range for selection of bands to fit 
double SourceLikelihood::s_emin(100.); 
double SourceLikelihood::s_emax(1e6);
void SourceLikelihood::set_energy_range(double emin, double emax){
  s_emin = emin; s_emax=emax;
}

// //  ----- static (class) variables -----
// std::vector<double> SourceLikelihood::s_gamma_level(gamma_list,gamma_list+14);
// std::vector<double> SourceLikelihood::s_sigma_level(sigma_list,sigma_list+14);


// double pointlike::SourceLikelihood::s_radius(7.0);
// int    pointlike::SourceLikelihood::s_minlevel(6);
// int    pointlike::SourceLikelihood::s_maxlevel(13);
double SourceLikelihood::s_minalpha(0.05);
int    SourceLikelihood::s_skip1(1);
int    SourceLikelihood::s_skip2(2);
int    SourceLikelihood::s_itermax(1);
double SourceLikelihood::s_TSmin(5.0);
int    SourceLikelihood::s_useMinuit(0);
int    SourceLikelihood::s_minuitLevelSkip(0);
int    SourceLikelihood::s_verbose(0);
double SourceLikelihood::s_maxstep(0.25);  
// if calculated step is larger then this (deg), abort localization

void SourceLikelihood::set_verbose(bool verbosity){s_verbose=verbosity;}
bool SourceLikelihood::verbose(){return s_verbose;}

void SourceLikelihood::setParameters(const embed_python::Module& par){
  static std::string prefix("SourceLikelihood.");
  
  par.getValue(prefix+"emin",     s_emin, s_emin);
  par.getValue(prefix+"emax",     s_emax, s_emax);
  //     par.getValue(prefix+"pslradius",   s_radius,   s_radius);
  //     par.getValue(prefix+"minlevel", s_minlevel, s_minlevel);
  //     par.getValue(prefix+"maxlevel", s_maxlevel, s_maxlevel);
  par.getValue(prefix+"minalpha", s_minalpha, s_minalpha);
  
  par.getValue(prefix+"skip1",    s_skip1, s_skip1);
  par.getValue(prefix+"skip2",    s_skip2, s_skip2);
  par.getValue(prefix+"itermax",  s_itermax, s_itermax);
  par.getValue(prefix+"TSmin",    s_TSmin, s_TSmin);
  par.getValue(prefix+"UseMinuit",s_useMinuit, s_useMinuit);
  par.getValue(prefix+"verbose",  s_verbose, s_verbose);
  par.getValue(prefix+"maxstep",  s_maxstep, s_maxstep); // override with global
  
  // needed by ExtendedLikelihood
  double umax(pointlike::ExtendedLikelihood::defaultUmax());
  par.getValue(prefix+"umax", umax, umax);
  pointlike::ExtendedLikelihood::setDefaultUmax(umax);
  
  double tolerance(pointlike::ExtendedLikelihood::tolerance());
  par.getValue("Diffuse.tolerance",  tolerance, tolerance);
  pointlike::ExtendedLikelihood::setTolerance(tolerance);
  
  // SF: no more needed
  //     // load parameters from the setup .3
  //     s_gamma_level.clear(); s_sigma_level.clear();
  //     par.getList(prefix+"gamma_list", s_gamma_level);
  //     par.getList(prefix+"sigma_list", s_sigma_level);
  //     // require that all  levels were set
  //     if( s_gamma_level.size() !=14 || s_sigma_level.size() !=14){
  //         throw std::invalid_argument("SourceLikelihood::setParameters: gamma or sigma parameter not set properly");
  //     }
  
  std::string diffusefile;
  par.getValue("Diffuse.file", diffusefile);
  double exposure(1.0);
  par.getValue("Diffuse.exposure", exposure);
  int interpolate(0);
  par.getValue("interpolate", interpolate, interpolate);
  if( ! diffusefile.empty() ) {
    //         ExtendedLikelihood::setDiffuse(new CompositeSkySpectrum(
    set_diffuse(new skymaps::CompositeSkySpectrum
		(new skymaps::DiffuseFunction(diffusefile, interpolate!=0), exposure) );
    
    std::cout << "Using diffuse definition "<< diffusefile 
	      << " with exposure factor " << exposure << std::endl; 
  }
}

SourceLikelihood::SourceLikelihood(const skymaps::BinnedPhotonData& data,    
				   std::string name,
				   const astro::SkyDir& dir,
				   const std::string type,
				   const std::vector<double>& src_param)
  //     : m_energies( data.energyBins() ) // load energies from the data object
  //     , m_minlevel (data.minLevel() )   // and the minimum level
  //     , m_nlevels( data.levels() )      // and the number of levels
  : m_name(name)
  , m_dir(dir)
  , m_type(type)
  , m_npar(src_param.size())
  , m_sourceParameters(src_param)
  , m_sourceParErrors(src_param.size(),0)
  , m_out(&std::cout)
  , m_background(0){
  //     if( s_gamma_level.size()==0){
  //         s_gamma_level.resize(14,2.2);
  //         s_sigma_level.resize(14,0.4); 
  //         std::cerr << "Warning, SourceLikelihood: PSF not set up, setting default gamma, sigma to 2.2, 0.4" 
  //             << std::endl;
  //     }
  
  if( s_diffuse !=0){
    m_background = new skymaps::CompositeSkySpectrum(s_diffuse);
  }else {
    // may not be valid?
    m_background = 0; //new skymaps::CompositeSkySpectrum();
  }
  
  //   m_verbose = s_verbose!=0;
  //   if( s_minlevel< m_minlevel || s_maxlevel>m_minlevel+m_nlevels ){
  //     throw std::invalid_argument("SourceLikelihood: invalid levels for data");
  //   }
  //     m_energies.push_back(1e5); // put guard at 100GeV
  //   setup( data, s_radius, s_minlevel, s_maxlevel);
  setup( data);
}

// void SourceLikelihood::setup(const skymaps::PhotonMap& data,double radius, int minlevel, int maxlevel)

void SourceLikelihood::setup(const skymaps::BinnedPhotonData& data){
  
  for( skymaps::BinnedPhotonData::const_iterator bit = 
	 data.begin(); bit!=data.end(); ++bit){
    const skymaps::Band& b = *bit;
    
    double emin(floor(b.emin()+0.5) ), emax(floor(b.emax()+0.5));
    if( emin < s_emin && emax < s_emin ) continue;
    if( emax > s_emax ) break;
    
    
    //         // get PSF parameters from fits
    //         double gamma( gamma_level(level) ),
    //             sigma ( scale_factor(level)* sigma_level(level));
    
    //         double emin( m_energies[level-m_minlevel]), emax( m_energies[level-m_minlevel+1]);
    
    //	std::cout << "   Creating ExtendedLikelihood object for level: "
    //		  << level << " " << m_minlevel 
    //		  << " at Energy " << emin << "GeV to " << emax << "GeV" << std::endl;
    // and create the simple likelihood object
    //         ExtendedLikelihood* sl = new ExtendedLikelihood(m_data_vec[level], m_dir, m_type,m_sourceParameters,
    // 							gamma, sigma,
    // 							-1, // background (not used now)
    // 							ExtendedLikelihood::defaultUmax(), 
    // 							emin, emax);
    pointlike::ExtendedLikelihood* sl 
      = new pointlike::ExtendedLikelihood(b, m_dir, 
					  m_type,m_sourceParameters,
					  pointlike::ExtendedLikelihood::defaultUmax(), 
					  m_background);
    this->push_back(sl);
    //         (*this)[level] = sl;
    
    //         bool debug_print(false);
    //         if( debug_print ) { // make table of parameters
    //             out() << std::setw(6) << level 
    //                 << " " << std::setw(10) << std::left << gamma 
    //                 << std::setw(10)<< std::setprecision(5) << sigma  
    //                 << std::right << std::endl;
    //         }
    
  }
  if( this->empty()){
    throw std::invalid_argument("SourceLikelihood::setup: no bands to fit1");
  }
  
}

pointlike::SourceLikelihood::~SourceLikelihood()
{
  for( iterator it = begin(); it!=end(); ++it){
    delete *it;
  }
  delete m_background;
}

double pointlike::SourceLikelihood::maximize()
{
  
#if 0
  std::cout << "**** Testing position: " 
	    << std::setw(10)<< std::setprecision(5)
	    << m_dir.ra() << " " << m_dir.dec() 
	    << std::endl; 
  
  std::cout << std::left << std::setw(20) 
	    <<"  iteration" << " alpha   alpha'   alpha''   loglike\n";
#endif
  
  m_TS = 0;
  m_loglike = 0;
  int photons=0;
  iterator it = begin();
  for( ; it!=end(); ++it){
    pointlike::ExtendedLikelihood& like = **it;
    std::pair<double,double> a(like.maximize());
#if 0	
    std::cout << "  --  Summary ------- Level: " << it->first  
	      << " (" << like.GetEmin() << ", " << like.GetEmax() << ") " 
	      << " alpha: " << a.first << " +- " << a.second << std::endl;
#endif        
    if( a.first > s_minalpha ) {
      m_loglike += logL();
      m_TS+= like.TS();
      photons+=like.photons();
    }	
  }
  
#if 0
  std::cout << "**** Iteration: " << gFitCounter 
	    << "**** Position: " 
	    << std::setw(10)<< std::setprecision(5)
	    << m_dir.ra() << " " << m_dir.dec() 
	    << " - TS: " << m_TS <<" photons:" << photons <<" b="<<avb<<" u="
	    <<avu<<" F="<<avf<<std::endl;
#endif
  
  return m_TS;
}

#if 0 // obsolete?
void pointlike::SourceLikelihood::setBackgroundDensity(const std::vector<double>& density)
{
  std::vector<double>::const_iterator id = density.begin();
  for( iterator it = begin(); it!=end(); ++it, ++id){
    double bk(*id);
    it->second->setBackgroundDensity(bk);
  }
}
#endif

void pointlike::SourceLikelihood::setDir(const astro::SkyDir& dir, bool subset){
  for( iterator it = begin(); it!=end(); ++it){
    (*it)->setDir(dir,subset);
  }
  m_dir = dir;
}

void pointlike::SourceLikelihood::setDir(const astro::SkyDir& dir, 
					 const std::vector<double>& srcparam, 
					 bool subset){
  for( iterator it = begin(); it!=end(); ++it){
    (*it)->setDir(dir,srcparam,subset);
  }
  m_dir = dir;
}

const std::vector<double> pointlike::SourceLikelihood::gradient() const{
  std::vector<double> gradient(2+m_npar,0);  
  const_iterator it = begin();
  for( ; it!=end(); ++it){
    //        if( (*it)->TS()< s_TScut) continue;
    const std::vector<double>& level_grad=(*it)->gradient(gFitDeltaX,gFitDeltaY);
    for(int j =0;j<gradient.size(); j++) gradient[j]+= level_grad[j];
  }
  //   for(unsigned int k=0;k<gradient.size();k++) std::cout<<"gradient("<<k<<")="<< gradient[k]<<std::endl;
  return gradient;
}


const Hep3Vector& pointlike::SourceLikelihood::ps_gradient() const{
  m_gradient=Hep3Vector(0);  
  const_iterator it = begin();
  for( ; it!=end(); ++it){
    if( (*it)->TS()< s_TScut) continue;
    Hep3Vector grad((*it)->ps_gradient());
    double curv((*it)->ps_curvature());
    if( curv > 0 ) m_gradient+= grad;
  }
  return m_gradient;
}


double pointlike::SourceLikelihood::ps_curvature() const{
  double t(0);
  const_iterator it = begin();
  for( ; it!=end(); ++it){
    if( (*it)->TS()< s_TScut) continue;
#if 0 // Marshall?
    it->second->ps_gradient();
#endif
    double curv((*it)->ps_curvature());
    if( curv>0 )  t+= curv;
  }
  return t;
}

void pointlike::SourceLikelihood::printSpectrum()
{
  
  using std::setw; using std::left; using std::setprecision; 
  out() << "\nSpectrum of source " << m_name << " at ra, dec=" 
        << setprecision(6) << m_dir.ra() << ", "<< m_dir.dec() << std::endl;
  
  out() << "  emin eclass events   signal_fract    TS " << std::right << std::endl;
  
  m_TS =0;
  for( const_iterator it = begin(); it!=end(); ++it){
    
    pointlike::ExtendedLikelihood& levellike = **it;
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
      // 	  out() << "     -    ";
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
    out() << setprecision(2) << std::scientific << " " 
	  <<levellike.average_b()<< std::fixed ;
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

std::vector<double> pointlike::SourceLikelihood::energyList() const {
  
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

double pointlike::SourceLikelihood::localize(int skip1, int skip2){
  double t(100);
  gFitCounter = 0;
  for( int skip=skip1; skip<=skip2; ++skip){
    if (s_useMinuit) 
#ifdef USE_MINPACK
      t = localizeMinpack(skip);
#else
    t = localizeMinuit(skip);
#endif
    else 
      t = localize(skip);
    if (t<1) break;
  }
  return t;
}
 
 
#ifdef USE_MINPACK
double pointlike::SourceLikelihood::localizeMinpack(int skip){
  using std::setw; using std::left; using std::setprecision; using std::right;
  using std::fixed;
  int wd(10);
  
  astro::SkyDir last_dir(dir()); // save current direction
  setDir(dir(), m_sourceParameters, false);    // initialize
  gFitStartDir = dir()();
  
  gFitDeltaX = CLHEP::Hep3Vector(gFitStartDir.y(),-gFitStartDir.x(),0.);
  if( gFitDeltaX.mag2()<1e-10) 
    gFitDeltaX = gFitStartDir.z()>0 ? Hep3Vector(1.,0.,0.) : Hep3Vector(-1.,0.,0.);
  gFitDeltaX.setMag(1.);	     
  gFitDeltaY = gFitStartDir.cross(gFitDeltaX);
  gFitDeltaY.setMag(1.);
  
  std::cout<<"Fit start direction  : "<<gFitStartDir<<" "<<gFitStartDir.mag()<<std::endl;
  std::cout<<"Fit coordinate system: x="<<gFitDeltaX<<" y="<<gFitDeltaY<<std::endl;
  
  int npar = 2 + m_sourceParameters.size();
  
  minpack::Minpack& mini(minpack::Minpack::Minimizer());
  mini.init(20000,npar,1e-7);
  
  std::cout << "Setting likelihood for "<<npar<<" parameters : " << this << std::endl;
  
  // Set the pointer for access in the minuit function
  gSourcePointer = this;
  mini.setResidualFunc(minpack_likelihood_wrapper);
  
  setMinuitLevelSkip(skip);
  
  std::vector<double> parameter(npar);
  parameter[0]=0;
  parameter[1]=0;
  for (int i=2; i<npar; i++) parameter[i] = m_sourceParameters[i-2]*180./M_PI;
  mini.setX(parameter);
  
  mini.minimize();
  
  std::vector<double> result = mini.getX(); 
  
  for (int i=2; i<npar; i++) {
    m_sourceParameters[i-2]=result[i];
    m_sourceParErrors[i-2]=0.;
  };
  
  CLHEP::Hep3Vector newDir = gFitStartDir+ result[0]* gFitDeltaX + result[1]* gFitDeltaY;
  newDir.setMag(1.);
  astro::SkyDir dirNew(newDir);
  m_errorX=0;
  m_errorY=0;
  
  setDir(dirNew,m_sourceParameters,false);
  
  return sqrt(m_errorX*m_errorX + m_errorY*m_errorY);
}  

#endif


double pointlike::SourceLikelihood::localizeMinuit(int skip)
{
  using std::setw; using std::left; using std::setprecision; using std::right;
  using std::fixed;
  int wd(10);
  
  if( verbose()){
    out() 
      << "      Searching for best position, start at level "
      // 	    << skip+s_minlevel<<"\n"
      << setw(wd) << left<< "Gradient   " 
      << setw(wd) << left<< "delta  "   
      << setw(wd) << left<< "ra"
      << setw(wd) << left<< "dec "
      << setw(wd) << left<< "error "
      << setw(wd) << left<< "Ts "
      <<std::endl;
  }
  
  astro::SkyDir last_dir(dir()); // save current direction
  setDir(dir(), m_sourceParameters, false);    // initialize
  
  gFitStartDir = dir()();
  
  std::cout<<"Fit start direction: "<<gFitStartDir<<" "<<gFitStartDir.mag()<<std::endl;
  
  gFitDeltaX = CLHEP::Hep3Vector(gFitStartDir.y(),-gFitStartDir.x(),0.);
  if( gFitDeltaX.mag2()<1e-10) 
    gFitDeltaX = gFitStartDir.z()>0 ? Hep3Vector(1.,0.,0.) : Hep3Vector(-1.,0.,0.);
  gFitDeltaX.setMag(1.);	     
  gFitDeltaY = gFitStartDir.cross(gFitDeltaX);
  gFitDeltaY.setMag(1.);
  
  std::cout<<"Fit coordinate system: x="<<gFitDeltaX<<" y="<<gFitDeltaY<<std::endl;
  
  int npar = 2 + m_sourceParameters.size();
  //    int npar = 3;
  TMinuit gMinuit(npar);
  std::cout << "Setting likelihood for "<<npar<<" parameters : " << this << std::endl;
  
  // Set the pointer for access in the minuit function
  gSourcePointer = this;
  gMinuit.SetFCN(minuit_likelihood_wrapper);
#ifndef WIN32
  double par[npar];
  double stepSize[npar];
  double minVal[npar];
  double maxVal[npar];
  std::string parName[npar];
#else
  std::vector<double> par(npar), stepSize(npar),minVal(npar),maxVal(npar);
  std::vector<std::string> parName(npar);
#endif  
  par[0] = 0.;            par[1] = 0.;
  stepSize[0] = 0.01;     stepSize[1] = 0.01;
  minVal[0] = -0.5;       minVal[1] = -0.5;
  maxVal[0] = 0.5;        maxVal[1] = 0.5;
  
  parName[0] = std::string("Delta(theta)"); 
  parName[1] = std::string("Delta(phi)"); 
  
  for(int i=2; i<npar; i++){ 
    par[i]      = m_sourceParameters[i-2];
    stepSize[i] = 0.01;
    minVal[i]   = 0;
    maxVal[i]   = 0.5;
    std::stringstream nameStream(parName[i]);
    nameStream<<"srcparam("<<i-2<<")"; 
    parName[i]=nameStream.str();
  };   
  
  setMinuitLevelSkip(skip);
  
  for (int i = 0; i < npar; ++i)
    gMinuit.DefineParameter(i, parName[i].c_str(), par[i], 
			    stepSize[i], minVal[i], maxVal[i]);
  
  Int_t ierflag = 0; // the minuit output flag. Can be queried after each command
  // = 0: command executed normally
  //   1: command is blank, ignored
  //   2: command line unreadable, ignored
  //   3: unknown command, ignored
  //   4: abnormal termination (e.g., MIGRAD not converged)
  //   9: reserved
  //   10: END commandw
  //   11: EXIT or STOP command
  //   12: RETURN command
  
  Double_t arglist[2];
  arglist[0] = 0.5; 
  int nargs = 1;
  // up. Defines parameter errors
  // Minuit defines parameter errors as the change
  // in parameter value required to change the function
  // value by up. For negative log-like: up = 0.5
  gMinuit.mnexcm("SET ERR", arglist, nargs, ierflag);
  
#ifdef USE_GRADIENT
  nargs=1; arglist[0] = 1; 
  gMinuit.mnexcm("SET GRA", arglist, nargs, ierflag);
  gMinuit.SetPrintLevel(0);
#endif
  
  nargs=1; arglist[0] = 1; 
  gMinuit.mnexcm("SET STR", arglist, nargs, ierflag);
  gMinuit.SetPrintLevel(0);
  
  // SF: eventually read itermax from the config-file
  //     arglist[0] = itermax; 
  arglist[0] = 500; 
  // approximate maximum number of function calls
  // even if the minimisation hanot converged
  arglist[1] = 1.; 
  // tolerance (in units of 0.001*UP)
  // minimisation will stop if estimated vertical
  // distance to minimum is less than 0.001*tolerance*UP
  nargs = 2;
  gMinuit.mnexcm("MIGRAD", arglist, nargs, ierflag);
  if (ierflag == 4) {
    gMinuit.mnexcm("HESSE", arglist, nargs, ierflag);
  };
  
  double x = 0;
  double y = 0;
  gMinuit.GetParameter(0, x, m_errorX);
  gMinuit.GetParameter(1, y, m_errorY);
  for (int i=2; i<npar; i++) gMinuit.GetParameter(i, m_sourceParameters[i-2],m_sourceParErrors[i-2]);
  
  CLHEP::Hep3Vector newDir = gFitStartDir+ x* gFitDeltaX + y* gFitDeltaY;
  newDir.setMag(1.);
  astro::SkyDir dirNew(newDir);
  
  //    astro::SkyDir dirNew = astro::SkyDir(x, y, astro::SkyDir::EQUATORIAL);
  setDir(dirNew,m_sourceParameters,false);
  
  if (ierflag == 4) {
    // fitting did not converge
    std::cerr<<"WARNING: Minuit returned ierflag=4: Fit did not converge."<<std::endl;
    //      setDir(last_dir()); // restore position  
    return -1;    
  }
  
  return sqrt(m_errorX*m_errorX + m_errorY*m_errorY);
}  


double pointlike::SourceLikelihood::localize(int skip)
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
      << "      Searching for best position, start at band "<< skip <<"\n"
      << setw(wd) << left<< "Gradient   " 
      << setw(wd) << left<< "delta  "   
      << setw(wd) << left<< "ra"
      << setw(wd) << left<< "dec "
      << setw(wd) << left<< "error "
      << setw(wd) << left<< "Ts "
      <<std::endl;
  }
  
  astro::SkyDir last_dir(dir()); // save current direction
  setDir(dir(), true);    // initialize
  double oldTs( maximize()); // initial (partial) TS
  bool backingoff;  // keep track of backing
  
  
  for( ; iter<maxiter; ++iter){
    Hep3Vector grad( ps_gradient() );
    double     curv( ps_curvature() );
    
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

double pointlike::SourceLikelihood::localize()
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
	//	      int style = useMinuit();
	// Final maximisation with Toby's code ...
	//	      useMinuit() = 0;
	maximize();
	printSpectrum();
	///	      useMinuit() = style;
      }
    }
    if( TS() < currentTS+0.1 ) break; // done if didn't improve
    currentTS = TS();
  }
  return sig;
}

double pointlike::SourceLikelihood::value(const astro::SkyDir& dir, double energy) const
{
  double result(0);
  const_iterator it = begin();
  for( ; it!=end(); ++it){
    const skymaps::Band& band ( (*it)->band() );
    if( energy >= band.emin() && energy < band.emax() ){
      result += (**it)(dir);
    }
  }
  return result;
  
}

double pointlike::SourceLikelihood::display(const astro::SkyDir& dir, double energy, int mode) const
{
    const_iterator it = begin();
    for( ; it!=end(); ++it){
        const skymaps::Band& band ((*it)->band());
        if( energy >= band.emin() && energy < band.emax() )break;
    }
    if( it==end() ){
        throw std::invalid_argument("SourceLikelihood::display--no fit for the requested energy");
    }
    return (*it)->display(dir, mode);
}

///@brief integral for the energy limits, in the given direction
double pointlike::SourceLikelihood::integral(const astro::SkyDir& dir, double emin, double emax)const
{
  // implement by just finding the right bin
  return value(dir, sqrt(emin*emax) );
}

// SF: no more needed
// void SourceLikelihood::recalc(int level) {
//     // get PSF parameters from fits
//     double gamma( gamma_level(level) ),
//         sigma ( scale_factor(level)* sigma_level(level));
//     find(level)->second->setgamma(gamma);
//     find(level)->second->setsigma(sigma);
//     find(level)->second->recalc();
// }

// double SourceLikelihood::sigma(int level)const
// {
//     std::map<int, ExtendedLikelihood*>::const_iterator it = find(level);
//     if( it==end() ){
//         throw std::invalid_argument("SourceLikelihood::sigma--no fit for the requested level");
//     }
//     return it->second->sigma();
// }

skymaps::SkySpectrum* SourceLikelihood::set_diffuse(skymaps::SkySpectrum* diffuse, 
						    double exposure)
{  
  skymaps::SkySpectrum* ret = s_diffuse;
  s_diffuse = diffuse;
  return ret;
}


void pointlike::SourceLikelihood::addBackgroundPointSource(const pointlike::SourceLikelihood* fit)
{
  if( fit==this){
    throw std::invalid_argument("SourceLikelihood::setBackgroundFit: cannot add self as background");
  }
  if( m_background==0){
    throw std::invalid_argument("SourceLikelihood::setBackgroundFit: no diffuse background");
  }
  
  if( s_verbose>0 ) {
    std::cout << "Adding source " << fit->name() << " to background" << std::endl;
  }
  m_background->add(fit);
  setDir(dir()); // recomputes background for each SimpleLikelihood object
}

void pointlike::SourceLikelihood::clearBackgroundPointSource()
{
  if( m_background==0){
    throw std::invalid_argument("SourceLikelihood::setBackgroundFit: no diffuse to add to");
  }
  while( m_background->size()>1) {
    m_background->pop_back();
  }
  setDir(dir()); // recomputes background for each SimpleLikelihood object
}

const skymaps::SkySpectrum * pointlike::SourceLikelihood::background()const
{
    return m_background;
}

/// @brief set radius for individual fits
void pointlike::SourceLikelihood::setDefaultUmax(double umax)
{ 
  pointlike::ExtendedLikelihood::setDefaultUmax(umax); 
}


// double SourceLikelihood::set_gamma_level(int level, double v)
// {
//     double t = s_gamma_level[level]; s_gamma_level[level]=v; 
//     return t;
// }

// double SourceLikelihood::set_sigma_level(int level, double v)
// {
//     double t = s_sigma_level[level]; s_sigma_level[level]=v; 
//     return t;
// }



// /// @brief get the starting, ending levels used
// int SourceLikelihood::minlevel(){return s_minlevel;}

// int SourceLikelihood::maxlevel(){return s_maxlevel;}

// void SourceLikelihood::set_levels(int minlevel, int maxlevel)
// { s_minlevel= minlevel; s_maxlevel = maxlevel;
// }

double pointlike::SourceLikelihood::set_tolerance(double tol)
{
  double old(pointlike::ExtendedLikelihood::tolerance());
  pointlike::ExtendedLikelihood::setTolerance(tol);
  return old;
}

// double SourceLikelihood::gamma_level(int i)
// {
//     return s_gamma_level.at(i);
// }
// double SourceLikelihood::sigma_level(int i)
// {
//     return s_sigma_level.at(i);
// }

//=======================================================================
//         SLdisplay implementation
pointlike::SLdisplay::SLdisplay(const pointlike::SourceLikelihood & psl, int mode)
  : m_psl(psl)
    , m_mode(mode)
{}

double SLdisplay::value(const astro::SkyDir& dir, double e)const{
  return m_psl.display(dir, e, m_mode);
}

///@brief integral for the energy limits, in the given direction -- not impleme
double SLdisplay::integral(const astro::SkyDir& dir, double a, double b)const{
  return value(dir, sqrt(a*b));
}

std::string SLdisplay::name()const
{
  static std::string type[]={"density", "data", "background", "fit", "residual"};
    if( m_mode<0 || m_mode>4) return "illegal";
    return m_psl.name()+"--"+type[m_mode];
    
}
