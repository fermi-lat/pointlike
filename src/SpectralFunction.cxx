/** @file SpectralFunction.cxx
@brief implementation of class SpectralFunction

$Header$
*/
#include "pointlike/SpectralFunction.h"

#include "skymaps/Band.h"
#include "skymaps/EffectiveArea.h"
#include "skymaps/Exposure.h"

#include <cmath>
#include <stdexcept>

using namespace pointlike;

namespace {
    std::string names[] = {"PowerLaw", "ExpCutoff", "BrokenPowerLaw"};
}
int SpectralFunction::s_n(8);
double SpectralFunction::s_e0(1000.);
double SpectralFunction::s_flux_scale(1.);
std::vector<const skymaps::Exposure*> SpectralFunction::s_exposures;

void SpectralFunction::set_simpson(int n){
    s_n=n;
    if( s_n<2 || (s_n&2) !=0 || s_n>20){
        throw std::invalid_argument("SpectralFunction::set_simpson--bad Simpsons rule count: n must be even, <20");
    }
}


void SpectralFunction::set_exposures(const skymaps::Exposure* front, const skymaps::Exposure* back)
{
    s_exposures.clear();
    s_exposures.push_back(front);
    s_exposures.push_back(back);
}
const skymaps::Exposure* SpectralFunction::get_exposure(int n){
    if( s_exposures.size() <2 ) {
        throw std::invalid_argument("SpectralFunction: exposures not set");
    }
    return s_exposures.at(n);
}

SpectralFunction::SpectralFunction(std::string name, const std::vector<double>& pars)
: m_name(name)
, m_type(-1)
, m_pars(pars)
{
    if( s_exposures.size() <2 ) {
        throw std::invalid_argument("SpectralFunction: exposures not set");
    }
    for(int i=0; i< sizeof(names)/sizeof(std::string); ++i){
        if( name== names[i]) {
            m_type=i;
            return;
        }
    }
    std::cerr << "unrecognized name: valid names are: ";
    for(int i=0; i< sizeof(names)/sizeof(std::string); ++i){ std::cerr << " "<<names[i];}
    std::cerr << std::endl;
    throw std::invalid_argument("SpectralFunction: name " + name +" is not recognized ");
}

double SpectralFunction::value(double e)const
{
    switch( m_type ){
        case PowerLaw:
/*          n0         log10 differential flux at e0 MeV
            gamma      (absolute value of) spectral index
            return (10**n0/self.flux_scale)*(self.e0/e)**gamma
  */
            {
            double n0(m_pars[0]), gamma(m_pars[1]);
            return std::pow(10.,n0)/s_flux_scale*std::pow(s_e0/e, gamma); 
            }
        case ExpCutoff:
/*              n0         differential flux at e0 MeV
  gamma      (absolute value of) spectral index
  cutoff     e-folding cutoff energy (MeV)
      """
   def __call__(self,e):
      n0,gamma,cutoff=self.p
      if cutoff < 0: return 0
      return (10**n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-e/cutoff)
      */
            {
            double n0(m_pars[0]), gamma(m_pars[1]), cutoff(m_pars[2]);
            return std::pow(10.,n0)/s_flux_scale*std::pow(s_e0/e, gamma)*exp(-e/cutoff);
            }
        default:
            return 0; // temporary
    }
}

// self.psf_correction   = 1 - (1+sl.umax()/sl.gamma())**(1-sl.gamma())

double SpectralFunction::expected(const astro::SkyDir& dir,const skymaps::Band& band)const
{
    double a(band.emin()), b(band.emax()); // range of integration
    const skymaps::Exposure& exposure(* get_exposure(band.event_class()) ); // exposure object to use

    double step( log(b/a)/s_n ), // step in log scale
           ratio( exp(step) ), // ratio of energies
           c(2.),              // initial simpsons
           e(a);              //  inital energy

    double result( a* exposure(dir,a)*value(a) ); // value at low end
    for( int i = 1; i< s_n; ++i ){
        e *= ratio; // next energy
        c = 6-c;   // toggle Simpsons coefficient
        result += c * e* exposure(dir,e)*value(  e); 
    }
    
    result += b*exposure(dir,b)*value(b);// value at high end.
    return result*step/3.;
}

