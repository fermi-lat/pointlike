/** @file HealpixDiffuseFunc.cxx
    @brief implement HealpixDiffuseFunc
$Header: /nfs/slac/g/glast/ground/cvs/skymaps/src/HealpixDiffuseFunc.cxx,v 1.2 2008/03/10 16:38:43 burnett Exp $
*/

#include "pointlike/HealpixDiffuseFunc.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"
#include "tip/Table.h"
#include "tip/IFileSvc.h"
#include <cmath>
#include <map>
#include <stdexcept>

using namespace skymaps;


HealpixDiffuseFunc::HealpixDiffuseFunc(std::string diffuse_healpix_file, double energy, bool interpolate)
  : SkySpectrum(energy)
  , m_name(diffuse_healpix_file)
  , m_fitio(diffuse_healpix_file)
  , m_skymap(m_fitio.skymap())
  , m_energies(m_fitio.energies()) {

    m_emin = m_energies.front();
    m_emax = m_energies.back();
    std::cout << "HealpixDiffuseFunc: read file "<< diffuse_healpix_file <<", with " 
        << layers() << " energies from " << m_emin << " to " << m_emax << std::endl;
    setEnergy(energy);
    
}

HealpixDiffuseFunc::~HealpixDiffuseFunc()
{
}

double HealpixDiffuseFunc::isotropicFlux(double energy)const
{
    static double E0(100), flux(1.5E-5), alpha(1.1);
    return flux*alpha*pow(E0/energy, alpha)/energy;
}

size_t HealpixDiffuseFunc::layer(double e)const
{
    if( e< m_emin  ){
        std::stringstream error; 
        error << "Diffuse function: energy out of range: "<< e ;
        throw std::invalid_argument(error.str());
    }
    if( e>m_emax ) {
        return m_energies.size();
    }
    size_t step(0);
    for( std::vector<double>::const_iterator it( m_energies.begin()); it!=m_energies.end(); ++it, ++step){
        if( (*it) > e ) break;
    }
    return step-1;
}


double HealpixDiffuseFunc::value(const astro::SkyDir& dir, double e)const
{
///never use the operator [] for accecssing the healpix_array with an index
///it is overloaded by astro::SkyDir and will pass on your index without warning or error
///as a parameter to astro::SkyDir. Always use at()

    size_t idx = m_skymap.healpix().pixel(dir).index();
    size_t l(layer(e)); 
    if( l==m_energies.size()) return (m_skymap.at(idx))[l-1];

    double e1(m_energies[l]), e2(m_energies[l+1]);
    double f1( (m_skymap.at(idx))[l] ), f2( (m_skymap.at(idx))[l+1] );
    double alpha ( log(f1/f2)/log(e2/e1) );

//    std::cout<<"HealpixDiffuseFunc::value: l="<<dir.l()<<" b="<<dir.b()<<" layer="<<l<<" idx="<<idx<<" f1="<<f1<<" f2="<<f2
//             <<" e="<<e<<" e1="<<e1<<" e2="<<e2<<" alpha="<<alpha<<" result="<<(f1*pow( e1/e, alpha))
//	     <<std::endl;
    
    return f1*pow( e1/e, alpha);
};

double HealpixDiffuseFunc::integral(const astro::SkyDir& dir, double a, double b)const
{
    // estimate integral by assuming power-law from start to end
    double fa(value(dir, a)), fb(value(dir,b));
    double q ( 1. - log(fa/fb)/log(b/a) );
//    std::cout<<"HealpixDiffuseFunc::integral: l="<<dir.l()<<" b="<<dir.b()<<" fa="<<fa<<" fb="<<fb
//    <<" emin="<<a<<" emax="<<b<<" result="<<(fa* a * (pow(b/a, q)-1)/q)<<std::endl;
    return fa* a * (pow(b/a, q)-1)/q;
}

std::vector<double> HealpixDiffuseFunc::integral(const astro::SkyDir& dir, const std::vector<double>&energies)const
{
    std::vector<double> result;
    static double infinity (300000);

    for( std::vector<double>::const_iterator it = energies.begin(); it!=energies.end(); ++it){
        std::vector<double>::const_iterator next(it+1);
        double a = *it;
        double b = next!=energies.end()? *next : infinity;
        result.push_back(integral(dir,a,b));
    }
    return result;
}



/**\brief Read skymap from a fits file*/

HealpixDiffuseFunc::FitsIO::FitsIO(const std::string& filename)
  : fileSvc(tip::IFileSvc::instance())
  , m_skymap(0)
  , m_energies(0)
  {
     read(filename);
  };

HealpixDiffuseFunc::FitsIO::~FitsIO(){
  if (m_skymap) delete m_skymap;
};


void HealpixDiffuseFunc::FitsIO::read(const std::string &filename){
     int nRows, nEnergies, nSide;
     std::string ordering;

     tip::Table& eTab = *(fileSvc.editTable(filename,"ENERGIES"));
     const std::vector<std::string>& fields=eTab.getValidFields();
    
     tip::Table::Iterator itor = eTab.begin();
     for (; itor != eTab.end(); ++itor) m_energies.push_back((*itor)[fields.front()].get());

     tip::Table& skymapTab = *(fileSvc.editTable(filename,"SKYMAP"));
     tip::Header& skymapHeader = skymapTab.getHeader();

     skymapHeader["NAXIS2"].get<int>(nRows);
     skymapHeader["NBRBINS"].get<int>(nEnergies);
     skymapHeader["NSIDE"].get<int>(nSide);
     skymapHeader["ORDERING"].get<std::string>(ordering);

     healpix::Healpix::Ordering Ordering;
     if (ordering=="NEST") Ordering=healpix::Healpix::NEST;
     else if  (ordering=="RING") Ordering=healpix::Healpix::RING;
     else throw std::runtime_error("Unknown ordering scheme encountered in Healpix FitsIO.");

     healpix::Healpix hpix(nSide,Ordering,astro::SkyDir::GALACTIC);
     m_skymap=new healpix::HealpixArray<std::vector<double> > (hpix);
     healpix::HealpixArray<std::vector<double> >::iterator hpx_itor;
     double * bins=new double[nEnergies];

     hpx_itor=m_skymap->begin();
     
     itor = skymapTab.begin(); hpx_itor=m_skymap->begin(); size_t i=0;
     for (; itor != skymapTab.end() && hpx_itor!=m_skymap->end(); ++itor,++hpx_itor,++i) {
         (*itor)["spectra"].get<double>(0,nEnergies,bins);
	 (*hpx_itor).resize(nEnergies);
	 (*hpx_itor)=std::vector<double>(bins,bins+nEnergies);
     };
     	 
/**/	 
     delete bins;
     std::cout<<"Read "<<i<<" pixels."<<std::endl;
     
};

