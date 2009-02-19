/** @file Draw.cxx
@brief implementation of Draw

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Draw.cxx,v 1.21 2008/12/05 20:20:40 funk Exp $

*/


#include <iomanip>
#include "pointlike/Draw.h"
#include "pointlike/Data.h"

#include "pointlike/PointSourceLikelihood.h"
#ifndef WIN32
#include "pointlike/SourceLikelihood.h"
#endif

#include "astro/SkyDir.h"

#include "skymaps/SkyImage.h"
#include "skymaps/BinnedPhotonData.h"
#include "skymaps/Band.h"
#include "skymaps/SkySpectrum.h"

using namespace pointlike;
using astro::SkyDir;
using skymaps::BinnedPhotonData;
using skymaps::SkyImage;

namespace pointlike {
  struct lt{
    bool operator()(astro::SkyDir d1, astro::SkyDir d2) const {
      if (fabs(d1.b() - d2.b()) > 1e-5) return (d1.b() < d2.b());
      else return (d1.l() < d2.l());
    }
  };

  typedef std::map<astro::SkyDir, std::vector<double>, pointlike::lt > TSCache;


  /// @class SkyTS
  /// @brief adapt a Binned Photon Data to give TS for bins with given energy
  template <class Likelihood>
  class SkyTS : public astro::SkyFunction {

  private: 
    skymaps::BinnedPhotonData& m_data;
    const skymaps::SkySpectrum& m_background;
    unsigned int m_nenergybins;
    unsigned int m_bin;
    double m_emin;
    double m_minalpha;
    std::vector<TSCache>& m_Cache;  // vector of TSCaches
    // [0] TS-values
    // [1] alphas
    // [2] nCounts
    // [3] exposure
    // [4] model
    int m_returnValue;
    // see cache Index
    Likelihood * m_ps;

  public: 
    SkyTS(BinnedPhotonData& data, Likelihood* pslike,
	  std::vector<TSCache>& cacheVec, int energyBin, int retVal, double emin=100., double minalpha=0.):
      m_data(data),
      m_nenergybins(0),
      m_bin(energyBin),
      m_ps(pslike),
      m_emin(emin),
      m_minalpha(minalpha),
      m_Cache(cacheVec),
      m_returnValue(retVal){

      for( skymaps::BinnedPhotonData::const_iterator bit = m_data.begin(); 
	   bit!=m_data.end(); ++bit){
	if (bit->emax() < m_emin) continue;
	m_nenergybins++;
//	std::cout << "SkyTS> Added energyBin: " << m_nenergybins << " " 
//		  << bit->emin() << " " << bit->emax() << std::endl;
      }
       m_ps = pslike;
    }


    ~SkyTS() {};  

    double operator() (const astro::SkyDir& sd) const {

	if (m_Cache[0].find(sd) == m_Cache[0].end()){

	  m_ps->setDir(sd,true);
	  // 	ps->setup(m_data);
	  m_ps->maximize();

	  std::vector<double> ts(m_nenergybins+1,0);
	  std::vector<double> alphas(m_nenergybins+1,0);
	  std::vector<double> alphaErrors(m_nenergybins+1,0);
	  std::vector<double> counts(m_nenergybins+1,0);
	  std::vector<double> exposure(m_nenergybins+1,0);

	  double t = m_ps->TS();
	  ts[0] = t;

	  double totalCounts   = 0;
	  for (unsigned int i = 0; i < m_nenergybins; ++i){
	    double e_ts = m_ps->at(i)->TS();
	    double e_alpha = m_ps->at(i)->alpha();
	    double e_alphaerror = m_ps->at(i)->sigma_alpha();
	    double e_counts = m_ps->at(i)->photons();
	    double e_exposure =m_ps->at(i)->exposure();
	    
	    if (e_alpha < m_minalpha) {
	      e_ts = 0;
	      e_alpha = 0;
	      e_alphaerror = 0;
	    }
	    ts[i+1]=e_ts;
	    alphas[i+1]=e_alpha;
	    alphaErrors[i+1]=e_alphaerror;
	    counts[i+1]=e_counts;
	    exposure[i+1]=e_exposure;
	    totalCounts += e_counts;
// 	      if (e_ts > 0)
// 	    std::cout << "Using: " << std::setprecision(2) << e_alpha 
// 		      << " "  << e_ts << " " <<  std::setprecision(0);

//	    std::cout << "    Bin: " << i << " ts: " << e_ts << " alpha: " << e_alpha
//		      << " +- " << e_alphaerror << " cnts: " << e_counts << " exposure: "
//		      << e_exposure << " " << totalCounts << std::endl;
	  }

	  counts[0] = totalCounts;
	  m_Cache[0][sd] = ts;
	  m_Cache[1][sd] = alphas;
	  m_Cache[2][sd] = alphaErrors;
	  m_Cache[3][sd] = counts;
	  m_Cache[4][sd] = exposure;

   	  std::cout<<".";std::cout.flush();
	  
//	    m_ps->printSpectrum();
//	  std::cout << "l="<<sd.l()<<" b="<<sd.b()<<" TS=" << t << " Counts=" 
//		    << totalCounts <<" exposure=" << exposure[m_nenergybins-1] << std::endl;
// 	    std::cout << "-----" << std::endl;
	  return m_Cache[m_returnValue][sd][m_bin];
	    //	  return ts[m_bin];
	} else {
	  return (*m_Cache[m_returnValue].find(sd)).second[m_bin];
	  //	  return (*m_cache[0].find(sd)).second[m_bin];
	};
      };
   };
//--------------------------------------------------------------------------------------


}

Draw::Draw(BinnedPhotonData& map, const skymaps::SkySpectrum* background, 
	   bool ts, double eMin, double minalpha, bool sourcelike)
: m_map(map)
, m_background(background)
, m_galactic(true)
, m_zenith(false)
, m_proj("")
, m_exposure(0) // default: do not apply
  , m_layers(1)
  , m_emin(eMin)
  , m_minalpha(minalpha)
  , m_ts(ts)
  , m_sourcelike(sourcelike)
{ pointlike::PointSourceLikelihood::set_energy_range(m_emin, 1e6); }

Draw::Draw(Data& data)
: m_map(data.map())
, m_background(0)
, m_galactic(true)
, m_zenith(false)
, m_proj("")
, m_exposure(0) // default: do not apply
, m_layers(1)
  , m_ts(false)
{}

void Draw::density(const astro::SkyDir& dir, std::string outputFile, double pixel, 
                      double fov, bool smooth, int mincount) {

    std::string proj (fov>90? "AIT":"ZEA");
    if( !m_proj.empty()){ proj = m_proj;}
                
    SkyImage image(dir, outputFile, pixel, fov, m_layers, proj,  m_galactic, m_zenith);
    image.fill(SkyDensity(m_map, smooth, mincount, m_exposure), 0); // PhotonMap is a SkyFunction of the density 

    std::cout 
        <<   "\t minimum "<< image.minimum()
        << "\n\t maximum "<< image.maximum()
        << "\n\t average "<< (image.total()/image.count())
        << std::endl;
		   
};		   
		   
template <class Likelihood>		   
void Draw::mapTS(Likelihood* like, std::string outputFile, double pixel, 
                 double fov, bool smooth, int mincount) {
#ifndef WIN32
    std::string proj (fov>90? "AIT":"ZEA");
    if( !m_proj.empty()){ proj = m_proj;}

    std::vector<double> eMin;
    eMin.push_back(1.); // sum of all energies 
    for( skymaps::BinnedPhotonData::const_iterator bit = m_map.begin(); 
	 bit!=m_map.end(); ++bit){
      if (bit->emax() < m_emin) continue;
//      std::cout << "SKYTS>  energy: " << bit->emin() << " " << bit->emax() 
//		<< std::endl;
      eMin.push_back(bit->emin());
    }

    std::cout << "SkyTS> Constructing SkyImage with: " << eMin.size() << " layers" 
	      << std::endl;

    TSCache tsCache;
    TSCache alphaCache;
    TSCache countsCache;
    TSCache exposureCache;
    TSCache modelCache;
    
    std::vector<TSCache> caches;
    caches.push_back(tsCache);
    caches.push_back(alphaCache);
    caches.push_back(countsCache);
    caches.push_back(exposureCache);
    caches.push_back(modelCache);
    
    astro::SkyDir image_center=like->dir();
    
    for (int returnValue = 0; returnValue < caches.size(); ++returnValue){
      std::string outname = outputFile;
      std::string replacestring("TS.fits");
      if (outname.find("TS.fits") == std::string::npos) replacestring = std::string(".fits");
	
      if (returnValue == 1) outname.replace(outname.rfind(replacestring), replacestring.size(), "alphas.fits");
      if (returnValue == 2) outname.replace(outname.rfind(replacestring), replacestring.size(), "alphasErrors.fits");
      else if (returnValue == 3) outname.replace(outname.rfind(replacestring), replacestring.size(), "counts.fits");
      else if (returnValue == 4) outname.replace(outname.rfind(replacestring), replacestring.size(), "exposure.fits");

      std::cout << "Generating map " << outname << " centered at l="<<image_center.l()<<" b="<<image_center.b()<<std::endl;
      SkyImage image2(image_center, outname, pixel, fov, eMin.size(), proj,  m_galactic, m_zenith); 
      image2.setEnergies(eMin);

      for (unsigned int i = 0; i <eMin.size(); ++i){
        std::cout << " SkyTS> Filling layer: " << i+1 << " emin: " << m_emin;
	image2.fill(SkyTS<Likelihood>(m_map, like, caches, i, returnValue, m_emin,m_minalpha), i);
	std::cout<<"o"<<std::endl;  
      }
      like->setDir(image_center,false);
    }	 
#endif
};		 

void Draw::TS(SourceLikelihood* like, std::string outputFile, double pixel, 
                 double fov, bool smooth, int mincount){
    mapTS<SourceLikelihood>(like,outputFile,pixel,fov,smooth,mincount);		 
};		 

//void Draw::TS(PointSourceLikelihood* like, std::string outputFile, double pixel, 
//                 double fov, bool smooth, int mincount){
//    mapTS<PointSourceLikelihood>(like,outputFile,pixel,fov,smooth,mincount);		 
//};

void Draw::region(const astro::SkyDir& dir, std::string outputFile, double pixel, 
                  double fov, bool smooth, int mincount)
{
    density(dir,outputFile,pixel,fov,smooth,mincount);
    
/*
    if (m_ts && m_background){
      std::string outputTS = outputFile;
      if (outputTS.find(".fits") != std::string::npos)
	outputTS.replace(outputTS.find(".fits"), 0, "_TS");
      else 
	outputTS = std::string("tsvalues.fits");
      TS(dir,outputTS,pixel,fov,smooth,mincount); 
    }
*/

//--------------------------------------------------------------------------------------

    /// @class SkyCount
    /// @brief adapt a BinnedPhotonData to give counts for bins with given energy
    class SkyCount : public astro::SkyFunction {
    public:
        SkyCount(const BinnedPhotonData& data, double energy): 
          m_data(data), m_energy(energy){}

          double operator()(const astro::SkyDir & sd) const {
              double  value = m_data.value(sd, m_energy); 
              return value;    
          }
    private:
        const BinnedPhotonData& m_data;
        double m_energy;
    };
#if 0

    std::cout << "Filling layers "<< layer << " and above with ... ";
    int points(0);
    for (BinnedPhotonData::const_iterator it = m_map.begin();
        it!= m_map.end(); ++it)
    {
        if (it->first.level() >= level) {
            int layer = it->first.level() - m_map.minLevel() + 1;
            if(image.addPoint((it->first)(), it->second, layer)) points+= it->second;
        }
    }
    std::cout <<  points << " hit display pixels" << std::endl;

    for(int layer=1; layer!=m_layers; ++m_layers){
       image.fill(SkyCount(m_map), 100.);
    }
#endif

    std::cout << "Writing image to file \""
        << outputFile << "\""<<std::endl;

}
void Draw::sky(std::string outputfile, double pixel, bool smooth, int mincount)
{
    region(SkyDir(0,0, SkyDir::GALACTIC), outputfile, pixel, 180., smooth, mincount);
}

void Draw::googleSky(std::string outfile, double pixelsize, bool smooth, int mincount)
{
    std::string proj = m_proj;    m_proj = "CAR";
    bool gal = m_galactic; m_galactic = false;
    region(SkyDir(180,0), outfile, pixelsize, 180, smooth, mincount);
    m_proj = proj;
    m_galactic = gal;
}


//-------------------------------------------------------------------
//             SkyDensity methods
//-------------------------------------------------------------------
SkyDensity::SkyDensity(const skymaps::BinnedPhotonData& data, bool smooth, int mincount
                       ,const astro::SkyFunction* exposure)
: m_data(&data)
, m_band(0)
, m_smooth(smooth)
, m_mincount(mincount)
, m_exposure(exposure)
{}

SkyDensity::SkyDensity(const skymaps::Band& band)
: m_data(0)
, m_band(&band)
, m_smooth(false)
, m_mincount(0)
, m_exposure(0)
{}

double SkyDensity::operator()(const astro::SkyDir & sd) const 
{
    double  value;

    if( m_band!=0) {
        // only selected band
        return (*m_band)(sd);
    }

    if (m_smooth)
        value = m_data->smoothDensity(sd, m_mincount);
    else
        value = m_data->density(sd);

    if(m_exposure!=0){
        // note we are not using energy dependence here
        double exposure( (*m_exposure)(sd) );
        if( exposure>0.) value /= exposure;
    }
    return value;    
}
