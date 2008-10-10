/** @file PhotonBinner.cxx
@brief implement class BinnedPhotonData 

$Header: /nfs/slac/g/glast/ground/cvs/skymaps/src/PhotonBinner.cxx,v 1.12 2008/08/06 19:42:24 mar0 Exp $
*/

#include "pointlike/FlexibleBinner.h"
#include "skymaps/Band.h"
#include "skymaps/IParams.h"
#include "astro/Photon.h"
#include <algorithm>
#include <functional>
#include <map>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <iterator>

typedef std::vector<double>::iterator vdIter;
typedef std::vector<int>::iterator viIter;
 

namespace {
// Bin scheme id: pass6/spectrum:0
  const double bins_pass6_s0[13]={1.,50.,100.,200.,400.,800.,1600.,3200.,6400.,12800.,25600.,51200.,1e6};
  const double gf_pass6_s0[12]={2., 1.885, 1.933, 2.252 , 2.198 , 2.375 , 2.572 , 2.782 , 2.664 , 2.311 , 1.947 , 1.841 };
  const double sf_pass6_s0[12]={4.,1.76 ,1.14 ,0.728 , 0.4096 , 0.2399 , 0.1477 , 0.0949 , 0.0593 , 0.0344 , 0.0152 , 0.015 };
  const double gb_pass6_s0[12]={2., 1.1, 1.213, 2.201 , 2.120 , 2.188 , 2.258 , 2.418 , 2.427 , 2.173 , 2.002 , 1.842 };
  const double sb_pass6_s0[12]={4., 3.57, 1.84, 1.264 , 0.7116 , 0.4096 , 0.2446 , 0.1541 , 0.0988 , 0.0620 , 0.0421 , 0.0291 }; 

// Bin scheme id: pass6/spectrum:-
  const double bins_pass6_sm[9]={1.,25.,100.,400.,1600.,6400.,25600.,102400.,1e6};
  const double gf_pass6_sm[8]={2.,1.598,1.522,1.881,2.362,2.466,1.888,1.837};
  const double sf_pass6_sm[8]={4.,1.77e+00,7.82e-01,2.98e-01,1.19e-01,4.89e-02,1.10e-02,5.e-03}; 
  const double gb_pass6_sm[8]={2.,1.1,1.335,1.629,2.033,2.286,1.951,1.852};
  const double sb_pass6_sm[8]={4.,3.55e+00,1.35e+00,4.90e-01,1.93e-01,8.36e-02,3.83e-02,2.62e-02};

// Old binning scheme: pass6/classic
  const double bins_pass6_old_front[11] = {10, 42.555, 100, 235,   552.26, 1297.79 , 3050,  7167.2, 16842.7, 39580., 1e6};
  const double bins_pass6_old_back [9]  = {10, 79.14,  186, 437.1, 1027.19, 2413.9,  6500., 21000, 1e6};
  const double gf_pass6_old[10] ={2., 1.785, 1.874, 1.981, 2.003, 2.225, 2.619, 2.527, 2.134, 1.923};
  const double sf_pass6_old[10] ={4., 1.77e+00, 1.06e+00, 5.69e-01, 2.85e-01, 1.54e-01, 9.16e-02, 5.15e-02, 2.52e-02, 1.60e-02};
  const double gb_pass6_old[8]  ={2., 1.737, 1.742, 1.911, 2.008, 2.009, 2.207, 1.939};
  const double sb_pass6_old[8]  ={4., 2.18e+00, 1.17e+00, 6.00e-01, 3.09e-01, 1.61e-01, 8.43e-02, 3.90e-02 };
  const int    level_pass6_old[10]  ={4,5,6,7,8,9,10,11,12,13};
};

namespace pointlike {

FlexibleBinner::FlexibleBinner(const std::string& id, const int pixel_density):classic_mode(false){

   if(id=="p6_v1/classic"){
      m_bins              =std::vector<double>(bins_pass6_old_front,bins_pass6_old_front+11);
      m_bins_classic_back =std::vector<double>(bins_pass6_old_back,bins_pass6_old_back+9);
      classic_mode=true;
      m_gammaFront= std::vector<double>(gf_pass6_old,gf_pass6_old+10);
      m_sigmaFront= std::vector<double>(sf_pass6_old,sf_pass6_old+10);
      m_gammaBack = std::vector<double>(gb_pass6_old,gb_pass6_old+8);
      m_sigmaBack = std::vector<double>(sb_pass6_old,sb_pass6_old+8);
      m_level     = std::vector<int>(level_pass6_old,level_pass6_old+10);
      return;
   };
   if(id=="p6_v1/spectrum:0"){
      m_bins	  =std::vector<double>(bins_pass6_s0,bins_pass6_s0+13);
      m_gammaFront=std::vector<double>(gf_pass6_s0,gf_pass6_s0+12);
      m_sigmaFront=std::vector<double>(sf_pass6_s0,sf_pass6_s0+12);
      m_gammaBack =std::vector<double>(gb_pass6_s0,gb_pass6_s0+12);
      m_sigmaBack =std::vector<double>(sb_pass6_s0,sb_pass6_s0+12);
      calc_healpix_level(pixel_density);
      return;
   };
   if(id=="p6_v1/spectrum:-"){
      m_bins      =std::vector<double>(bins_pass6_sm,bins_pass6_sm+9);
      m_gammaFront=std::vector<double>(gf_pass6_sm,gf_pass6_sm+8);
      m_sigmaFront=std::vector<double>(sf_pass6_sm,sf_pass6_sm+8);
      m_gammaBack =std::vector<double>(gb_pass6_sm,gb_pass6_sm+8);
      m_sigmaBack =std::vector<double>(sb_pass6_sm,sb_pass6_sm+8);
      calc_healpix_level(pixel_density);
      return;
   };
   throw std::runtime_error("Did not find a binning with the id specified.");
};


skymaps::Band FlexibleBinner::operator()(const astro::Photon& p)const
{
    double energy ( p.energy() );

    if( energy>=m_bins.back()) energy=0.9999 * m_bins.back();
    if( energy<=m_bins.front()) energy=1.0001 * m_bins.front();

    int event_class (  p.eventClass() );
    if( event_class<0 || event_class>15) event_class=0; // should not happen?

    const std::vector<double>& m_gamma ( event_class<=0? m_gammaFront : m_gammaBack );
    const std::vector<double>& m_sigma ( event_class<=0? m_sigmaFront : m_sigmaBack );

    const std::vector<double>& m_binsf ( (classic_mode && event_class>0) ? m_bins_classic_back : m_bins);

    std::vector<double>::const_iterator bin_it=
        std::lower_bound(m_binsf.begin(), m_binsf.end(), energy, std::less<double>());
    int i =  bin_it-m_binsf.begin()-1;
    double sigma = m_sigma[i]*M_PI/180.;
    double gamma = m_gamma[i];
    unsigned int nside = 1<<m_level[i];

    double elow  = m_binsf[i];
    double ehigh = m_binsf[i+1];

//    std::cout<<"Flexible binner: e="<<energy<<" eclass="<<event_class<<" elow="<<elow<<" ehigh="<<ehigh
//             <<" sigma="<<sigma<<" gamma="<<gamma<<" nside="<<nside<<" level="<<m_level[i]<<std::endl;
    return skymaps::Band(nside, event_class, elow, ehigh, sigma, gamma);
}

    void FlexibleBinner::calc_healpix_level(const int density){
       m_level.resize(m_sigmaFront.size());
       vdIter si=m_sigmaFront.begin();       
       viIter li=m_level.begin();     
       for (;si!=m_sigmaFront.end();si++,li++){
          double nside = 2*180./(3*(*si));
	  *li = int(log(nside)/log(2))+1+density; 
          *li=std::min(13,*li);
	  *li=std::max(1,*li);
       };	  	   
    };

};
