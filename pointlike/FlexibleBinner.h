/** @file FlexibleBinner.h
@brief declare class FlexibleBinner

$Header: /nfs/slac/g/glast/ground/cvs/skymaps/skymaps/FlexibleBinner.h,v 1.6 2008/05/02 23:29:13 burnett Exp $

*/
#ifndef pointlike_FlexibleBinner_h
#define pointlike_FlexibleBinner_h

//#include "skymaps/BinnedPhoton.h"

namespace astro {class Photon;}
//namespace skymaps {class BinnedPhoton;}

#include <string>
#include <vector>


#include "skymaps/Band.h"
#include "skymaps/PhotonBinner.h"


namespace pointlike {
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    /** @class FlexibleBinner
    @brief manage binning of photons 
    */
    class FlexibleBinner: public skymaps::PhotonBinner { 
    public:
         /**@brief ctor (default) 
        */
        FlexibleBinner(const std::string& id,const int pixel_density=0);

        /** @brief ctor  takes arguments of the left edge bin energy in MeV
        */
        FlexibleBinner(const std::vector<double>& bins,
		       const std::vector<double>& gamma_front,
		       const std::vector<double>& sigma_front,
		       const std::vector<double>& gamma_back,
		       const std::vector<double>& sigma_back,
		       const std::vector<int>& level):
	   classic_mode(false),
	   m_bins(bins),
	   m_bins_classic_back(bins),
	   m_gammaFront(gamma_front),
	   m_sigmaFront(sigma_front),
	   m_gammaBack(gamma_back),
	   m_sigmaBack(sigma_back),
	   m_level(level)
	{};

        ~FlexibleBinner(){};
	
        ///@brief bin a photon by returning an appropriate Band object
        skymaps::Band operator()(const astro::Photon& photon)const;

    private:
        /**@brief setupbins  sets up bin to pixel connection with current bin set
        */
        void calc_healpix_level(const int density);
	
	bool classic_mode;
	
        std::vector<double> m_bins;               //the energy of each left bin edge
        std::vector<double> m_bins_classic_back;  //the energy of each left bin edge (classic front/back)
        std::vector<double> m_gammaFront;         //the PSF gamma for front bins
        std::vector<double> m_sigmaFront;         //the PSF sigma for front bins
        std::vector<double> m_gammaBack;          //the PSF gamma for back bins
        std::vector<double> m_sigmaBack;          //the PSF sigma for back bins
        std::vector<int>    m_level;              //the mapping between energy bins and healpix levels
    };
}

#endif
