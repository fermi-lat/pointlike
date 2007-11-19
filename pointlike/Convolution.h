/** @file Convolution.h
@brief Convolves healpix maps with a skyfunction

@author M. Roth 

$Header: /nfs/slac/g/glast/ground/cvs/healpix/healpix/Map.h,v 1.4 2007/06/04 22:14:25 mar0 Exp $
*/

#include "astro/SkyFunction.h"
#include "map_tools/PhotonMap.h"


namespace healpix {

    class Convolution {

    public:
        /*
        * Convolution of a sky map with arbitrary sky function
        * @param sf Map as a sky function
        * @param ker convolution function
        * @param level resolution of healpix map, level = log2(nside)
        */
        Convolution(const astro::SkyFunction &sf, const astro::SkyFunction &ker, int level);

        /*
        * Convolution of a sky map with the LAT psf a particular energy band (assumes azithmuthal symmetry of psf)
        */
        Convolution(const astro::SkyFunction &sf, int level);
        
        /// returns photonmap
        map_tools::PhotonMap&  map() {return m_map;}
    
    private:

        map_tools::PhotonMap m_map;

    };


}