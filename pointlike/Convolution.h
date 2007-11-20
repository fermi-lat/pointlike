/** @file Convolution.h
@brief Convolves healpix maps with a skyfunction

@author M. Roth 

$Header: /nfs/slac/g/glast/ground/cvs/healpix/healpix/Map.h,v 1.4 2007/06/04 22:14:25 mar0 Exp $
*/

#include "pointlike/SkySpectrum.h"
#include "healpix/Map.h"
#include <map>

namespace pointlike {

    class Convolution : public pointlike::SkySpectrum, public std::map<int,healpix::Map<double> >{

       
    public:
        /*
        * Convolution of a sky map with arbitrary sky function
        * @param sf Map as a sky function
        * @param ker convolution function
        * @param level resolution of healpix map, level = log2(nside)
        */
        Convolution(const pointlike::SkySpectrum &sf, const pointlike::SkySpectrum &ker, double energy, int level);

        /*
        * Convolution of a sky map with the LAT psf a particular energy band (assumes azithmuthal symmetry of psf)
        * @param sf SkyFunction
        * @param energy determines shape of PSF (from IRF)
        * @param level max resolution (healpix bins) 6:fast-coarse -> 11:slow-fine
        */
        Convolution(const pointlike::SkySpectrum& ss, double energy, int level);

        void createConv(const pointlike::SkySpectrum&sf, const pointlike::SkySpectrum& ker, double energy);
        
        void createConv(const pointlike::SkySpectrum&sf, double energy);

        int layer(double e, bool front) const;

        double getvalue(const astro::SkyDir& dir, double e);

        /// returns photonmap
        //map_tools::PhotonMap&  map(int energy) {return m_map[energy];}

        ///@brief interpolate table 
        ///@param e energy in MeV
        virtual double value(const astro::SkyDir& dir, double e)const;

        ///@brief integral for the energy limits, in the given direction
        virtual double integral(const astro::SkyDir& dir, double a, double b)const;

        virtual std::string name()const{return std::string("Convolution");}

    private:
        int m_level;

    };


}