/** @file SimpleTSmap.h
    @brief declare class SimpleTSmap

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SimpleTSmap.h,v 1.4 2008/05/02 23:31:04 burnett Exp $

*/
#ifndef pointlike_SimpleTSmap_h
#define pointlike_SimpleTSmap_h

#include "skymaps/SkySpectrum.h"
#include "astro/SkyDir.h"
#include <map>
#include <vector>

namespace skymaps {
    class BinnedPhotonData;
}

namespace pointlike {

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/** @class SimpleTSmap
    @brief a SkySpectrum that represents the TS calculated at the center of each pixel

 
*/

class SimpleTSmap : public skymaps::SkySpectrum {
public:

    /** @brief ctor
    @param fits_file create from a FITS exposure cube
    @param tablename [SimpleTSmap]
    
    (need to connect to a discription of the effective area function or functions)

    */
    SimpleTSmap(const skymaps::BinnedPhotonData& pmap);
    
    ~SimpleTSmap();

    ///@brief a single energy 
    ///@param e energy in MeV (ignored here)
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    ///@param a lower limit
    ///@param b upper limit
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    virtual std::string name()const;

    /// @brief make a map of TS values for the level only, in the selected region
    void run(astro::SkyDir& center, double radius, double emin=500, bool verbose=false);

    /// @brief save ascii file of level, followed by (index, TS) values
    void save(std::string filename);

    /// @read back saved file
    void restore(std::string filename);

    size_t size()const{return m_tsmap.size();}

    int nside() const {return m_nside;}

    void clear(){m_tsmap.clear();}

    /// @brief return value of a healpix index, or zero if not set
    float operator[](int index)const; 

private:
    const skymaps::BinnedPhotonData& m_pmap;
    std::map<int, std::vector<float> > m_tsmap; ///< the data: a sparse map of a vector of floats
    int m_nside;
    
};


}
#endif

