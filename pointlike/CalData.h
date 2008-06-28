/** @file CalData.h 
    @brief declaration of the CalData wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/CalData.h,v 1.23 2008/05/02 23:31:38 burnett Exp $
*/


#ifndef pointlike_CalData_h
#define pointlike_CalData_h

namespace astro {
}

namespace skymaps {
class BinnedPhotonData;
}

#include "skymaps/PhotonMap.h"
#include "embed_python/Module.h"
#include "CLHEP/Vector/Rotation.h"
#include <string>
#include <vector>

namespace pointlike {
/***
wrapper for BinnedPhotonData-- maybe move there
*/
class CalData {
public:

    //! constructor loads data from a fits FT1 or root file (MeritTuple) to make a BinnedPhotonData
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source (MC only, of course)
    CalData(const std::string& file, int event_type, double emin, double emax, int source_id=-1);
    
    //! behave like a BinnedPhotonData object
    operator const skymaps::BinnedPhotonData&() const {return *m_data;}

    //! same as above, for python use
    const skymaps::BinnedPhotonData& map()const{return *m_data;}

   ///@brief access to static that defines class level cut for ROOT
    static int class_level();

    ~CalData();

private:
    void lroot(const std::string& infile,int type, double emin, double emax);
    static double s_scale[4]; // scale factors
    static int s_class_level; // set to 1,2,3 for transient, source, diffuse

    skymaps::BinnedPhotonData * m_data;
};

}
#endif

