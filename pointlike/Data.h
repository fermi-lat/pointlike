/** @file Data.h 
    @brief declaration of the Data wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Data.h,v 1.3 2007/07/14 03:50:54 burnett Exp $
*/


#ifndef pointlike_Data_h
#define pointlike_Data_h

namespace astro {
class SkyDir;
}

namespace map_tools {
class PhotonMap;
}

#include "embed_python/Module.h"
#include <string>
#include <vector>

namespace pointlike {
/***
wrapper for PhotonMap-- maybe move there
*/
class Data {
public:

    //! constructor loads data from a fits FT1 or root file (MeritTuple) to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(const std::string& file, int event_type, int source_id=-1)
        ;
    //! constructor loads data from a list of fits or root files to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(std::vector<std::string> files, int event_type, int source_id=-1, 
        std::string ft2file=""
        );
#if 1 // not currently implemented?
    //! constructor loads a PhotonMap that was saved in a fits file
    //! @param inputFile the fits file name
    //! @param tablename ["PHOTONMAP"] the fits table name
    Data(const std::string & inputFile, const std::string & tablename="PHOTONMAP");
#endif

    Data(embed_python::Module& setup);

    //! add  data from the file to current set
    //! @param file Either FT1 or  MeritTuple ROOT file
    //! @param event_type 0 for class A front, etc
    //! @param source_id select given source
    void add(const std::string& file, int event_type=-1, int source_id=-1);

    //! behave like a PhotonMap object
    operator const map_tools::PhotonMap&() const {return *m_data;}

    //fits image output
    void draw_sky(std::string outputfile, double pixel);
    void draw_region(const astro::SkyDir& dir, std::string outputFile, double pixel, double fov);

    //! same as above, for python use
    const map_tools::PhotonMap& map()const{return *m_data;}

    static double s_scale[4]; // scale factors
    static double set_scale(int i, double s){double t(s_scale[i]);  s_scale[i]=s; return t;}

private:
    map_tools::PhotonMap * m_data;
    std::string m_ft2file;
};

}
#endif

