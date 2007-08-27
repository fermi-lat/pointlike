/** @file CalData.h 
    @brief declaration of the CalData wrapper class for use with the alignment method, returns a PhotonMap with Ra and Dec equivalent
           to the first two rotation angles @thetax and @thetay

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/CalData.h,v 1.3 2007/07/14 03:50:54 burnett Exp $
*/


#ifndef pointlike_CalData_h
#define pointlike_CalData_h

namespace map_tools {
class PhotonMap;
}

namespace astro {
class SkyDir;
}
#include <string>
#include <vector>

namespace pointlike {
/***
wrapper for PhotonMap-- maybe move there
*/
class CalData {
public:

    //! constructor loads data from a fits FT1 or root file (MeritTuple) to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    CalData(const std::string& file, std::vector<astro::SkyDir> sources, int start, int stop, int event_type, int source_id=-1)
        ;
    //! constructor loads data from a list of fits or root files to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    CalData(std::vector<std::string> files, std::vector<astro::SkyDir> sources, int start, int stop, int event_type, int source_id=-1, 
        std::string ft2file=""
        );
#if 0 // not currently implemented?
    //! constructor loads a PhotonMap that was saved in a fits file
    //! @param inputFile the fits file name
    //! @param tablename ["PHOTONMAP"] the fits table name
    CalData(const std::string & inputFile, const std::string & tablename="PHOTONMAP");
#endif


    //! add  data from the file to current set
    //! @param file Either FT1 or  MeritTuple ROOT file
    //! @param event_type 0 for class A front, etc
    //! @param source_id select given source
    void add(const std::string& file, int event_type=-1, int source_id=-1);

    //! behave like a PhotonMap object
    operator const map_tools::PhotonMap&() const {return *m_data;}

    //! same as above, for python use
    const map_tools::PhotonMap& map()const{return *m_data;}

    static double s_scale[4]; // scale factors
    static double set_scale(int i, double s){double t(s_scale[i]);  s_scale[i]=s; return t;}

private:
    map_tools::PhotonMap * m_data;
    std::string m_ft2file;
    std::vector<astro::SkyDir> m_sources;
    int m_start;
    int m_stop;
};

}
#endif

