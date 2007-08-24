/** @file Data.h 
    @brief declaration of the Data wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Data.h,v 1.4 2007/08/12 04:18:57 burnett Exp $
*/


#ifndef pointlike_Data_h
#define pointlike_Data_h

namespace map_tools {class PhotonMap;}
namespace astro { class SkyDir; }
#include <string>
#include <vector>
#include "embed_python/Module.h"


namespace pointlike {
/***
wrapper for PhotonMap-- maybe move there
*/
class Data {
public:
    /** constructor uses python module, look for following attributes:
        pixelfile  if set, use it to find table PHOTONMAP, otherwise, expect files
        files      list of files, either FITS or ROOT
        event_class [-1]
        source_id [-1]

    **/
    Data( embed_python::Module& setup);

    //! constructor loads data from a fits FT1 or root file (MeritTuple) to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(const std::string& file, int event_type, int source_id=-1)
        ;
    //! constructor loads data from a list of fits or root files to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(std::vector<std::string> files, int event_type, int source_id=-1 );

    //! constructor loads a PhotonMap that was saved in a fits file
    //! @param inputFile the fits file name
    //! @param tablename ["PHOTONMAP"] the fits table name
    Data(const std::string & inputFile, const std::string & tablename="PHOTONMAP");

    //! add  data from the file to current set
    //! @param file Either FT1 or  MeritTuple ROOT file
    //! @param event_type 0 for class A front, etc
    //! @param source_id select given source
    void add(const std::string& file, int event_type=-1, int source_id=-1);

    //! behave like a PhotonMap object
    operator const map_tools::PhotonMap&() const {return *m_data;}

    //! same as above, for python use
    const map_tools::PhotonMap& map()const{return *m_data;}

    
    //! create FITS image file using the data
    //! @param dir center
    //! @param dir file to write
    //! @param pixel
    //! @param fov

    void draw_region(const astro::SkyDir& dir, std::string outputFile, double pixel, double fov);
    void draw_sky(std::string outputfile, double pixel);
    ~Data();

    static double s_scale[4]; // scale factors
    static double set_scale(int i, double s){double t(s_scale[i]);  s_scale[i]=s; return t;}

private:
    map_tools::PhotonMap * m_data;
};

}
#endif

