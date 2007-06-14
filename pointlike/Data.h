/** @file Data.h 
    @brief declaration of the Data wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/users/burnett/pointlike/pointlike/Data.h,v 1.1.1.1 2007/06/10 01:05:26 burnett Exp $
*/


#ifndef pointlike_Data_h
#define pointlike_Data_h

namespace map_tools {
class PhotonMap;
}
#include <string>
#include <vector>

namespace pointlike {
/**
@brief wrapper for PhotonData -- maybe move there
*/
class Data {
public:

    //! constructor loads data from a fits file to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(const std::string& file, int event_type, int source_id=-1)
        ;
    //! constructor loads data from a list of fits files to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(std::vector<std::string> files, int event_type, int source_id=-1);

    //! constructor loads a PhotonMap that was saved in a fits file
    //! @param inputFile the fits file name
    //! @param tablename ["PHOTONMAP"] the fits table name
    Data(const std::string & inputFile, const std::string & tablename="PHOTONMAP");


    //! add  data from the file to current set
    //! @param file Either FT1 or special ROOT file
    //! @param event_type 0 for class A front, etc
    //! @param source_id select given source
    void add(const std::string& file, int event_type=-1, int source_id=-1);

    //! behave like a PhotonMap object
    operator const map_tools::PhotonMap&() const {return *m_data;}

private:
    map_tools::PhotonMap * m_data;
};

}
#endif

