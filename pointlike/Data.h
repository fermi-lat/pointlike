/** @file Data.h 
    @brief declaration of the Data wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Data.h,v 1.16 2008/01/29 19:17:52 burnett Exp $
*/


#ifndef pointlike_Data_h
#define pointlike_Data_h

namespace astro {
class SkyDir;
class PointingHistory;
}

namespace skymaps {
class PhotonMap;
}

#include "embed_python/Module.h"
#include "CLHEP/Vector/Rotation.h"
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
    Data(const std::string& file, int event_type, double tstart, double tstop, int source_id=-1)
        ;
    //! constructor loads data from a list of fits or root files to make a PhotonMap
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(std::vector<std::string> files, int event_type, double tstart, double tstop,int source_id=-1, 
        std::string ft2file=""
        );
    //! constructor loads a PhotonMap that was saved in a fits file
    //! @param inputFile the fits file name
    //! @param tablename ["PHOTONMAP"] the fits table name
    Data(const std::string & inputFile, const std::string & tablename="PHOTONMAP");


    //! constructor configure from a python "data" file
    //! @param setup the Python module--the following parameters can be set with a "Data." preceding.
    /** @verbatum
    pixelfile=''    # set to a photon data FITS file, or use a list of FT1 files, If set, ignore alignment, times, output
    files = []      # set to a list of FT1-like FITS files (if no pixefile)
    event_class = -1 # 0, select front only; -1 no selection
    source_id =-1   # -1: all sources -- select according to Monte Carlo source id, if present
    output_pixelfile = '' # set to create an output pixel file (if reading FT1 or ROOT files)
    start_time=0.   # select interval if greater than zero
    stop_time=0.    # "
    history = ''    # optional history or FT2 file, needed to correct for misalignment if readign FT1
    Latalignment=[] # alignment correction angles about x,y,z axes, in arcseconds
    @endverbatum
    */
    Data(const embed_python::Module& setup);

    //! add  data from the file to current set
    //! @param file Either FT1 or  MeritTuple ROOT file
    //! @param event_type 0 for class A front, etc
    //! @param source_id select given source
    void add(const std::string& file, int event_type=-1, int source_id=-1);

    //! behave like a PhotonMap object
    operator const skymaps::PhotonMap&() const {return *m_data;}

    //! same as above, for python use
    const skymaps::PhotonMap& map()const{return *m_data;}

    //! @brief define FT2 file to use for rotation
    //! Needed for FT1 file.
    void setHistoryFile(const std::string& history);

    ~Data();

    double minTime()const { return m_start;}
    double maxTime()const { return m_stop; }
    static double scale(int i);
    static double set_scale(int i, double s);
    static double class_level();

    //! set corrections to fixed rotation in GLAST frame, default is (0,0,0)
    static CLHEP::HepRotation set_rot(double arcsecx, double arcsecy, double arcsecz);
    static void set_rot(std::vector<double> align);
    static const CLHEP::HepRotation& get_rot();

private:
    void lroot(const std::string& infile);
    static double s_scale[4]; // scale factors
    static int s_class_level; // set to 1,2,3 for transient, source, diffuse

    skymaps::PhotonMap * m_data;
    static CLHEP::HepRotation s_rot;
    std::string m_ft2file;
    double m_start, m_stop;  ///< overall time
    astro::PointingHistory* m_history; ///< pointer to optional FT2 info.
    std::vector<std::pair<double,double> > m_gti; ///< time intervals (Good Time Interval)
};

}
#endif

