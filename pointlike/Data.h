/** @file Data.h 
    @brief declaration of the Data wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Data.h,v 1.21 2008/04/04 09:53:54 burnett Exp $
*/


#ifndef pointlike_Data_h
#define pointlike_Data_h

namespace astro {
class SkyDir;
class PointingHistory;
class PointingInfo;
}

namespace skymaps {
class BinnedPhotonData;
}

#include "embed_python/Module.h"
#include "CLHEP/Vector/Rotation.h"
#include <string>
#include <vector>

namespace pointlike {
/***
wrapper for BinnedPhotonData-- maybe move there
*/
class Data {
public:

    //! constructor loads data from a fits FT1 or root file (MeritTuple) to make a BinnedPhotonData
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(const std::string& file, int event_type, double tstart, double tstop, int source_id=-1)
        ;
    //! constructor loads data from a list of fits or root files to make a BinnedPhotonData
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(std::vector<std::string> files, int event_type, double tstart, double tstop,int source_id=-1, 
        std::string ft2file=""
        );
    //! constructor loads a BinnedPhotonData that was saved in a fits file
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

    //! add  gti info from the file to current set
    //! @param file Either FT1 or  MeritTuple ROOT file
    void addgti(const std::string& file);

    //! behave like a BinnedPhotonData object
    operator const skymaps::BinnedPhotonData&() const {return *m_data;}

    //! same as above, for python use
    const skymaps::BinnedPhotonData& map()const{return *m_data;}

    //! @brief define FT2 file to use for rotation
    //! Needed for FT1 file.
    void setHistoryFile(const std::string& history);

    ~Data();

    double minTime()const { return m_start;}
    double maxTime()const { return m_stop; }

    /// @brief summary printout of the BinnedPhotonData object
    void info(std::ostream& out = std::cout);

    static double scale(int i);
    static double set_scale(int i, double s);
    static double class_level();

    //! set corrections to fixed rotation in GLAST frame, default is (0,0,0)
    static CLHEP::HepRotation set_rot(double arcsecx, double arcsecy, double arcsecz);
    static void set_rot(std::vector<double> align);
    static const CLHEP::HepRotation& get_rot();
    static const std::string& historyfile();
    static const astro::PointingInfo& get_pointing(double time);
    static bool inTimeRange(double time);

private:
    void lroot(const std::string& infile);
    static double s_scale[4]; // scale factors
    static int s_class_level; // set to 1,2,3 for transient, source, diffuse

    skymaps::BinnedPhotonData * m_data;
    static CLHEP::HepRotation s_rot;
    static std::string s_ft2file;
    double m_start, m_stop;  ///< overall time
    static astro::PointingHistory* s_history; ///< pointer to optional FT2 info.
    std::vector<std::pair<double,double> > m_gti; ///< time intervals (Good Time Interval)
};

}
#endif

