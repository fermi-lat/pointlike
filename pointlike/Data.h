/** @file Data.h 
    @brief declaration of the Data wrapper class

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/Data.h,v 1.38 2009/02/24 23:07:10 kerrm Exp $
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
class PhotonBinner;
}


#include "embed_python/Module.h"
#include "CLHEP/Vector/Rotation.h"
#include <string>
#include <vector>
#include <iostream>

namespace pointlike {
class Alignment; //  forward declaration
    
    /***
wrapper for BinnedPhotonData-- maybe move there
*/
class Data {
public:

    //! constructor loads data from a fits FT1 or root file (MeritTuple) to make a BinnedPhotonData
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source (MC only, of course)
    Data(const std::string& file, int event_type, double tstart, double tstop, int source_id=-1)
        ;
    //! constructor loads data from a list of fits or root files to make a BinnedPhotonData
    //! @param event_type 0 for class A front, etc, -1 for all
    //! @param source_id select given source
    Data(std::vector<std::string> files, int event_type=-1, double tstart=0, double tstop=0,int source_id=-1, 
        std::string ft2file=""
        );



    //! constructor loads a BinnedPhotonData that was saved in a fits file
    //! @param inputFile the fits file name
    //! @param tablename ["BANDS"] the fits table name
    Data(const std::string & inputFile, const std::string & tablename="BANDS");


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
    alignment_data  # name of an ascii file, with columns  (MET, x, y, z); MET is start of applicability, rotation numbers in arcsec. 
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
    //! @return true if ok, false if GTI not included in range
    bool addgti(const std::string& file);

    //! behave like a BinnedPhotonData object
    operator const skymaps::BinnedPhotonData&() const {return *m_data;}

    //! same as above, for python use
    const skymaps::BinnedPhotonData& map()const{return *m_data;}
    skymaps::BinnedPhotonData& map() {return *m_data;}


    ~Data();

    double minTime()const { return m_start;}
    double maxTime()const { return m_stop; }

    /// @brief summary printout of the BinnedPhotonData object
    void info(std::ostream& out = std::cout);

    /// @brief combine similar bands after read in
    void combine_bands(); 

    /**@brief Write  to a fits file
    @param outputFile Fully qualified fits output file name
    @param clobber Whether to delete an existing file first 
    */
    void write(const std::string & outputFile, bool clobber = true) const;

    //! @brief define FT2 file to use for rotation
    //! Needed for FT1 file, if applying alignment correction
    static void setHistoryFile(const std::string& history);

    ///@brief access to static that defines minimum class level cut
    static int class_level();
    static void set_class_level(int cut);

    ///@brief access to boolean switch to use MC_ENERGY instead of ENERGY
    static bool use_mc_energy();
    static void set_use_mc_energy(bool use_it);

    //! set corrections to fixed rotation in GLAST frame, default is (0,0,0)
    static void set_rot(double arcsecx, double arcsecy, double arcsecz);
    static void set_rot(std::vector<double> align);
    
    ///@brief set the time-dependent alignment constants from a file
    ///@param filename name of a file of constants, or "default"
    static void set_alignment(const std::string& filename);

    static const CLHEP::HepRotation& get_rot(double time);

    static const std::string& historyfile();
    static const astro::PointingInfo& get_pointing(double time);
    static bool inTimeRange(double time);

    /// @brief change default binning: must be done before loading data files
    static void setEnergyBins(const std::vector<double>& bins);
    /// @brief change energy binning class: must be done before loading data files
    static void setPhotonBinner(skymaps::PhotonBinner* binner);

    /// @brief return  value of the cut
    static double zenith_angle_cut();
    static void set_zenith_angle_cut(double cut);

    /// @brief return  value of FOV cut
    static double theta_cut();
    static void set_theta_cut(double cut);

    /// @brief set Earth coordinate system
    static void useEarthCoordinates();
    


private:
    skymaps::BinnedPhotonData * m_data; ///< manage this object
    double m_start, m_stop;  ///< overall time interval (default 0 for no cut)
    std::ostream* m_log;

    std::ostream& log(){return m_log!=0? *m_log: std::cout; }

    void load_filelist(const std::vector<std::string>& filelist, int event_class=-1, int source_id=-1);

    void lroot(const std::string& infile,int event_class);
    static double s_scale[4]; // scale factors
    static int s_class_level; // set to 1,2,3 for transient, source, diffuse

    static pointlike::Alignment* s_alignment; ///< LAT alignment manager
    static std::string s_ft2file;
    static astro::PointingHistory* s_history; ///< pointer to optional FT2 info.
    static double s_zenith_angle_cut;     ///< static value for cut on earth photons
    static double s_theta_cut;        ///< static value for theta, or FOV cut
    static bool s_use_mc_energy;     ///< static value for mc energy flag
};

}
#endif

