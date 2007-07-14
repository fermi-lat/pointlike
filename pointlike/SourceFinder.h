/** @file SourceFinder.h
@brief declare class SourceFinder

$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/tools/tools/SourceFinder.h,v 1.9 2007/06/03 18:00:23 burnett Exp $
*/

#ifndef pointlike_SourceFinder_h
#define pointlike_SourceFinder_h

//#include "tools/PowerLawFilter.h"
#include "pointlike/Data.h"

#include "map_tools/PhotonMap.h"

#include "astro/SkyDir.h"
#include "astro/HealPixel.h"

#include <vector>
#include <map>


namespace pointlike {
    //----------------------------------------------
    /** @class CanInfo

    */
    class CanInfo // Candidate info.  
    {
    public:
        CanInfo(double value = 0, double sigma = 0, astro::SkyDir dir = astro::SkyDir(0,0)):
          m_value(value), m_sigma(sigma), m_dir(dir), m_2bdeleted(false),
              m_isSource(false) {}

              double value () const {return m_value;}
              double sigma () const {return m_sigma;}
              astro::SkyDir dir () const {return m_dir;}
              double ra() const {return m_dir.ra();}
              double dec() const {return m_dir.dec();}
              bool is2bdeleted () const {return m_2bdeleted;}
              bool isSource () const {return m_isSource;}
              void setDelete () {m_2bdeleted = true;}
              void setSource (bool value = true) {m_isSource = value;}

    private:
        double m_value; ///< TS value.
        double m_sigma; ///< localization sigma
        astro::SkyDir m_dir;
        bool   m_2bdeleted; // True means this is flagged to be deleted later.
        bool   m_isSource;  // True if this corresponds to a confirmed source
    }; 

    class DiffuseCounts; // forward declaration: disabled for now

    /**
    @class SourceFinder
    @brief Find point source Candidates in a data field

    */
    class SourceFinder {
    public:

        typedef std::map<astro::HealPixel, CanInfo> Candidates; 

        //! Region selection
        typedef enum  
        {
            ALL = 0, ///< Select all regions.
            EQUATORIAL = 1,  ///< Select equatorial region only.
            MIDDLE = 2, ///< Select middle region only.
            POLAR = 3, ///< Select polar region only.
        } RegionSelector;

        /** @brief ctor sets up search
        @param datafile the root file containing the data, to be ingested to a PhotonMap
        */
        SourceFinder(const std::string& datafile, DiffuseCounts* dc);

        SourceFinder(const std::string& rootfile, int event_type=-1, int source_id=-1 );

#if 0 // not currently implemented
        /** @brief ctor sets up search
        @param inputFile fits file containing stored PhotonMap structure
        @param tablename fits table name
        */
        SourceFinder(const std::string & inputFile, const std::string & tablename,
            DiffuseCounts* dc);
#endif

        //! add  data from the file to current set
        //! @param event_type 0 for class A front, etc
        //! @param source_id select given source
        void add(const std::string& file, int event_type=-1, int source_id=-1)
        {
            m_data.add(file, event_type, source_id);
        }


        /** @brief return modifiable reference to candidates map
        */
        Candidates & getCandidates() {return m_can;}

        /** @brief
        Analyze range of likelihood significance values for all pixels at a particular level  
        */
        void examineRegion(
            const astro::SkyDir& dir, 
            double radius = 30, 
            double eq_TS_min = 25.0, // Equatorial TS cutoff
            double mid_TS_min = 19.0, // Mid-latitude TS cutoff
            double polar_TS_min = 18.0, // Polar TS cutoff
            int    pix_level = 8, 
            int    count_threshold = 16,
            bool   includeChildren = true, 
            bool   weighted = true,
            bool   background_filter = true,
            int	   skip_TS_levels = 0,
            RegionSelector region = ALL,
            double equator_boundary = 10.0, // abs(b)in degrees for equatorial region < this number.
            double polar_boundary = 40.0 // abs(b)in degrees for polar region > this number.
            ) ;

        /** @brief
        Analyze likelihood significance for a particular direction  
        */

        void checkDir(astro::SkyDir & sd,
            double eq_TS_min = 25.0, // Equatorial TS cutoff,
            double mid_TS_min = 19.0, // Mid-latitude TS cutoff
            double polar_TS_min = 18.0, // Polar TS cutoff
            int    pix_level = 8, 
            int count_threshold = 16,
            bool   background_filter = false,
            int	skip_TS_levels = 2,
            double equator_boundary = 10.0,
            double polar_boundary = 40.0);

        // List selected pixels
        void list_pixels();

        //! Eliminate neighbors within cone
        void prune_neighbors( double degrees);

        //! Eliminate weaker adjacent neighbors
        void prune_neighbors();

        //! summarize results in a ds9 region file
        void createReg(const std::string& filename, double radius = -1., const std::string& color = "white");

        void createTable(const std::string& filename, bool get_background = false, int skip_TS = 0);

        //! allow access to map
        const map_tools::PhotonMap& getMap() {return(m_pmap);}

        //! return vector of candidates, copy of current list

        std::vector<CanInfo> candidateList()const;


    private:
        Data m_data;
        const map_tools::PhotonMap& m_pmap;
        Candidates m_can;
        DiffuseCounts* m_counts;

    };



} // namespace tools

#endif

