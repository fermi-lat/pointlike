/** @file SourceFinder.h
@brief declare class SourceFinder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SourceFinder.h,v 1.25 2008/05/02 23:31:04 burnett Exp $
*/

#ifndef pointlike_SourceFinder_h
#define pointlike_SourceFinder_h

//#include "tools/PowerLawFilter.h"
#include "pointlike/Data.h"
#include "skymaps/BinnedPhotonData.h"
#include "pointlike/PointSourceLikelihood.h"

#include "astro/SkyDir.h"

#include "embed_python/Module.h"

#include <vector>
#include <map>


namespace pointlike {
    //----------------------------------------------
    /** @class CanInfo

    */


    class CanInfo // Candidate info.  
    {
    public:
        CanInfo(double value = 0, double sigma = 0, 
            astro::SkyDir dir = astro::SkyDir(0,0)):
        m_value(value),
            m_sigma(sigma), 
            m_dir(dir), 
            m_2bdeleted(false),
            m_isSource(false),
            m_weighted_count(0),
            m_hasStrongNeighbor(false),
            m_strongNeighbor( 0) //healpix::HealPixel(0,0))
        {
            m_values.clear();
            m_photons.clear();
            m_sigalph.clear();
        }

        double value () const {return m_value;}
        double values (int level) {return m_values[level];}
        double photons (int level) {return m_photons[level];}
        double sigalph (int level)  {return m_sigalph[level];}
        double sigma () const {return m_sigma;}
        astro::SkyDir dir () const {return m_dir;}
        double ra() const {return m_dir.ra();}
        double dec() const {return m_dir.dec();}
        bool is2bdeleted () const {return m_2bdeleted;}
        bool isSource () const {return m_isSource;}
        double pl_slope () const {return m_pl_slope;}
        double pl_constant () const {return m_pl_constant;}
        double pl_confidence () const {return m_pl_confidence;}
        int weighted_count () const {return m_weighted_count;}
        int skipped () const {return m_skipped;}
        bool hasStrongNeighbor() const {return m_hasStrongNeighbor;}
        int strongNeighbor() const {return m_strongNeighbor;}

        void setDelete () {m_2bdeleted = true;}
        void setSource (bool value = true) {m_isSource = value;}
        void set_total_value (double value = 0.0) {m_value = value;}
        void set_sigma (double value = 0.0) {m_sigma = value;}
        void set_dir (astro::SkyDir value = astro::SkyDir(0,0)) {m_dir = value;}
        void setValue(int level, double val) {m_values[level] = val;}
        void setPhotons(int level, double photons) {m_photons[level] = photons;}
        void setSigalph(int level, double sigalph) {m_sigalph[level] = sigalph;}
        void set_pl_slope (double value = 0.0) {m_pl_slope = value;}
        void set_pl_constant (double value = 0.0) {m_pl_constant = value;}
        void set_pl_confidence (double value = 0.0) {m_pl_confidence = value;}
        void set_weighted_count (int value = 0) {m_weighted_count = value;}
        void set_skipped (int value = 0) {m_skipped = value;}
        void setHasStrongNeighbor (bool value = true) {m_hasStrongNeighbor = value;}
        void setStrongNeighbor (int value) {m_strongNeighbor = value;}

    private:
        double m_value; ///< TS value.
        std::map<int,double> m_values;
        std::map<int,double> m_photons;
        std::map<int,double> m_sigalph;
        double m_sigma; ///< localization sigma
        astro::SkyDir m_dir;
        bool   m_2bdeleted; // True means this is flagged to be deleted later.
        bool   m_isSource;  // True if this corresponds to a confirmed source
        double m_pl_slope; // Slope of power law fit.
        double m_pl_constant; // b from (y = mx + b) for power law fit.
        double m_pl_confidence; // Confidence of the power law fit. 1 == perfect fit.
        int m_weighted_count; // weighted count of photons in enclosing pixel.  level of enclosing pixel is determined by pointfind_setup.py
        int m_skipped; // number of candidates rejected before this one was accepted.  Count is reset each time a candidate is accepted.
        bool m_hasStrongNeighbor;  // Is there a stronger nearby candidate?
        int m_strongNeighbor;  // Location of strongest nearby candidate.

    }; 

    class DiffuseCounts; // forward declaration: disabled for now

    /**
    @class SourceFinder
    @brief Find point source Candidates in a data field

    */
    class SourceFinder {
    public:

        SourceFinder(const pointlike::Data& data,  const embed_python::Module & Mod);
        typedef std::map<int, CanInfo> Candidates;
        typedef std::map<int, pointlike::PointSourceLikelihood > LikelihoodMap;
        typedef std::map<double, CanInfo> Prelim; // Preliminary candidates

        //! Region selection
        typedef enum  
        {
            ALL = 0, ///< Select all regions.
            EQUATORIAL = 1,  ///< Select equatorial region only.
            MIDDLE = 2, ///< Select middle region only.
            POLAR = 3, ///< Select polar region only.
        } RegionSelector;

#if 0 // not currently implemented

        /** @brief ctor sets up search
        @param datafile the root file containing the data, to be ingested to a BinnedPhotonData
        */
        SourceFinder(const std::string& datafile, DiffuseCounts* dc, const Module & Mod);

        SourceFinder(const std::string& rootfile, int event_type=-1, int source_id=-1, 
            const Module & Mod );

        /** @brief ctor sets up search
        @param inputFile fits file containing stored BinnedPhotonData structure
        @param tablename fits table name
        */
        SourceFinder(const std::string & inputFile, const std::string & tablename,
            DiffuseCounts* dc, const Module & Mod);

        //! add  data from the file to current set
        //! @param event_type 0 for class A front, etc
        //! @param source_id select given source
        void add(const std::string& file, int event_type=-1, int source_id=-1)
        {
            m_data.add(file, event_type, source_id);
        }
#endif


        /** @brief return modifiable reference to candidates map
        */
        Candidates & getCandidates() {return m_can;}

        /** @brief Analyze range of likelihood significance values for all pixels at a particular level  
        */
        void examineRegion(void) ;

        /** @brief Analyze likelihood again for candidates that had a strong neighbor.  Strong neighbor is added to background.  
        */
        void reExamine(void) ;

        /** @brief Analyze likelihood significance for a particular direction  
        */

        void checkDir(astro::SkyDir & sd,
            double eq_TS_min = 25.0, // Equatorial TS cutoff,
            double mid_TS_min = 19.0, // Mid-latitude TS cutoff
            double polar_TS_min = 18.0, // Polar TS cutoff
            int    pix_level = 8, 
            int count_threshold = 16);

        // List selected pixels
        void list_pixels();

        //! Eliminate candidates that don't meet power law tests
        void prune_power_law(void);

        //! Group nearby candidates to facilitate further examination
        //! Any candidate that is within "prune_radius" of another is matched with its strongest neighbor
        void group_neighbors(void);

        //! Eliminate neighbors within cone
        void prune_neighbors(void);

        //! Eliminate weaker adjacent neighbors
        void prune_adjacent_neighbors();

        //! summarize results in a ds9 region file
        void createReg(const std::string& filename, double radius = -1.,
            const std::string& color = "white");

        void createTable(const std::string& filename, bool get_background = false, int skip_TS = 0);

        //! allow access to map
        const skymaps::BinnedPhotonData& getMap() {return(m_pmap);}

        //! return vector of candidates, copy of current list

        std::vector<CanInfo> candidateList()const;

        //! write contents of CanInfo to fits file
        void createFitsFile(const std::string & outputFile,
            const std::string & tablename="PNTFIND", bool clobber= true) const;

        //! write a reg file
        void createRegFile(std::string filename, std::string color="white")const;


        //! run the current set of steps
        void run();

    private:
        const skymaps::BinnedPhotonData& m_pmap;
        Candidates m_can;
        DiffuseCounts* m_counts;
        const embed_python::Module & m_module;

    };



} // namespace tools

#endif

