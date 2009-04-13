/** @file SourceList.h
@brief declaration of classes Source and SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SourceList.h,v 1.14 2009/03/17 23:33:12 burnett Exp $
*/

#ifndef pointlike_SourceList_h
#define pointlike_SourceList_h

#include "astro/SkyDir.h"

#include <string>
#include <vector>
#include <list>
#include <map>

namespace skymaps{
    class BinnedPhotonData;
}

namespace pointlike{

    class PointSourceLikelihood;


    /** @class Source
    @brief manage a single point source
    */
    class Source {
    public:

        /** @brief ctor
        */
        Source(const std::string& name, const astro::SkyDir& seed_dir, double TS=0);
        Source(const std::string& name, double ra, double dec, double TS=0);

        Source():m_name("default"){}; ///< default ctor
        ~Source();

        double localize();

        const std::string& name()const { return m_name;}
        const astro::SkyDir& dir()const{return m_dir;}
        const astro::SkyDir& seed_dir()const{return m_seed_dir;}

        /// @brief if localized, the distance moved (in degrees)
        double moved()const;

        /// @brief dump a line of info
        void info(std::ostream& out=std::cout)const;
        static void header(std::ostream& out=std::cout); ///< header line for info

        double TS()const{return m_TS;}
        double& TS(){return m_TS;} 
        double seedTS()const{return m_seedTS;}
        double sigma()const{return m_sigma;}
        PointSourceLikelihood* fit(){return m_fit;}
        const PointSourceLikelihood* fit()const{return m_fit;}
        const Source* neighbor()const{return m_neighbor;}

        void set_neighbor(const Source* other){m_neighbor = other;}
        const std::vector<double>& fit_params()const{return m_fitparams;}

        /// @brief return an array of delta TS values, nominally a circle at 3 sigma
        std::vector<double> TScircle(double radius=3.)const;
        
    private:
        std::string m_name;
        astro::SkyDir m_dir; ///< current dir
        astro::SkyDir m_seed_dir; ///< initial seed direction
        double m_seedTS;     ///< initial seed tS
        PointSourceLikelihood * m_fit; ///< pointer to the fitter
        double m_TS;
        double m_sigma;  ///< localization rms error
        const Source * m_neighbor; ///< pointer to strong neighbor (zero if none)
        void setup();  ///< called by ctors
        std::vector<double> m_fitparams; ///< parameters of lsqfit
    };

    /** @class SourceList
    @brief manage a list of Source objects
    */
    class SourceList: public std::list<pointlike::Source> {
    public:
        /** @brief ctor, initialized with map

        values in the vector are (ra, dec, [ts])
        */
        SourceList(std::map<std::string, std::vector<double> > input_list);

        /** @brief ctor, from text file

        Format of the text file is:
          name ra dec [ts]

        lines starting with '#' are ignored.
        The optional TS field is for initial sorting.
        */
        SourceList(const std::string& filename, bool verbose=true, std::ostream* log=&std::cout);

        SourceList(){}; ///< default ctor

        void sort_TS(); ///< sort the list in decreasing TS order
        void sort_ra(); ///< sort the list in increasing ra order

        /// @brief refit all sources, in TS order, taking nearby sources into account
        void refit(); 

        void filter_TS(double threshold=0); ///< remove entries with TS < this level

        /// @brief formatted dump to an open stream
        /// Format is consistent with text file ctor
        void dump(std::ostream& out=std::cout)const;

        ///@brief formatted dump to nameed file
        void dump(const std::string& outfilename)const;
        ///@brief write a VOtable
        void dump_xml(std::ostream& out, std::string name="pointfit")const;

        void createRegFile(std::string filename, std::string color, double tsmin)const;

        bool verbose()const{return m_verbose;}
        static void set_log(std::ostream* logstream);

        static void set_data(const skymaps::BinnedPhotonData* data);
        static const skymaps::BinnedPhotonData* data();

        /// @brief set the radius used to select nearby sources for inclusion in a fit
        static double set_group_radius(double value);

    private:
        bool m_verbose;
        static const skymaps::BinnedPhotonData* s_data;

    };

}

#endif


