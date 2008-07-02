/** @file PointSourceLikelihood.h
@brief declaration of classes Source and SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/SourceList.h,v 1.2 2008/06/29 21:07:17 burnett Exp $
*/

#ifndef pointlike_SourceList_h
#define pointlike_SourceList_h

#include "astro/SkyDir.h"

#include <string>
#include <vector>
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
        double sigma()const{return m_sigma;}
        PointSourceLikelihood* fit(){return m_fit;}
        const PointSourceLikelihood* fit()const{return m_fit;}
        const Source* neighbor()const{return m_neighbor;}

        void set_neighbor(const Source* other){m_neighbor = other;}
        
    private:
        std::string m_name;
        astro::SkyDir m_dir; ///< current dir
        astro::SkyDir m_seed_dir; ///< initial seed direction
        double m_TS;
        PointSourceLikelihood * m_fit; ///< pointer to the fitter
        double m_sigma;  ///< localization rms error
        const Source * m_neighbor; ///< pointer to strong neighbor (zero if none)
    };

    /** @class SourceList
    @brief manage a list of Source objects
    */
    class SourceList: public std::vector<pointlike::Source> {
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
        SourceList(const std::string& filename);

        SourceList(){}; ///< default ctor

        void sort_TS();

        /// @brief refit all sources, in TS order, taking nearby sources into account
        void refit(); 

        /// @brief formatted dump to an open stream
        /// Format is consistent with text file ctor
        void dump(std::ostream& out=std::cout)const;

        ///@brief formatted dump to nameed file
        void dump(const std::string& outfilename)const;

        void createRegFile(std::string filename, std::string color, double tsmin)const;

        static void set_data(const skymaps::BinnedPhotonData* data);
        static const skymaps::BinnedPhotonData* data();

    private:
        static const skymaps::BinnedPhotonData* s_data;
    };

}

#endif

