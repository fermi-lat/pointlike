/** @file PointSourceLikelihood.h
@brief declaration of classes Source and SourceList

$Header$
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
        Source(const std::string& name, const astro::SkyDir& sdir, double TS=0);

        ~Source();

        double localize();

        const std::string& name()const { return m_name;}
        const astro::SkyDir& dir()const{return m_dir;}
        double TS()const{return m_TS;}
        double sigma()const{return m_sigma;}
        PointSourceLikelihood* fit(){return m_fit;}
        const PointSourceLikelihood* fit()const{return m_fit;}
        const Source* neighbor()const{return m_neighbor;}

        void set_neighbor(const Source* other){m_neighbor = other;}
        
    private:
        std::string m_name;
        astro::SkyDir m_dir;
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
        */
        SourceList(std::map<std::string, std::vector<double> > input_list);

        /** @brief ctor, from text file
        */
        SourceList(const std::string& filename);

        void sort_TS();

        /// @brief refit all sources, in TS order, taking nearby sources into account
        void refit(); 

        /// @brief formatted dump
        void dump(std::ostream& out=std::cout)const;

        void createRegFile(std::string filename, std::string color, double tsmin)const;

        static void set_data(const skymaps::BinnedPhotonData* data);
        static const skymaps::BinnedPhotonData* data();

    private:
        static const skymaps::BinnedPhotonData* s_data;
    };

}

#endif

