/** @file PointSourceLikelihood.h

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/PointSourceLikelihood.h,v 1.2 2007/06/14 20:01:45 burnett Exp $
*/

#ifndef tools_PointSourceLikelihood_h
#define tools_PointSourceLikelihood_h

#include "pointlike/SimpleLikelihood.h"
#include "astro/SkyDir.h"
#include <iostream>
#include <map>

namespace map_tools { class PhotonMap; }

namespace pointlike {


    /** @class PointSourceLikelihood
    @brief manage a set of SimpleLikelihood objects, one for each energy band

    Note that it is a map of the SimpleLikelihood objects, with the key the Healpix level,
    usually starting at 6, for 0.9 degree bins.


    */
    class PointSourceLikelihood : public  std::map<int, SimpleLikelihood*> {
    public:
        /** ctor
        @param data   source of the data to fit
        @param name   descriptive name for printout
        @param dir    initial direction
        @param radius selection radius (degrees) [5.]
        @param minlevel minimum healpix level [6]
        @param maxlevel maximum healpix level [13]
        */
        PointSourceLikelihood(const map_tools::PhotonMap& data,
            std::string name,
            const astro::SkyDir& dir,  
            double radius=5.0,
            int minlevel=6, int maxlevel=13);

        ~PointSourceLikelihood();

        //! fit to signal fraction for each level
        /// @return total TS
        /// @param skip levels to skip
        double  maximize(int skip=0);

        //! change the current direction
        void setDir(const astro::SkyDir& dir);

        //! set the background density for all the layers
        void setBackgroundDensity(const std::vector<double>& density);

        /// @return the gradient, summed over all levels, skiping skip
        Hep3Vector gradient(int skip=0) const;

        ///@return the curvature, summed over all levels
        double curvature(int skip=0) const;

        /// @brief 
        void printSpectrum();

        /// @brief perform localization fit, maximizing joint likelihood
        /// @param skip [0] number of levels to skip
        /// @return error circle radius (deg) or negative if bad or no fit.
        double localize(int skip=0);

        /// @brief invoke localize with skip values from skip1 to skip 2 or until good fit
        double localize(int skip1, int skip2);

        /// @brief localate with iteration to refit the levels
        double localize(int skip1, int skip2, int itermax, double TSmin=5.0 );

        const std::string& name()const{return m_name;}

        const astro::SkyDir& dir()const{return m_dir;}

        double TS()const { return m_TS; } 

        /// @param level
        /// @return the invidual TS for the level
        double levelTS(int level)  { return (*this)[level]->TS();}

        double logL(int level){ return (*this)[level]->operator()();}

        double errorCircle(int skip=0)const{return  sqrt(1./curvature(skip))*180/M_PI;}

        void set_ostream(std::ostream& out){m_out=&out;}

        void set_verbose(bool verbosity=true){m_verbose=verbosity;}

        bool verbose()const{return m_verbose;}

        /// @brief set radius for individual fits
        static void setDefaultUmax(double umax){ SimpleLikelihood::s_defaultUmax=umax; }

        static std::vector<double> gamma_level;
        static std::vector<double> sigma_level;

        static double set_gamma_level(int level, double v){
            double t = gamma_level[level]; gamma_level[level]=v; return t;}

        static double set_sigma_level(int level, double v){
            double t = sigma_level[level]; sigma_level[level]=v; return t;}

    private:
        std::string m_name;
        astro::SkyDir m_dir; ///< common direction
        double m_dir_sigma;  ///< error circle from fit (radians)
        double m_TS;         ///< total TS value

        bool m_verbose;
        std::ostream * m_out;
        std::ostream& out()const{return *m_out;}

        // the data to feed each guy, extracted from the database
        std::map<int, std::vector<std::pair<astro::HealPixel,int> > >m_data_vec;
    };

}
#endif

