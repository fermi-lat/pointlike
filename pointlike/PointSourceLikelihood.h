/** @file PointSourceLikelihood.h

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/pointlike/PointSourceLikelihood.h,v 1.58 2009/05/25 18:53:13 burnett Exp $
*/

#ifndef pointlike_PointSourceLikelihood_h
#define pointlike_PointSourceLikelihood_h

#include "pointlike/SimpleLikelihood.h"

#include "skymaps/SkySpectrum.h"
#include "skymaps/Background.h"


#include "astro/SkyDir.h"
#include <iostream>
#include <vector>

namespace embed_python { class Module; }
namespace skymaps { 
    class BinnedPhotonData;
    class CompositeSkySpectrum;
    class BandBackground;
    class SpectralFunction;
}

namespace pointlike {


    /** @class PointSourceLikelihood
    @brief manage a set of SimpleLikelihood objects, one for each energy band / event type

    Note that it is a vector of  SimpleLikelihood objects, as well as a SkySpectrum.


    */
    class PointSourceLikelihood : public  std::vector<SimpleLikelihood*>, public skymaps::SkySpectrum{
    public:


        /** ctor
        @param data   source of the data to fit
        @param name   source name for printout
        @param dir    initial direction
        */
        PointSourceLikelihood(
            const skymaps::BinnedPhotonData& data,
            std::string name,
            const astro::SkyDir& dir);;

        ~PointSourceLikelihood();

        //! fit to signal fractions
        /// @return total TS (selected bands )
        double  maximize();

        //! change the current direction -- resets data and refits
        void setDir(const astro::SkyDir& dir,bool subset=false);

        /// @return the gradient, summed over all bands
        /// @param skip starting band (or number to skip)
        /// @param count [0] number to add: 0 means all
        const CLHEP::Hep3Vector& gradient() const;

        ///@return the curvature, summed over all bamds
        double curvature() const;

        /// @brief Make a table of the fit spectrum, to the ostream
        /// See set_ostream to direct to a file.
        void printSpectrum();

        /// @brief access to a list of energies for the bands
        std::vector<double> energyList()const;

        /// @brief perform localization fit, maximizing joint likelihood
        /// @param skip [0] number of bands to skip
        /// @return error circle radius (deg) or large number corresponding to error condition
        double localize();

        /// @brief invoke localize with skip values from skip1 to skip 2 or until good fit
        double localize(int skip1, int skip2);

        /// @brief calculate -log likelihood
        /// @param model a spectrum to compare with data
        /// @param extended [true] use extended likelihood, assuming a good model for background flux level
        double logLikelihood(const skymaps::SpectralFunction& model, bool extended=true)const;

        std::string name()const{return m_name;}

        const astro::SkyDir& dir()const{return m_dir;}

        double TS()const;
        double TS(int band) const{return (*this).at(band)->TS();}
        double alpha(int band) const{return (*this).at(band)->alpha();}

        double errorCircle()const{return  sqrt(1./curvature())*180/M_PI;}

        void set_ostream(std::ostream& out){m_out=&out;}

        static void set_verbose(bool verbosity=true); //{s_verbose=verbosity;}

        static bool verbose(); //{return s_verbose;}

        ///! implement the SkyFunction interface
        ///@brief return differential value 
        ///@param e energy in MeV
        virtual double value(const astro::SkyDir& dir, double e)const;

        ///@brief use a band to select interval. 
        ///@param dir direction
        ///@param band use band to select energy range, and event class
        virtual double band_value(const astro::SkyDir& dir, const skymaps::Band& band)const;


        ///@brief integral for the energy limits, in the given direction
        virtual double integral(const astro::SkyDir& dir, double a, double b)const;

        /// @brief set all parameters using the embedded python module
        static void setParameters(const embed_python::Module& par);

        /// @brief set radius for individual fits
        static void setDefaultUmax(double umax);

        ///! Set the global diffuse background function, return current value 
        /// @param diffuse any sky spectrum, presumably a DiffuseFunction
        /// @param exposure [1.0] multiplicative factor, presumably the exposure 
        static  const skymaps::SkySpectrum* set_diffuse(const skymaps::SkySpectrum* diffuse, double exposure=1.0);

        ///! Set the global diffuse background function, return current value 
        /// @param diffuse any sky spectrum, presumably a DiffuseFunction
        /// @param exposures vector of exposure opjects, one for each event type 
        static  const skymaps::SkySpectrum* set_diffuse(const skymaps::SkySpectrum* diffuse, 
            std::vector<const skymaps::SkySpectrum*> exposures);

        ///! Set the global diffuse background function
        /// @param background a Background object, with diffuse and effective areas
        static  const skymaps::Background* set_background(const skymaps::Background* background); 

        ///! clear global diffuse background.
        static const skymaps::Background* clear_background(); 

        ///! add a point source fit to the background, for this object only, for subsequent fits
        void addBackgroundPointSource(const PointSourceLikelihood* fit);

        ///! remove all such
        void clearBackgroundPointSource();

        ///! access to background model 
        const skymaps::SkySpectrum * background()const;

        /// @brief set the integration tolerance for the background, return present value
        static double set_tolerance(double tol);

        /// @brief set the range of energy to fit
        static void set_energy_range(double emin, double emax=1e6);
        static void get_energy_range(double& emin, double& emax);

        /// @brief set/get the minimum alpha for TS definition
        static double set_min_alpha(double a);

        /// @brief special display function
        /// @param dir direction
        /// @param energy selects energy band
        /// @mode 0: same as the operator; 1: data; 2: background; 3:fit; 4:residual
        ///     
        double display(const astro::SkyDir& dir, double energy, int mode)const;


        /// @brief evaluate the TS, using current data set, for any direction
        /// @param sdir the direction
        /// @band [-1] return TS for given band index, or all for -1
        double TSmap(const astro::SkyDir& sdir, int band=-1)const;

        static void set_merge(bool merge);
        static bool merge();
        static void set_fitlsq(bool fit);
        static bool fitlsq();


        static void set_maxROI(double roi); ///< set the maximum ROI (degrees, default 180)
        static double maxROI();   ///< return the maximum ROI (degrees)

        static void set_minROI(double roi); ///< set the minimum ROI (degrees, default 0)
        static double minROI();   ///< return the minimum ROI (degrees)

        static void set_maxSig(double r); ///< set the maximum value of the psf width to be used in localization/TS (degrees, default 0.25)
        static double maxSig(); ///< return the maximum sigma

        static void set_conf(bool conf);
        static bool conf();

        void setup(const skymaps::BinnedPhotonData& data);

        std::ostream& out()const{return *m_out;}




    private:
        iterator begin_skip(int skip);
        const_iterator begin_skip(int skip)const;

        std::string m_name;
        astro::SkyDir m_dir; ///< common direction
        double m_dir_sigma;  ///< error circle from fit (radians)
        double m_TS;         ///< total TS value
        std::ostream * m_out;
        mutable CLHEP::Hep3Vector m_gradient; ///< current gradient

        skymaps::CompositeSkySpectrum * m_background;  ///< background spectrum to use
        std::vector<skymaps::BandBackground*> m_backlist; ///< list of wrapped background objects to delete

        static const skymaps::Background* s_diffuse; ///< global diffuse used by all PSL objects

        static double s_emin, s_emax, s_minalpha, s_TSmin, s_tolerance, 
            s_maxstep; //
        static int s_skip1, s_skip2, s_itermax, s_verbose;
        static int s_merge, s_fitlsq, s_conf; ///<least squares localization flag
        static double s_maxROI; ///< the maximum ROI, in degrees: if set, a limit on umax
        static double s_minROI; ///< the minimum ROI, in degrees: if set, a limit on umax
        static double s_maxSig;

    };

    /** @class PSLdisplay

    */
    class PSLdisplay :  public skymaps::SkySpectrum {
    public:
        PSLdisplay(const PointSourceLikelihood & psl, int mode);
        virtual double value(const astro::SkyDir& dir, double e)const;

        ///@brief integral for the energy limits, in the given direction
        virtual double integral(const astro::SkyDir& dir, double a, double b)const;

        std::string name()const;
    private:
        const PointSourceLikelihood& m_psl;
        int m_mode;

    };

    /** class TSmap
    Create a SkyFunction with the TS defined by the data set for 

    */
    class TSmap : public astro::SkyFunction {
    public:
        ///! @brief ctor to create function based on given fit
        ///! @param psl 
        ///! @param band [-1] index of band to use, -1 for sum
        ///! if doing sum, subtract value at dir()
        TSmap(const PointSourceLikelihood& psl, int band=-1);

        /// @brief ctor to compure TS independently at the given position
        /// @param data data to use
        /// @param band [-1] index of band to use, -1 for sum
        TSmap(const skymaps::BinnedPhotonData& data, int band=-1);

        /// @brief set the point source to be used as background for a data scan
        /// Note that it requires that a diffuse background be defined
        void setPointSource(const PointSourceLikelihood& psl);

        virtual ~TSmap(){};
        virtual double operator()(const astro::SkyDir& sdir)const;
    private:
        const skymaps::BinnedPhotonData * m_data;
        const PointSourceLikelihood* m_psl;
        int m_band;
        double m_offset;
    };
}

#endif
