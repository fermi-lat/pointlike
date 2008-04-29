/** @file PointSourceLikelihood.h

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/PointSourceLikelihood.h,v 1.29 2008/04/28 03:42:10 burnett Exp $
*/

#ifndef tools_PointSourceLikelihood_h
#define tools_PointSourceLikelihood_h

#include "pointlike/SimpleLikelihood.h"

#include "skymaps/SkySpectrum.h"

#include "astro/SkyDir.h"
#ifdef OLD
#include "healpix/HealPixel.h"
#endif
#include <iostream>
#include <vector>

namespace embed_python { class Module; }
namespace skymaps { 
    class BinnedPhotonData;
    class CompositeSkySpectrum;}

namespace pointlike {


/** @class PointSourceLikelihood
@brief manage a set of SimpleLikelihood objects, one for each energy band / event type

Note that it is a map of the SimpleLikelihood objects, with the key an index to a corresponding Band 
object, as maintained by the corresponding BinnedPhotonData.


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
    /// @return total TS
    double  maximize();

    //! change the current direction -- resets data and refits
    void setDir(const astro::SkyDir& dir,bool subset=false);

    /// @return the gradient, summed over all bands
    const CLHEP::Hep3Vector& gradient() const;

    ///@return the curvature, summed over all bamds
    double curvature() const;

    /// @brief 
    void printSpectrum();

    /// @brief perform localization fit, maximizing joint likelihood
    /// @param skip [0] number of bands to skip
    /// @return error circle radius (deg) or large number corresponding to error condition
    double localize(int skip);

    /// @brief invoke localize with skip values from skip1 to skip 2 or until good fit
    double localize(int skip1, int skip2);

    /// @brief localate with iteration to refit the levels, using parameters set in ctor
    double localize();

    std::string name()const{return m_name;}

    const astro::SkyDir& dir()const{return m_dir;}

    double TS()const { return m_TS; } 

    /// @param level
    /// @return the invidual TS for the level
    double levelTS(int level)  { return (*this)[level]->TS();}

    double logL(int level){ return (*this)[level]->operator()();}

    double errorCircle()const{return  sqrt(1./curvature())*180/M_PI;}

    void set_ostream(std::ostream& out){m_out=&out;}

    void set_verbose(bool verbosity=true){m_verbose=verbosity;}

    bool verbose()const{return m_verbose;}

    ///! implement the SkyFunction interface
    ///@brief return differential value 
    ///@param e energy in MeV
    virtual double value(const astro::SkyDir& dir, double e)const;

    ///@brief integral for the energy limits, in the given direction
    virtual double integral(const astro::SkyDir& dir, double a, double b)const;

    /// @brief set all parameters using the embedded python module
    static void PointSourceLikelihood::setParameters(const embed_python::Module& par);

    /// @brief set radius for individual fits
    static void setDefaultUmax(double umax);

    /// @brief access to the sigma (radians) used for the individual SimpleLikelihood objects
    double sigma(int level)const;

    ///! Set the global diffuse background function, return current value 
    /// @param diffuse any sky spectrum, presumably a DiffuseFunction
    /// @param exposure [1.0] multiplicative factor, presumably the exposure 
    static  skymaps::SkySpectrum* set_diffuse( skymaps::SkySpectrum* diffuse, double exposure=1.0);

    ///! add a point source fit to the background for subsequent fits
    void addBackgroundPointSource(const PointSourceLikelihood* fit);

    ///! remove all such
    void clearBackgroundPointSource();

    ///! access to background model 
    const skymaps::SkySpectrum * background()const;

    /// @brief set the integration tolerance for the background, return present value
    static double set_tolerance(double tol);
    /// @brief special display function
    /// @param dir direction
    /// @param energy selects energy band
    /// @mode 0: same as the operator; 1: data; 2: background; 3:fit; 4:residual
    ///     
    double display(const astro::SkyDir& dir, double energy, int mode)const;


private:
    void setup(const skymaps::BinnedPhotonData& data);
    std::string m_name;
    astro::SkyDir m_dir; ///< common direction
    double m_dir_sigma;  ///< error circle from fit (radians)
    double m_TS;         ///< total TS value
    bool m_verbose;
    std::ostream * m_out;
    std::ostream& out()const{return *m_out;}
    mutable CLHEP::Hep3Vector m_gradient; ///< current gradient

    static skymaps::SkySpectrum* s_diffuse; ///< global diffuse used by all PSL objects
    

    skymaps::CompositeSkySpectrum * m_background;  ///< background spectrum to use
    
    static double s_emin, s_minalpha, s_TSmin, s_tolerance, 
        s_maxstep; //
    static int s_skip1, s_skip2, s_itermax, s_verbose;

 
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

}
#endif

