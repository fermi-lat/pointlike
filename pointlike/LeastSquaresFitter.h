/** @file LeastSquaresFitter.h 
@brief declaration of the least squares fitter class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/LeastSquaresFitter.h,v 1.3 2009/04/08 16:55:38 mar0 Exp $
*/
#ifndef pointlike__LeastSquaresFitter_h
#define pointlike__LeastSquaresFitter_h
#include "pointlike/PointSourceLikelihood.h"
#include <vector>
using namespace pointlike;

namespace pointlike {

    //! class LeastSquaresFitter
    //! Tries to find the maximum TS of a PointSourceLikelihood with a least squares fit
    //! to a quadratic function
    //! eg: TS(x,y) = a[0]*x**2 + a[1]*x + a[2]*y**2 + a[3]*y + a[4]*x*y + a[5]
    class LeastSquaresFitter {

    public:
        //! constructor 
        //! @param psl PointSourceLikelihood Object
        //! @param sigma[0] grid size to use, default get from the psl
        LeastSquaresFitter(PointSourceLikelihood& psl, double sigma=0);
        
        ///@brief return statistical error in degrees from curvature matrix 
        double err() {return m_err;}

        ///@brief return fit parameters, a[i]
        std::vector<double> params() {return m_fitparams;}        

        /** @brief return ellipse parameters:
        semi-major axis
        semi-minor axis
        rotation angle
        relative ra center
        relative dec center
        chi-squared fit to surface

         */
        std::vector<double> ellipse() {return m_ellipse;}        


        //@brief returns quadratic function evaluated at x and y coordinates
        double func(std::vector<double> params,double x, double y);

        ///@header for ellipse parameters.
        static std::string header();

    private:

        ///! fitting routine
        ///@brief returns the statistical error in degrees from curvature matrix
        ///@param values  likelihood values evaluated on ring
        ///@param err   error in degrees of initial fit
        double fit(std::vector<double> values, double err);

        ///! TS evaluation routine
        ///@brief returns values around a ring about the SkyDir sd
        std::vector<double> ring(astro::SkyDir& sd, double err);

        ///! determines the new position from parameters in m_ellipse
        astro::SkyDir maxDir(astro::SkyDir& iDir);

        std::vector<double> m_fitparams; //parameterization of TS surface
        PointSourceLikelihood* m_psl; // Pointer to PSL object, can be modified

        std::vector<double> m_ellipse; //0-7 [a,b,ecc,phi,x0,y0,err,chisq,R**2]

        double m_err; //statistical error value

    };

}//namespace

#endif

