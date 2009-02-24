/** @file LeastSquaresFitter.h 
@brief declaration of the least squares fitter class

$Header$
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


    private:

        ///! fitting routine
        ///@brief returns the statistical error in degrees from curvature matrix
        ///@param err error in degrees of initial fit
        double fit(double err);

        std::vector<double> m_fitparams; //parameterization of TS surface
        PointSourceLikelihood* m_psl; // Pointer to PSL object, can be modified
        const static int s_points = 10; //grid points: will crash matrix multiplication if larger than 3
        double m_err; //statistical error value

    };

}//namespace

#endif

