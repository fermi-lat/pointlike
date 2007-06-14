/** @file PsfFunction.h
@brief declear PsfFunction.

$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/pointlike/pointlike/PsfFunction.h,v 1.1.1.1 2007/06/10 01:05:26 burnett Exp $
*/

#ifndef wavelet_PsfFunction_h
#define wavelet_PsfFunction_h

#include "astro/SkyDir.h"
#include "astro/HealPixel.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace pointlike{
/**
    @class PsfFunction
    @brief Define the power-law psf function 

*/

class PsfFunction 
{
    public:
        /// @param gamma the power-law factor
        PsfFunction(double g=2.0)
            : m_gamma(g)
            , m_norm(1.-1/m_gamma)
        {}

        /** @brief the function.
            @param r is the reference direction.
            @param r_prime points to a position to be evaluated relative to r. 
            @param sigma angular resolution scale factor 
        */
        double operator () (const astro::SkyDir & r, 
            const astro::SkyDir & r_prime,
            double sigma) ;

        ///@return the value as a function of the scaled deviation squared
        double operator()(double u)const;
        double gamma()const{return m_gamma;}

        /// @return integral of the function from 0 to umax
        double integral(double umax)const;

        /// @return integral of the square of the function from 0 to umax
        double integralSquare(double umax)const;

    private:
        double m_gamma;
        double m_norm;
};

}

#endif
