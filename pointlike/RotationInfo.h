/** @file RotationInfo.h 
@brief declaration of the rotation likelihood class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/RotationInfo.h,v 1.4 2008/01/25 23:09:43 mar0 Exp $
*/
#ifndef pointlike__Rotation_h
#define pointlike__Rotation_h
#include "CLHEP/Vector/Rotation.h"
#include "skymaps/PsfFunction.h"
#include <vector>
using namespace CLHEP;

namespace pointlike {

    //! class RotationInfo
    //! Stores information about the log likelihood of a rotation about each axis in GLAST coordinates
    class RotationInfo {

    public:
        //! constructor makes a grid of rotations {-s_points*arcsec,...0,...+s_points*arcsec}+offset about the X-Y-Z axes  
        //! @param arcsec rotation in arcseconds
        RotationInfo(double arcsec):
        m_psf(2.25)
        {
            for(int x=-s_points;x<=s_points;++x) {
                for(int y=-s_points;y<=s_points;++y) {
                    for(int z=-s_points;z<=s_points;++z) {
                        m_matrices.push_back(HepRotationX(arcsec*x)*HepRotationY(arcsec*y)*HepRotationZ(arcsec*z));
                        m_likelihood.push_back(0);
                    }
                }
            }
        }

        //! returns the loglikelihood of one of the (2*s_points+1)**3 grid points
        //! @param x -s_points,...0,...s_points
        //! @param y -s_points,...0,...s_points
        //! @param z -s_points,...0,...s_points
        double likelihood(int x,int y,int z) {return (fabs(x*1.)<=s_points&&fabs(y*1.)<=s_points&&fabs(z*1.)<=s_points)?m_likelihood[(2*s_points+1)*(2*s_points+1)*(x+s_points)+(2*s_points+1)*(y+s_points)+z+s_points]:0;}

        //! accumulates the loglikelihood
        //! @param tru source direction
        //! @param meas measured direction
        //! @param sigmasq sigma-squared (angular resolution)
        void acc(const Hep3Vector& tru, const Hep3Vector& meas, double sigmasq, double alpha, int level);

        static int points() {return s_points;}
    private:
        skymaps::PsfFunction m_psf;
        const static int s_points = 2; //grid points: will crash matrix multiplication if larger than 3
        std::vector<HepRotation> m_matrices; //matrix storage
        std::vector<double> m_likelihood; //loglikelihood for each rotation

    };

}//namespace

#endif
