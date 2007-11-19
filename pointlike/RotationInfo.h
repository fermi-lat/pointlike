/** @file RotationInfo.h 
@brief declaration of the rotation likelihood class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/pointlike/RotationInfo.h,v 1.3 2007/07/14 03:50:54 burnett Exp $
*/

#include "CLHEP/Vector/Rotation.h"
#include <vector>
using namespace CLHEP;

namespace pointlike {

    //! class RotationInfo
    //! Stores information about the log likelihood of a rotation about each axis in GLAST coordinates
    class RotationInfo {

    public:
        //! constructor makes a grid of rotations +2,+1,0,-1,-2 arcmin +offset about the X-Y-Z axes  
        //! @param arcmin rotation in arcminutes
        RotationInfo(double arcmin, double offx, double offy, double offz){
            for(int x=-2;x<=2;++x) {
                for(int y=-2;y<=2;++y) {
                    for(int z=-2;z<=2;++z) {
                        m_matrices.push_back(HepRotationX(arcmin*x+offx)*HepRotationY(arcmin*y+offy)*HepRotationZ(arcmin*z+offz));
                        m_likelihood.push_back(0);
                    }
                }
            }
        }

        //! returns the loglikelihood of one of the 125 grid points
        //! @param x +2,+1,0,-1,-2
        //! @param y +2,+1,0,-1,-2
        //! @param z +2,+1,0,-1,-2
        double likelihood(int x,int y,int z) {return (fabs(x*1.)<=2&&fabs(y*1.)<=2&&fabs(z*1.)<=2)?m_likelihood[25*(x+2)+5*(y+2)+z+2]:0;}

        //! accumulates the loglikelihood
        //! @param tru source direction
        //! @param meas measured direction
        //! @param sigmasq sigma-squared (angular resolution)
        void acc(Hep3Vector& tru, Hep3Vector& meas, double sigmasq, int level);

        static std::vector<double> s_alphas;
        static void setalphas(std::vector<double>& newalpha) {s_alphas=newalpha;}

    private:
        std::vector<HepRotation> m_matrices; //matrix storage
        std::vector<double> m_likelihood; //loglikelihood for each rotation

    };

}//namespace