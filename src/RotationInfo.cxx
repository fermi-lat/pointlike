/** @file RotationInfo.cxx 
@brief Methods for rotation information

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/RotationInfo.cxx,v 1.8 2008/04/28 03:42:11 burnett Exp $
*/

#include "pointlike/AlignProc.h"

using namespace pointlike;

void RotationInfo::acc(const Hep3Vector& tru, const Hep3Vector& meas, double sigmasq, double alpha) {
    std::vector<double>::iterator il = m_likelihood.begin();
    for(std::vector<HepRotation>::iterator im = m_matrices.begin();il!=m_likelihood.end()&&im!=m_matrices.end();++il,++im) {
        //From psf(u) = (1-1/gamma)*(1+u/gamma)**-gamma
        double u = (tru - (*im).inverse()*meas).mag2();
        u/=sigmasq*2;
        double Fint = m_psf.integral(AlignProc::umax());
        double fpsf = m_psf(u);
        double like = alpha*fpsf/Fint+(1-alpha)/AlignProc::umax();
        *il += log(like);
    }
}

