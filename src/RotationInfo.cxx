/** @file RotationInfo.cxx 
@brief Methods for rotation information

$Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/src/RotationInfo.cxx,v 1.9 2008/06/17 05:26:04 mar0 Exp $
*/

#include "pointlike/AlignProc.h"

using namespace pointlike;

void RotationInfo::acc(const CLHEP::Hep3Vector& tru, const CLHEP::Hep3Vector& meas, double sigmasq, double alpha) {
    std::vector<double>::iterator il = m_likelihood.begin();
    for(std::vector<CLHEP::HepRotation>::iterator im = m_matrices.begin();il!=m_likelihood.end()&&im!=m_matrices.end();++il,++im) {
        //From psf(u) = (1-1/gamma)*(1+u/gamma)**-gamma
        double u = (tru - (*im).inverse()*meas).mag2();
        u/=sigmasq*2;
        double Fint = m_psf.integral(AlignProc::umax());
        double fpsf = m_psf(u);
        double like = alpha*fpsf/Fint+(1-alpha)/AlignProc::umax();
        *il += log(like);
    }
}

