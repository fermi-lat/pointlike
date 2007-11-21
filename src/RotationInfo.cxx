/** @file RotationInfo.cxx 
@brief Methods for rotation information

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/RotationInfo.cxx,v 1.1 2007/11/19 20:40:39 mar0 Exp $
*/

#include "pointlike/AlignProc.h"

using namespace pointlike;

namespace {
    double fit_alpha[]={1.,1.,1.,1.,1.,1.,1.,1.,1.};
}

std::vector<double> RotationInfo::s_alphas(fit_alpha,fit_alpha+sizeof(fit_alpha)/sizeof(double));

void RotationInfo::acc(const Hep3Vector& tru, Hep3Vector& meas, double sigmasq, int level) {
    std::vector<double>::iterator il = m_likelihood.begin();
    for(std::vector<HepRotation>::iterator im = m_matrices.begin();il!=m_likelihood.end()&&im!=m_matrices.end();++il,++im) {
        //From psf(u) = (1-1/gamma)*(1+u/gamma)**-gamma
        double u = (tru - (*im)*meas).mag2();
        u/=pointlike::PointSourceLikelihood::sigma_level[level]*pointlike::PointSourceLikelihood::sigma_level[level]*sigmasq*2;
        u = log(RotationInfo::s_alphas[level-6]*pow(1+u/2.25,2.25)/(1-pow(1+AlignProc::s_umax/2.25,1-2.25))+(1-RotationInfo::s_alphas[level-6])/AlignProc::s_umax);
        *il -= u;
    }
}