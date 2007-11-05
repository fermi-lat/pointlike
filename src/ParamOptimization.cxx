/** @file ParamOptimization.cxx 
    @brief ParamOptimization member functions

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/ParamOptimization.cxx,v 1.1 2007/11/01 21:50:22 mar0 Exp $

*/

#include "pointlike/ParamOptimization.h"
#include <fstream>
//#define PRINT_DEBUG

using namespace pointlike;

#ifdef PRINT_DEBUG
std::ofstream goldout("bracket_output.xls");
#endif

ParamOptimization::ParamOptimization(const map_tools::PhotonMap &data, const std::vector<astro::SkyDir>& directions, std::ostream *out,int minlevel,int maxlevel) : 
m_data(data)
, m_minlevel(minlevel)
, m_maxlevel(maxlevel)
, m_out(out){
    m_minlevel<data.minLevel()?m_minlevel=data.minLevel():m_minlevel;
    m_maxlevel>data.minLevel()+data.levels()-1?m_maxlevel=data.minLevel()+data.levels()-1:m_maxlevel;
    for(std::vector<astro::SkyDir>::const_iterator it = directions.begin();it!=directions.end();++it) {
        PointSourceLikelihood *pl = new pointlike::PointSourceLikelihood(data,"",*it);
        pl->maximize();
        m_likes.push_back(pl);
    }
}

void ParamOptimization::compute(ParamOptimization::Param p) {
    bool sigma = (p==ParamOptimization::SIGMA);
    if(!sigma) {
        *m_out << "Computing optimum gamma values" << std::endl;
        *m_out << std::left << std::setw(10) << "Level" << "gamma          error          k (~1)    alpha         photons\n";
    } else {
        *m_out << "Computing optimum sigma values" << std::endl;
        *m_out << std::left << std::setw(10) << "Level" << "sigma          error          k (~1)    alpha         photons\n";
    }
    *m_out << "***************************************************************************\n";

    double num2look = m_likes.size();
    int timeout = 30;
    double tol = 1e-3;
    m_alphas.clear();
    //for each level, optimize
    for(int iter=m_minlevel;iter<=m_maxlevel;++iter) {
        int whileit =0;
        double maxfactor = 0;
        double osigma=0;
        //iterative method for finding best fit (k->1)
        while(maxfactor>=0&&fabs(maxfactor-1.)>tol&&whileit<timeout){
            maxfactor = goldensearch(num2look,iter,sigma);
            // param = param*maxfactor if maxfactor is an appropriate value
            if(maxfactor>0) {
                if(sigma) {
                    osigma = PointSourceLikelihood::set_sigma_level(iter,pow(maxfactor,-0.5));
                    osigma*= PointSourceLikelihood::set_sigma_level(iter,osigma*pow(maxfactor,-0.5));
                }else {
                    osigma = PointSourceLikelihood::set_gamma_level(iter,maxfactor);
                    osigma*= PointSourceLikelihood::set_gamma_level(iter,osigma*maxfactor);
                }
            }
            for(std::vector<PointSourceLikelihood*>::iterator it = m_likes.begin();it!=m_likes.end();++it) {
                (*it)->recalc(iter);
                (*it)->maximize();
            }
            whileit++;
        }
        int t_photons(0);
        double t_curvature(0.);
        double t_alpha(0.);
        //calculate fit statistics for all sources
        for(std::vector<PointSourceLikelihood*>::iterator it = m_likes.begin();it!=m_likes.end();++it) {
            PointSourceLikelihood::iterator ite = (*it)->find(iter);
            ite->second->maximize();
            if(ite->second->photons()>0) {
                double curv = sigma?ite->second->kcurvature(maxfactor)/2:ite->second->gcurvature(maxfactor)/2;
                t_curvature += pow(curv,-2);
                t_photons += ite->second->photons();
                t_alpha += ite->second->alpha()*pow(curv,-2.);
            }
        }
        m_alphas.push_back(t_alpha/t_curvature);
        *m_out << std::left << std::setw(10) << 
            iter << std::setw(15) << (maxfactor>0?osigma:-1) << 
            std::setw(15) << pow(t_curvature,-0.5) << std::setw(10) << maxfactor << std::setw(15) << t_alpha/t_curvature <<
            std::setw(10) << t_photons << std::endl;
    }
}

double ParamOptimization::goldensearch(int num2look, int level, bool sigma) {
    // Numerical Recipes section 10.1 - finding a minimum in one dimension with Golden Search
    double tol = 1e-2;
    double C = (3-sqrt(5.))/2;
    int iter2=0;
    double R = 1-C;
    double ax = 1e-2;
    double bx = 1.;
    double cx = 1.e1;
    double x0 = ax;
    double x3 = cx;
    double x1,x2;
    double xmin,fmin;
    if (fabs(cx-bx) > fabs(bx-ax)) {
        x1 = bx;
        x2 = bx + C*(cx-bx);
    }else {
        x2 = bx;
        x1 = bx - C*(bx-ax);
    }
    double f1 = 0;
    double f2 = 0;
    for(std::vector<PointSourceLikelihood*>::iterator it=m_likes.begin(); (it!=m_likes.end())&& (iter2<num2look);++it,++iter2) {
        (*it)->maximize();
        pointlike::PointSourceLikelihood::iterator ite = (*it)->find(level);
        if(ite->second->photons()==0) continue;
        if(sigma) {
            f1+=ite->second->feval(x1);
            f2+=ite->second->feval(x2);
        }
        else{
            f1+=ite->second->geval(x1);
            f2+=ite->second->geval(x2);
        }
    }
    if(f1==0||f2==0) return -1;
    int k = 1;
    while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
        iter2=0;
        double a1=0;
        double a2=0;
        for(std::vector<PointSourceLikelihood*>::iterator it=m_likes.begin(); (it!=m_likes.end())&& (iter2<num2look);++it,++iter2) {

            pointlike::PointSourceLikelihood::iterator ite = (*it)->find(level);
            if(ite->second->photons()==0) continue;
            if(sigma) {
                a1+=ite->second->feval(R*x1+C*x0);
                a2+=ite->second->feval(R*x2+C*x3);
            }
            else {
                a1+=ite->second->geval(R*x1+C*x0);
                a2+=ite->second->geval(R*x2+C*x3);
            }
        }
#ifdef PRINT_DEBUG
        goldout << R*x1+C*x0 << "\t" << a1 << std::endl;
        goldout << R*x2+C*x3 << "\t" << a2 << std::endl;
#endif
        if (f2 < f1){
            x0 = x1;
            x1 = x2;
            x2 = R*x1 + C*x3;   // x2 = x1+c*(x3-x1)
            f1 = f2;
            f2 = a2;
        }
        else {
            x3 = x2;
            x2 = x1;
            x1 = R*x2 + C*x0;   // x1 = x2+c*(x0-x2)
            f2 = f1;
            f1 = a1;
        }
        k = k+1;
    }
    if (f1 < f2) {
        xmin = x1;
        fmin = f1;
    }else {
        xmin = x2;
        fmin = f2;
    }
    return xmin;
}

