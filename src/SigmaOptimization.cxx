#include "pointlike/SigmaOptimization.h"
#include <fstream>
//#define PRINT_DEBUG

using namespace pointlike;

#ifdef PRINT_DEBUG
std::ofstream goldout("bracket_output.xls");
#endif

SigmaOptimization::SigmaOptimization(const map_tools::PhotonMap &data, const std::vector<astro::SkyDir>& directions, std::ostream *out,int minlevel,int maxlevel,double radius) : 
m_data(data)
, m_minlevel(minlevel)
, m_maxlevel(maxlevel) {
    m_minlevel<data.minLevel()?m_minlevel=data.minLevel():m_minlevel;
    m_maxlevel>data.minLevel()+data.levels()-1?m_maxlevel=data.minLevel()+data.levels()-1:m_maxlevel;
    (*out) << "Computing optimum sigma values" << std::endl;
    (*out) << std::left << std::setw(10) << "Level" << "sigma          error          k (~1)    alpha         photons\n";
    (*out) << "***************************************************************************\n";
    double num2look = directions.size();
    int timeout = 30;
    double tol = 1e-3;
    //for each level, optimize
    for(int iter=minlevel;iter<=maxlevel;++iter) {
        int whileit =0;
        double maxfactor = 0;
        double osigma=0;
        //iterative method for finding best fit (k->1)
        while(maxfactor>=0&&fabs(maxfactor-1.)>tol&&whileit<timeout){
            maxfactor = goldensearch(directions,num2look,iter,radius);
            // sigma = sigma*maxfactor if maxfactor is an appropriate value
            if(maxfactor>0) {
                osigma = PointSourceLikelihood::set_sigma_level(iter,pow(maxfactor,-0.5));
                osigma*= PointSourceLikelihood::set_sigma_level(iter,osigma*pow(maxfactor,-0.5));
            }
            whileit++;
        }
        int t_photons(0);
        double t_curvature(0.);
        double t_alpha(0.);
        //calculate fit statistics for all sources
        for(std::vector<astro::SkyDir>::const_iterator it = directions.begin();it!=directions.end();++it) {
            PointSourceLikelihood ps(m_data,"test",(*it)); //,radius,m_minlevel,m_maxlevel);
            PointSourceLikelihood::iterator ite = ps.find(iter);
            ite->second->maximize();
            if(ite->second->photons()>0) {
                t_curvature += ite->second->kcurvature(maxfactor)/2;
                t_photons += ite->second->photons();
                t_alpha += ite->second->alpha();
            }
        }
        (*out) << std::left << std::setw(10) << 
            iter << std::setw(15) << (maxfactor>0?osigma:-1) << 
            std::setw(15) << t_curvature/directions.size() << std::setw(10) << maxfactor << std::setw(15) << t_alpha/directions.size() <<
            std::setw(10) << t_photons << std::endl;
    }
}

double SigmaOptimization::goldensearch(const std::vector<astro::SkyDir>& directions, int num2look, int level, double radius) {
    // Numerical Recipes section 10.1 - finding a minimum in one dimension with Golden Search
    double tol = 1e-2;
    double C = (3-sqrt(5.))/2;
    int iter2=0;
    double R = 1-C;
    double ax = 1e-2;
    double bx = 1.;
    double cx = 1.e2;
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
    for(std::vector<astro::SkyDir>::const_iterator it=directions.begin(); (it!=directions.end())&& (iter2<num2look);++it,++iter2) {
        pointlike::PointSourceLikelihood ps(m_data, "test", (*it)); //,radius,m_minlevel,m_maxlevel);
        ps.maximize(2);
        pointlike::PointSourceLikelihood::iterator ite = ps.find(level);
        if(ite->second->photons()==0) continue;
        f1+=ite->second->feval(x1);
        f2+=ite->second->feval(x2);
    }
    if(f1==0||f2==0) return -1;
    int k = 1;
    while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
        iter2=0;
        double a1=0;
        double a2=0;
        for(std::vector<astro::SkyDir>::const_iterator it=directions.begin(); (it!=directions.end())&& (iter2<num2look);++it,++iter2) {
            pointlike::PointSourceLikelihood ps(m_data, "test", (*it)); //,radius,m_minlevel,m_maxlevel);
            ps.maximize(2);
            pointlike::PointSourceLikelihood::iterator ite = ps.find(level);
            if(ite->second->photons()==0) continue;
            a1+=ite->second->feval(R*x1+C*x0);
            a2+=ite->second->feval(R*x2+C*x3);

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
