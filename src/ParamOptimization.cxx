/** @file ParamOptimization.cxx 
@brief ParamOptimization member functions

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/ParamOptimization.cxx,v 1.10 2008/05/02 23:31:04 burnett Exp $

*/

#include "pointlike/ParamOptimization.h"
#include "skymaps/EnergyBinner.h"
#include "TMatrixD.h"
#include <fstream>
//#define PRINT_DEBUG

using namespace pointlike;

#ifdef PRINT_DEBUG
std::ofstream goldout("bracket_output.xls");
#endif

ParamOptimization::ParamOptimization(const skymaps::BinnedPhotonData &data, const std::vector<astro::SkyDir>& directions, std::ostream *out,int minlevel,int maxlevel) : 
m_data(data)
, m_out(out){
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

    int num2look = m_likes.size();
    int timeout = 30;
    double tol = 1e-3;
    //for each level, optimize

    for(int iter=0;iter<m_likes.front()->size();++iter) {
        int whileit =0;
        double maxfactor = 0;
        double osigma=sigma?(*(m_likes.front()))[iter]->sigma():(*(m_likes.front()))[iter]->gamma();
        std::vector<double> alphas;
        for(std::vector<PointSourceLikelihood*>::iterator it = m_likes.begin();it!=m_likes.end();++it) {
            alphas.push_back((*(*it))[iter]->alpha());
        }
        //iterative method for finding best fit (k->1)

        while(maxfactor>=0&&fabs(maxfactor-1.)>tol&&whileit<timeout){
            maxfactor = goldensearch(num2look,iter,sigma);
            // param = param*maxfactor if maxfactor is an appropriate value
            int i(0);
            for(std::vector<PointSourceLikelihood*>::iterator it = m_likes.begin();it!=m_likes.end();++it,++i) {
                if(maxfactor>0) {
                    sigma?(*(*it))[iter]->setsigma(osigma/sqrt(maxfactor)):(*(*it))[iter]->setgamma(osigma*maxfactor);
                }
                (*(*it))[iter]->setalpha(alphas[i]);
            }
            whileit++;
        }
        int t_photons(0);
        double t_curvature(0.);
        double t_alpha(0.);
        //calculate fit statistics for all sources
        for(std::vector<PointSourceLikelihood*>::iterator it = m_likes.begin();it!=m_likes.end();++it) {
            SimpleLikelihood* ite = (*(*it))[iter];
            if(ite->photons()>0) {
                double t_sa=ite->sigma_alpha();
                t_curvature+=1/(t_sa*t_sa);
                t_photons += ite->photons();
                t_alpha += ite->alpha()/(t_sa*t_sa);
            }
        }
        *m_out << std::left << std::setw(10) << std::setprecision(3)<<
            iter << std::setw(15) << (maxfactor>0?osigma*maxfactor:-1) << 
            std::setw(15) << (t_photons>0?curvature(sigma,iter,osigma*maxfactor):-1) << std::setw(10) << maxfactor << std::setw(15) << (t_photons>0?t_alpha/t_curvature:-1) <<
            std::setw(10) << t_photons << std::endl;

        std::cout << std::endl;
    }
}

double ParamOptimization::goldensearch(int num2look, int band, bool sigma) {

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
        SimpleLikelihood* ite = (*(*it))[band];
        if(ite->photons()==0) continue;
        if(sigma) {
            f1+=ite->feval(x1);
            f2+=ite->feval(x2);
        }
        else{
            f1+=ite->geval(x1);
            f2+=ite->geval(x2);
        }
    }
    if(f1==0||f2==0) return -1.0;
    int k = 1;
    while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
        iter2=0;
        double a1=0.;
        double a2=0.;

        for(std::vector<PointSourceLikelihood*>::iterator it=m_likes.begin(); (it!=m_likes.end())&& (iter2<num2look);++it,++iter2) {

            SimpleLikelihood* ite = (*(*it))[band];
            if(ite->photons()==0) continue;
            if(sigma) {
                a1+=ite->feval(R*x1+C*x0);
                a2+=ite->feval(R*x2+C*x3);
            }
            else {
                a1+=ite->geval(R*x1+C*x0);
                a2+=ite->geval(R*x2+C*x3);
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

double ParamOptimization::curvature(bool sigma,int band,double val)
{
    //Find least squares solution to a quadratic about the minimum likelihood position
    int npts = 15;        //number of points
    double sep = 0.003;   //spacing between sampling points
    TMatrixD A(npts,3);
    TMatrixD b(npts,1);
    for(int j(0);j<npts;++j) {
        double k = 1.-(npts/2-j)*sep;
        double like(0);
        A[j][0]= (k*val)*(k*val);
        A[j][1]= k*val;
        A[j][2]= 1;
        for(std::vector<PointSourceLikelihood*>::iterator it = m_likes.begin();it!=m_likes.end();++it) {
            if(sigma) {
                like+=(*(*it))[band]->feval(k);
            }
            else {
                like+=(*(*it))[band]->geval(k);
            }
        }
        b[j][0]=like;
    }
    Double_t* err(0);
    TMatrixD At(3,npts);
    At.Transpose(A);
    TMatrixD Cv(3,3);
    Cv = At*A;
    Cv = Cv.Invert(err);
    TMatrixD Cf(3,1);
    Cf = Cv*At*b;
    double curv = Cf[0][0];
    return(curv>0?1/sqrt(curv):-1.0);
}