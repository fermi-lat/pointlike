/** @file ParamOptimization.cxx 
@brief ParamOptimization member functions

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/ParamOptimization.cxx,v 1.12 2008/06/06 17:09:25 burnett Exp $

*/

#include "pointlike/ParamOptimization.h"
#include "skymaps/IParams.h"
#include "TMatrixD.h"
#include <fstream>
//#define PRINT_DEBUG
#define FIT_DEBUG

using namespace pointlike;

#ifdef PRINT_DEBUG
std::ofstream goldout("bracket_output.xls");
#endif

#ifdef FIT_DEBUG
std::ofstream fitout("fit_out.txt");
#endif

namespace {

    std::vector<double> m_params;
    std::vector<double> m_sigmas;
    std::vector<double> m_energy;

    double sigma_param(double a, double b, double c, double d, double e) {
        return sqrt(a*a+c*pow(e/100,b)+d*pow(e/100,2*b));
    }

    double sign(double a, double b) {
        return a*(b<0?-1:1);
    }

    double max(double a, double b) {
        return b>a?b:a;
    }

    void swap(double &a, double&b) {
        double c=a;
        a=b;
        b=c;
    }

    void shift3(double&a, double&b, double&c, const double d) {
        a=b;
        b=c;
        c=d;
    }

    //minimum bracketing from numerical recipes
    void brak(double& ax, double& bx, double &cx, double& fa, double& fb, double& fc, double func(const double)) {
        const double GOLD = 1.618034,GLIMIT=100.0,TINY=1.0e-20;
        double ulim,u,r,q,fu;
        fa = func(ax);
        fb = func(bx);
        if(fb > fa) {
            swap(ax,bx);
            swap(fb,fa);
        }
        cx=bx+GOLD*(bx-ax);
        fc=func(cx);
        while(fb >fc) {
            r=(bx-ax)*(fb-fc);
            q=(bx-cx)*(fb-fa);
            u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(fabs(q-r),TINY),q-r));
            ulim=bx+GLIMIT*(cx-bx);
            if((bx-u)*(u-cx)>0.0) {
                fu=func(u);
                if(fu<fc) {
                    ax=bx;
                    bx=u;
                    fa=fb;
                    fb=fu;
                    return;
                } else if (fu>fb) {
                    cx=u;
                    fc=fu;
                    return;
                }
                u=cx+GOLD*(cx-bx);
                fu=func(u);
            }
            else if((cx-u)*(u-ulim)>0.0) {
                fu=func(u);
                if(fu<fc) {
                    shift3(bx,cx,u,u+GOLD*(u-cx));
                    shift3(fb,fc,fu,func(u));
                }
            }
            else if((u-ulim)*(ulim-cx)>=0.0) {
                u=ulim;
                fu=func(u);
            } else {
                u=cx+GOLD*(cx-bx);
                fu=func(u);
            }
            shift3(ax,bx,cx,u);
            shift3(fa,fb,fc,fu);
        }
    }

    //brent's minimization method from numerical recipes
    double brent(const double ax, const double bx, const double cx, double f(const double), const double tol, double& xmin)
    {
        const int ITMAX=100;
        const double CGOLD=0.3819660;
        const double ZEPS=1e-20;
        double a,b,d=0.0,etemp,fu,fv,fw,fx;
        double p,q,r,tol1,tol2,u,v,w,x,xm;
        double e=0.0;

        a=(ax<cx?ax:cx);
        b=(ax>cx?ax:cx);
        x=w=v=bx;
        fw=fv=fx=f(x);
        for(int iter=0;iter<ITMAX;++iter) {
            xm=0.5*(a+b);
            tol2=2.*(tol1=tol*fabs(x)+ZEPS);
            if(fabs(x-xm)<=(tol2-0.5*(b-a))) {
                xmin=x;
                return fx;
            }
            if(fabs(e)>tol1) {
                r=(x-w)*(fx-fv);
                q=(x-v)*(fx-fw);
                p=(x-v)*q-(x-w)*r;
                q=2.*(q-r);
                if(q>0.0) p=-p;
                q=fabs(q);
                etemp=e;
                e=d;
                if(fabs(p)>=fabs(0.5*q*etemp)||p<=q*(a-x)||p>=q*(b-x)) d=CGOLD*(e=(x>=xm?a-x:b-x));
                else {
                    d=p/q;
                    u=x+d;
                    if(u-a<tol2||b-u<tol2)  d=sign(tol1,xm-x);
                }
            } else {
                d=CGOLD*(e=(x>=xm?a-x:b-x));
            }
            u=(fabs(d)>=tol1?x+d:x+sign(tol1,d));
            fu=f(u);
            if(fu<=fx) {
                if(u>=x) a=x; else b=u;
                shift3(v,w,x,u);
                shift3(fv,fw,fx,fu);
            } else {
                if(u<x) a=u; else b=u;
                if(fv<=fw||w==x) {
                    v=w;
                    w=u;
                    fv=fw;
                    fw=fu;
                } else if (fu<=fv||v==x||v==w) {
                    v=u;
                    fv=fu;
                }
            }
        }
        std::cout << "Too many iterations in brent method!" << std:: endl;
        xmin=x;
        return fx;
    }

    int ncom;
    double (*f1func)(std::vector<double>& );
    std::vector<double> *pcom_p,*xicom_p;

    //line minimization method for numerical recipes
    double f1dim(const double x) {
        int j;
        std::vector<double> xt(ncom);
        std::vector<double> &pcom=*pcom_p,&xicom=*xicom_p;
        for(j=0;j<ncom;++j)
            xt[j]=pcom[j]+x*xicom[j];
        return f1func(xt);
    }

    //line minimization method for numerical recipes
    void linmin(std::vector<double> &p, std::vector<double> &xi, double &fret, double func(std::vector<double>&)) {
        int j;
        const double TOL=1.0e-8;
        double xx,xmin,fx,fb,fa,bx,ax;
        int n=p.size();
        ncom=n;
        pcom_p= new std::vector<double>(n);
        xicom_p=new std::vector<double>(n);
        f1func=func;
        std::vector<double> &pcom =*pcom_p,&xicom=*xicom_p;
        for(j=0;j<n;++j) {
            pcom[j]=p[j];
            xicom[j]=xi[j];
        }
        ax=0.;
        xx=1.;
        brak(ax,xx,bx,fa,fx,fb,f1dim);
        fret=brent(ax,xx,bx,f1dim,TOL,xmin);
        for(j=0;j<n;++j) {
            xi[j]*=xmin;
            p[j]+=xi[j];
        }
        delete xicom_p;
        delete pcom_p;
    }

    //Powell's method for numerical recipes
    void powell(std::vector<double> &p,std::vector<std::vector<double> > &xi, const double ftol,int &iter,double &fret,double func(std::vector<double>&)) {
        const int ITMAX=200;
        const double TINY=1.e-25;
        int i,j,ibig;
        double del,fp,fptt,t;

        int n=p.size();
        std::vector<double> pt(n),ptt(n),xit(n);
        fret=func(p);
        for(j=0;j<n;++j) pt[j]=p[j];
        for(iter=0;;++iter) {
            fp=fret;
            ibig=0;
            del=0;
            for(i=0;i<n;++i) {
                for(j=0;j<n;++j) xit[j]=xi[j][i];
                fptt=fret;
                linmin(p,xit,fret,func);
                if(fptt-fret>del) {
                    del=fptt-fret;
                    ibig=i+1;
                }
            }
            if(2.*(fp-fret)<ftol*(fabs(fp)+fabs(fret))+TINY) {
                return;
            }
            if(iter==ITMAX) std::cout << "Too many iterations of Powell!" << std::endl;
            for(j=0;j<n;++j) {
                ptt[j]=2.*p[j]-pt[j];
                xit[j]=p[j]-pt[j];
                pt[j]=p[j];
            }
            fptt=func(ptt);
            if(fptt<fp) {
                t=2.*(fp-2.*fret+fptt)*(fp-fret-del)*(fp-fret-del)-del*(fp-fptt)*(fp-fptt);
                if(t<0.) {
                    linmin(p,xit,fret,func);
                    for(j=0;j<n;++j) {
                        xi[j][ibig-1]=xi[j][n-1];
                        xi[j][n-1]=xit[j];
                    }
                }
            }
        }
    }

    double chisq(std::vector<double> &params) {
        double chis=0;
        for(int i(0);i<m_sigmas.size();++i) {
            chis+=(sigma_param(params[0],params[1],params[2],params[3],m_energy[i])-m_params[i])
                *(sigma_param(params[0],params[1],params[2],params[3],m_energy[i])-m_params[i])/(m_sigmas[i]*m_sigmas[i]);
        }
        return chis;
    }
}

ParamOptimization::ParamOptimization(const skymaps::BinnedPhotonData &data, const std::vector<astro::SkyDir>& directions, std::ostream *out,int minlevel,int maxlevel) : 
m_data(data)
, m_out(out){
    for(std::vector<astro::SkyDir>::const_iterator it = directions.begin();it!=directions.end();++it) {
        PointSourceLikelihood *pl = new pointlike::PointSourceLikelihood(data,"",*it);
        pl->maximize();
        m_likes.push_back(pl);
    }
}

std::vector<double> ParamOptimization::compute(ParamOptimization::Param p) {
    bool sigma = (p==ParamOptimization::SIGMA);
    if(!sigma) {
        *m_out << "Computing optimum" << (m_data.front().event_class()==0?" front ":" back ") << "gamma values" << std::endl;
        *m_out << std::left << std::setw(10) << "Energy" << "gamma          error          k (~1)    alpha         photons\n";
    } else {
        *m_out << "Computing optimum" << (m_data.front().event_class()==0?" front ":" back ") << "sigma values" << std::endl;
        *m_out << std::left << std::setw(10) << "Energy" << "sigma(deg)     error(deg)     k (~1)    alpha         photons\n";
    }
    *m_out << "***************************************************************************\n";

    int num2look = m_likes.size();
    int timeout = 30;
    double tol = 1e-6;
    //for each level, optimize
    std::vector<double> params;
    std::vector<double> sigmas;
    std::vector<double> energies;
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
        double emin = m_likes.front()->energyList()[iter];
        double emax = m_likes.front()->energyList()[iter+1];
        double ebar = sqrt(emin*emax);
        energies.push_back(ebar);
        params.push_back(maxfactor>0?(sigma?osigma*maxfactor*180/M_PI:osigma*maxfactor):-1);
        sigmas.push_back(t_photons>0?curvature(sigma,iter,osigma*maxfactor)*(sigma?180/M_PI:1):-1);
        *m_out << std::left << std::setw(10) << std::setprecision(3)<<
            (int)ebar << std::setw(15) << (maxfactor>0?(sigma?osigma*maxfactor*180/M_PI:osigma*maxfactor):-1) << 
            std::setw(15) << (t_photons>0?curvature(sigma,iter,osigma*maxfactor)*(sigma?180/M_PI:1):-1) << std::setw(10) << maxfactor << std::setw(15) << (t_photons>0?t_alpha/t_curvature:-1) <<
            std::setw(10) << t_photons << std::endl;
    }
    m_energy=energies;
    m_params=params;
    m_sigmas=sigmas;
    return params;
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
    int npts = 6;        //number of points
    double sep = pow(2.,1./npts);   //spacing between sampling points
    TMatrixD A(2*npts+1,3);
    TMatrixD b(2*npts+1,1);
    double k = 0.5;
    for(int j(0);j<2*npts+1;++j) {
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
#ifdef FIT_DEBUG
        fitout << like << std::endl;
#endif
        k*=sep;
    }
    Double_t* err(0);
    TMatrixD At(3,2*npts+1);
    At.Transpose(A);
    TMatrixD Cv(3,3);
    Cv = At*A;
    Cv = Cv.Invert(err);
    TMatrixD Cf(3,1);
    Cf = Cv*At*b;
    double curv = Cf[0][0];
    return(curv>0?1/sqrt(curv):-1.0);
}



std::vector<double> ParamOptimization::fit_sigma() {
    std::vector<double> d = skymaps::IParams::params(m_data.front().event_class());
    std::cout << std::setprecision(3) << "Initial parameters: " ;
    std::cout << "a= " << d[0] << "\tb= " << d[1] << "\tc= " << d[2] << "\td= " << d[3];
    std::cout << std::endl;
    std::vector<std::vector<double> > unit;
    for(int i(0);i<d.size();++i) {
        std::vector<double> temp(d.size(),0);
        temp[i]=1;
        unit.push_back(temp);
    }
    double fval =0;
    int iters=0;
    powell(d,unit,1e-6,iters,fval,&chisq);
    std::cout << std::setprecision(3) << "Final parameters: " ;
    std::cout << "a= " << d[0] << "\tb= " << d[1] << "\tc= " << d[2] << "\td= " << d[3];
    std::cout << std::endl;
    return d;
}
