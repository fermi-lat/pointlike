/** @file LeastSquaresFitter.cxx 
@brief Methods for rotation information

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/LeastSquaresFitter.cxx,v 1.2 2009/02/25 22:29:26 mar0 Exp $
*/

//#define MIN_DEBUG

#include "pointlike/LeastSquaresFitter.h"
#include "TMatrixD.h"

#define SHAPE_DEBUG

#ifdef SHAPE_DEBUG
#include <fstream>
std::ofstream shape("shape.txt");
#endif
#include <iomanip>

using namespace pointlike;
using namespace astro;

namespace {

    double sqr(double d) {return d*d;}
    int s_points = 10; //grid points
    double file_num(0.);
    int itermax(3);
    bool flip(true);
    bool major(true);
    int npts=1;

#ifdef SHAPE_DEBUG
    bool output_flag = true;
#endif
}

LeastSquaresFitter::LeastSquaresFitter(PointSourceLikelihood& psl, double sigma):
m_psl(&psl),
m_err(sigma)
{

    //get current TS
    double iTS = m_psl->TSmap(m_psl->dir());
    double curTS = -1e10;

    if( sigma==0) sigma = m_psl->errorCircle();
    //save initial direction
    SkyDir iDir = m_psl->dir();

    if(m_psl->verbose()) m_psl->out() << "---------------------Least Squares Fit---------------------" << std::endl <<
        "delta       ra        dec        error         ecc        TS"<< std::endl;

    //try to fit 'itermax' times, or until the TS stops improving
    for(int iter(0);iter<itermax && iTS-curTS>0.1;++iter)
    {
        //try to fit with least squares
        std::vector<double> values = ring(iDir,sigma);

        m_err=fit(values,sigma);

        if(m_err>96) {
            m_psl->setDir(iDir);
            if(m_psl->verbose()) m_psl->out() << (m_err<98?">>>>>>>>step error, aborting":">>>>>>>>poor quadratic surface")
                << std::endl;
            break;
        }

        m_psl->setDir(maxDir(iDir),false);

        //calculate TS at new position

        double newTS = m_psl->TSmap(m_psl->dir());

        if(m_psl->verbose()&&m_fitparams.size()>0) m_psl->out() << std::setprecision(4) << sqrt(m_ellipse[4]*m_ellipse[4]+m_ellipse[5]*m_ellipse[5]) <<
            "     " << m_psl->dir().ra() << "   " << m_psl->dir().dec() << "        " << m_ellipse[6] << "       " 
            << m_ellipse[2] << "      " << std::setprecision(1) << newTS << std::endl;

        //if TS doesn't improve, return to previous position
        if(newTS-iTS<0) {
            m_psl->setDir(iDir,false);
            if(m_psl->verbose())
                m_psl->out() << "lsqfit -- Done" << std::endl;
            break;
            //m_err=m_psl->errorCircle();
        }

        //update position and TS comparisons
        iDir=m_psl->dir();
        curTS=iTS;
        iTS=newTS;
    }
#ifdef SHAPE_DEBUG
    if(m_err<97) {
        shape << m_psl->name() << "\t" << iTS << "\t" << m_psl->dir().ra() << "\t" << m_psl->dir().dec();
        for(int i(0); i<m_ellipse.size(); ++i)
            shape << "\t" << m_ellipse[i];
        shape << std::endl; 
    }
#endif 
}

double LeastSquaresFitter::fit(std::vector<double> values, double err)
{

    int rows=(2*npts+1)*(2*npts+1); //number of grid points = rows*rows
    TMatrixD A(rows,6);
    TMatrixD likes(rows,1);

    double chisq = 0;

    int maxstep(10);

    //create grid of TS values and setup 'A' matrix
    //the equation is of the form:
    //TS(x,y) = a[0]*x**2 + a[1]*x + a[2]*y**2 + a[3]*y + a[4]*x*y + a[5]
    int k =0;
    for(int i(-npts);i<npts+1;++i) {
        for(int j(-npts);j<npts+1;++j) {
            double norm = 1/3.;
            if (abs(i)+abs(j)>=npts)
                norm = sqrt(1.*i*i+1.*j*j)/3.;   

            A[k][0]=err*err*i*i/(norm*norm);
            A[k][1]=err*i/norm;
            A[k][2]=err*err*j*j/(norm*norm);
            A[k][3]=err*j/norm;
            A[k][4]=err*err*i*j/(norm*norm);
            A[k][5]=1;
            likes[k][0]=values[k];
            ++k;
        }
    }

    Double_t* merr(0);

    // Solve system of equations for least squares
    // x = (At*A)**-1*(At*b) 
    TMatrixD At(6,rows);
    TMatrixD AtA(6,6);
    TMatrixD Atb(6,1);
    At.Transpose(A);
    AtA = At*A;
    Atb = At*likes;
    AtA.Invert(merr);
    Atb = AtA*Atb;

    //solution to minimum from second order equation
    // [ [2*a0 a4] [a4 2*a2] ] [[x][y]] = [[-a1][-a3]]

    std::vector<double> p_sur;
    for(int i(0);i<6;++i)
        m_fitparams.push_back(Atb[i][0]);

    double a0 = Atb[0][0];
    double a1 = Atb[1][0];
    double a2 = Atb[2][0];
    double a3 = Atb[3][0];
    double a4 = Atb[4][0];


    double detA = 2*a0*a2-sqr(a4)/2.;
    double vx = fabs(2*a2/detA);
    double cxy = -a4/detA;
    double vy = fabs(2*a0/detA);
    double xc = (a4*a3-2*a1*a2)/(2*detA);
    double yc = (a4*a1-2*a0*a3)/(2*detA);

    //statistical error is roughly the geometric mean of the semi-major axes
    double nerr = sqrt(sqrt(vx*vy));

    //keep least squares fit from stepping to far
    if(sqrt(xc*xc+yc*yc)/nerr>maxstep) return 97;

    //check for a possible saddle-point
    if (a0*a2<=0) return 99;

    double x0 = a2-a0;
    double y0 = a4;
    if (fabs(x0)<1e-9) x0=0;
    double z0 = a0+a2;
    if (fabs(y0)<1e-9) y0=0;
    double d0 = sqrt(x0*x0+y0*y0);
    if (x0<0) d0=-d0;
    double angle = x0==0?M_PI/4:atan(y0/x0)/2;
    if(flip) angle = M_PI/2-angle;
    double alpha(fabs(z0-d0)/2),beta(fabs(z0+d0)/2);
    double a = 1/sqrt(alpha);
    double b = 1/sqrt(beta);
    if(major&& b>a) 
    {
        double temp = a;
        a=b;
        b=temp;
        angle += M_PI/2;
    }
    double si = sin(angle);
    double co = cos(angle);
    if(flip) {
        double temp = si;
        si=co;
        co=temp;
    }

    double ko = xc*xc*((co/a)*(co/a)+(si/b)*(si/b)) +
        yc*yc*((co/b)*(co/b)+(si/a)*(si/a)) - 
        2*si*co*xc*yc*(1/a/a-1/b/b);

    double ecc = sqrt(fabs(a*a-b*b))/fabs(a);
    if(ecc>1) ecc = 1/ecc;


    m_ellipse.clear();
    m_ellipse.push_back(a);             //semi-major axis
    m_ellipse.push_back(b);             //semi-minor axis
    m_ellipse.push_back(ecc);           //eccentricity
    m_ellipse.push_back(angle*180/M_PI);//rotation angle
    m_ellipse.push_back(xc);            //relative ra center
    m_ellipse.push_back(yc);            //relative dec center
    m_ellipse.push_back(nerr);          //curvature (geometric mean of axes)
    
    k=0;

    //check goodness of fit
    std::vector<double> chi_vec;
    double ts_sq(0),ts_ave(0);
    for(int i(-npts);i<npts+1;++i) {
        for(int j(-npts);j<npts+1;++j) {
            double norm = 1/3.;
            if (abs(i)+abs(j)>=npts)
                norm = sqrt(1.*i*i+1.*j*j)/3.;
            double chi = func(m_fitparams,i*err/norm,j*err/norm);
            ts_ave += likes[k][0];
            ts_sq += likes[k][0]*likes[k][0];
            chi = chi-likes[k][0];
            chi_vec.push_back(chi);
            chisq+= chi*chi;
            ++k;
        }
    }
    double r_sq = 1-chisq/(ts_sq-ts_ave*ts_ave/9);

    if(r_sq<0) {
        return 99;
    }

    m_ellipse.push_back(chisq); //chi-squared fit to surface
    m_ellipse.push_back(r_sq);  //R-squared value

    return nerr;
}

double LeastSquaresFitter::func(std::vector<double> p, double x, double y) {
    //return quadradtic function
    if(p.size()<6) return -1e40;
    return p[0]*x*x+p[1]*x+p[2]*y*y+p[3]*y+p[4]*x*y+p[5];
}

std::vector<double> LeastSquaresFitter::ring(astro::SkyDir& oldDir, double err) {

    //select new subset of pixels
    //m_psl->setDir(oldDir,false);

    //pick 2D coordinate system
    Hep3Vector rand_x = oldDir().orthogonal();
    rand_x=rand_x.unit();
    Hep3Vector rand_y = (oldDir().cross(rand_x)).unit();
    double iTS = m_psl->TSmap(oldDir);
    std::vector<double> likes;

    //evaluate TS around a ring 3 sigma from center
    for(int i(-npts);i<npts+1;++i) {
        for(int j(-npts);j<npts+1;++j) {
            double norm = 1/3.;
            if (abs(i)+abs(j)>=npts)
                norm = sqrt(1.*i*i+1.*j*j)/3.;  
            SkyDir newDir(M_PI/180*err*(i*rand_x+j*rand_y)/norm+oldDir());
            likes.push_back((i==0)&&(j==0)?0:iTS-m_psl->TSmap(newDir));
        }
    }
    return likes;
}

astro::SkyDir LeastSquaresFitter::maxDir(astro::SkyDir& sd) {
    //pick 2D coordinate system
    Hep3Vector rand_x = sd().orthogonal();
    rand_x=rand_x.unit();
    Hep3Vector rand_y = (sd().cross(rand_x)).unit();

    //set the position to the maximum likelihood
    return SkyDir((m_ellipse[4]*rand_x+m_ellipse[5]*rand_y)*M_PI/180+sd());
}
