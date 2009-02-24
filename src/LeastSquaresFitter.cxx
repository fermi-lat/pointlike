/** @file LeastSquaresFitter.cxx 
@brief Methods for rotation information

$Header$
*/

//#define MIN_DEBUG

#include "pointlike/LeastSquaresFitter.h"
#include "TMatrixD.h"

#ifdef MIN_DEBUG
#include <sstream>
#include <fstream>
#include <string>
#endif

using namespace pointlike;
using namespace astro;

namespace {

    double sqr(double d) {return d*d;}

    double file_num(0.);
    int itermax(3);
}

LeastSquaresFitter::LeastSquaresFitter(PointSourceLikelihood& psl, double sigma):
m_psl(&psl),
m_err(0)
{
    //get current TS
    double iTS = m_psl->TSmap(m_psl->dir());
    double curTS = -1e10;

    if( sigma==0) sigma = m_psl->errorCircle();
    //save initial direction
    SkyDir iDir = m_psl->dir();

    //try to fit 'itermax' times, or until the TS stops improving
    for(int iter(0);iter<itermax && iTS-curTS>0.1;++iter)
    {
        //try to fit with least squares
        m_err=fit(sigma);
        
        //calculate TS at new position
        double newTS = m_psl->TSmap(m_psl->dir());

        //if TS doesn't improve, return to previous position
        if(newTS-iTS<0) {
            m_psl->setDir(iDir,false);
            //m_err=m_psl->errorCircle();
        }

        //update position and TS comparisons
        iDir=m_psl->dir();
        curTS=iTS;
        iTS=newTS;
    }
}

double LeastSquaresFitter::fit(double err)
{
#ifdef MIN_DEBUG
    std::stringstream s;
    s << file_num;
    std::string name("source"+s.str()+".txt");
    const char * names(name.c_str());
    std::ofstream minf(names);
    ++file_num;
#endif

    //keep fit position around
    SkyDir oldDir = m_psl->dir();

    //select new subset of pixels
    m_psl->setDir(oldDir,false);

    //pick 2D coordinate system
    Hep3Vector rand_x = oldDir().orthogonal();
    rand_x=rand_x.unit();
    Hep3Vector rand_y = (oldDir().cross(rand_x)).unit();

    int npts=1;
    int rows=(2*npts+1)*(2*npts+1); //number of grid points = rows*rows
    TMatrixD A(rows,6);
    TMatrixD likes(rows,1);

    double chisq = 0;
    double iTS = m_psl->TSmap(oldDir);

    //create grid of TS values and setup 'A' matrix
    //the equation is of the form:
    //TS(x,y) = a[0]*x**2 + a[1]*x + a[2]*y**2 + a[3]*y + a[4]*x*y + a[5]
    int k =0;
    for(int i(-npts);i<npts+1;++i) {
        for(int j(-npts);j<npts+1;++j) {
            double norm = 1/3.;
            if (abs(i)+abs(j)>=npts)
                norm = sqrt(1.*i*i+1.*j*j)/3.;
            SkyDir newDir(M_PI/180*err*(i*rand_x+j*rand_y)/norm+oldDir());

            A[k][0]=err*err*i*i/(norm*norm);
            A[k][1]=err*i/norm;
            A[k][2]=err*err*j*j/(norm*norm);
            A[k][3]=err*j/norm;
            A[k][4]=err*err*i*j/(norm*norm);
            A[k][5]=1;
            likes[k][0]=m_psl->TSmap(newDir);
            if (abs(i)+abs(j)>=npts) {
                chisq+=sqr(iTS-likes[k][0]-1/sqr(norm))/(1/sqr(norm));
            }
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
    double a0 = Atb[0][0];
    double a1 = Atb[1][0];
    double a2 = Atb[2][0];
    double a3 = Atb[3][0];
    double a4 = Atb[4][0];
    double detA = 2*a0*a2-sqr(a4)/2.;
    double vx = fabs(2*a2/detA);
    double cxy = -a4/detA;
    double vy = fabs(2*a0/detA);
    double xc = (a4*a3-2*a1*a2)/detA;
    double yc = (a4*a1-2*a0*a3)/detA;

    //statistical error is roughly the geometric mean of the semi-major axes
    double nerr = sqrt(sqrt(vx*vy));

    double angle = atan(a4/(a2-a0))/2;
    double a = (a0+a2+(a0-a2)/cos(2*angle))/2;
    double b = (a0+a2+(a2-a0)/cos(2*angle))/2;
    double ecc = sqrt(fabs(a*a-b*b))/fabs(a);
    if(ecc>1) ecc = 1/ecc;


    //set the position to the minimum and grab new set of pixels
    oldDir=SkyDir((xc*rand_x+yc*rand_y)*M_PI/180+oldDir());
    m_fitparams.clear();
    m_fitparams.push_back(oldDir.ra());
    m_fitparams.push_back(oldDir.dec());
    m_fitparams.push_back(sqrt(vx>vy? vx:vy));
    m_fitparams.push_back(sqrt(vx>vy? vy:vx));
    m_fitparams.push_back(angle*180/M_PI);
    m_fitparams.push_back(chisq);


#ifdef MIN_DEBUG
    double cen_val = m_psl->TSmap(oldDir);
    npts = s_points;
    for(int i(-npts);i<npts+1;++i) {
        for(int j(-npts);j<npts+1;++j) {
            double norm = 1/2.;//(i==0&&j==0)?1:sqrt(1.*i*i+1.*j*j);
            SkyDir newDir(M_PI/180*err*(i*rand_x+j*rand_y)/norm+oldDir());
            //setDir(newDir,false);
            double acc = 0;
            acc+=err*err*i*i*a0/norm/norm;
            acc+=err*i*a1/norm;
            acc+=err*err*j*j*a2/norm/norm;
            acc+=err*j*a3/norm;
            acc+=err*err*i*j*a4/norm;
            acc+=Atb[5][0];
            double likes=cen_val-m_psl->TSmap(newDir);
            minf << i/norm << "\t" << j/norm << "\t" << likes  << "\t"<< cen_val-likes-acc << std::endl;;
        }
    }
#endif

    m_psl->setDir(oldDir,false);

    if( m_psl->verbose() ) {
        std::cout << "lsqfit ra, dec, error " 
            << oldDir.ra()<< ", " << oldDir.dec()
            << ", " << nerr
            <<std::endl;
    }

    //scale by root2 since TS = 2*logL
    return nerr;
}
