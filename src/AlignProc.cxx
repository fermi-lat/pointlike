/** 
Data Processing file, operates on a given Photon

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/AlignProc.cxx,v 1.5 2008/01/25 23:09:44 mar0 Exp $

*/

#include "pointlike/AlignProc.h"
#include "skymaps/PhotonMap.h"
#include "pointlike/Draw.h"
#include "TMatrixD.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TSystem.h"

#include <float.h>
#include <stdexcept>
#include <fstream>

using namespace pointlike;
using namespace CLHEP;
//#define PRINT
#ifdef PRINT
std::ofstream outfile("covar.txt");
std::ofstream surfacexyz("surfacexyz.txt");
#endif

namespace{

    //extra factor (kluge for now) for correcting scaled deviation
    double E6 =100;
    int minlevel = 8;
    int maxlevel = 13;
    void ShowPercent(int sofar, int total, int found)
    {
        static int toskip(50), skipped(0);
        if(++skipped<toskip) return; skipped=0;
        static int lastpercent(-1);
        int percent( static_cast<int>(100 * sofar / total +0.5) );
        if( percent==lastpercent) return;
        lastpercent=percent;
        char   s[50];
        sprintf(s, "%d%%, %d found.", percent, found);
        std::cout << s;
        if (sofar < total)
        {
            for (unsigned int j = 0; j < strlen(s); ++j)
                std::cout << "\b";
        }
        else
            std::cout << std::endl;
    }

    inline static void do_load (void)
    {
        static bool first = true;
        if( first) {
            gSystem->Load("libTree");
            first=false;
        }
    }

    bool isFinite(double val) {
        using namespace std; // should allow either std::isfinite or ::isfinite
#ifdef WIN32 
        return (_finite(val)!=0);  // Win32 call available in float.h
#else
        return (isfinite(val)!=0); // gcc call available in math.h 
#endif
    }

    std::string root_names[] = {"FT1Ra", "FT1Dec", 
        "CTBBestEnergy", "EvtElapsedTime",
        "FT1ConvLayer","PtRaz", 
        "PtDecz","PtRax",
        "PtDecx", "CTBClassLevel"
    };

    double scale[] = {1.,1.86,1.,1.};
}

//setup initally with no rotation
HepRotation AlignProc::s_hr(0,0,0);
int AlignProc::s_classlevel = 2;
int AlignProc::s_umax = 50;

AlignProc::AlignProc(std::vector<astro::SkyDir>& sources, std::vector<std::vector<double> >& alphas, std::vector<std::string>& files,double arcsecs, double start, double stop): 
m_photons(6,0),
m_start(start),
m_stop(stop),
m_skydir(sources),
m_alphas(alphas),
m_arcsec(arcsecs),
m_roti(arcsecs*M_PI/648000)//offx*M_PI/648000,offy*M_PI/648000,offz*M_PI/648000)
{
    for(std::vector<std::string>::const_iterator it = files.begin();it!=files.end();++it) {
        //either load through ROOT or cfitsio
        if( (*it).find(".root") != std::string::npos) {
            loadroot(*it);
        }else {
            loadfits(*it);
        }
    }
};


void AlignProc::loadroot(const std::string& file) {
    do_load(); //load necessary libraries, but only on first file
    TFile *tf = new TFile(file.c_str(),"READ");
    TTree *tt = static_cast<TTree*>(tf->Get("MeritTuple"));
    tt->SetBranchStatus("*", 0); // turn off all branches
    //turn on appropriate branches
    for(unsigned int j(0); j< sizeof(root_names)/sizeof(std::string); j++){
        tt->SetBranchStatus(root_names[j].c_str(), 1);
    }
    int entries = static_cast<int>(tt->GetEntries());
    std::vector<float> row;
    bool flag(true);
    //for each entry
    for(int i(0);i<entries&&(flag||m_start<0);++i) {

        tt->GetEvent(i);
        //for each
        for( int j(0); j< sizeof(root_names)/sizeof(std::string); j++){
            TLeaf * tl = tt->GetLeaf(root_names[j].c_str());
            if(0==tl) {
                tl = tt->GetLeaf(("_" + root_names[j]).c_str());
                if(0==tl) {
                    tt->Print();
                    throw std::invalid_argument(std::string("Tuple: could not find leaf ")+root_names[j]);
                }
            }
            float v = tl->GetValue();
            row.push_back(isFinite(v)?v:-1e8);
        }
        pointlike::AlignProc::Photona p = events(row);
        if(row[3]>=m_start) {
            add(p); 
            ShowPercent(i,entries,i);
        }
        if(row[3]>m_stop) {
            flag=false;
        }
        row.clear();
    }
}

void AlignProc::loadfits(const std::string& /*file */) {
    //todo : process fits file in same manner as ROOT

}

pointlike::AlignProc::Photona AlignProc::events(std::vector<float>& row) {
    float ra(0), dec(0), energy(0); // photon info
    float raz(0), decz(0), rax(90), decx(0); // sc orientation: default orthogonal
    double time(0);
    int event_class(99);
    int flag =1;
    int source_id(0);
    int class_level(0);
    for(unsigned int i = 0;i<row.size();++i) {
        if(row[i]<-1e7) flag=0;
    }
    if(flag) {
        time = row[3];
        event_class = static_cast<int>(row[4]);
        event_class = event_class>4? 0 : 1;  // front/back map to event class 0/1
        energy = row[2];
        class_level = static_cast<int>(row[9]);
        source_id = static_cast<int>(row[10]);
        if( class_level < AlignProc::class_level()) event_class=99;
        else{
            ra = row[0];
            dec = row[1];
            raz = row[5];
            decz = row[6];
            rax = row[7];
            decx = row[8];
            energy/=scale[event_class];
        }
    }
    pointlike::AlignProc::Photona p(astro::SkyDir(ra, dec), energy, time, event_class ,
        astro::SkyDir(raz,decz),astro::SkyDir(rax,decx));
    astro::Photon ap = p.transform(AlignProc::s_hr.inverse());
    return pointlike::AlignProc::Photona(ap.dir(),ap.energy(),ap.time(),ap.eventClass(),astro::SkyDir(raz,decz),astro::SkyDir(rax,decx));
}

int AlignProc::add(pointlike::AlignProc::Photona& p){
    int added(0);
    //unused int cl = p.eventClass();
    //above healpix level 8 (mostly signal photons)
    if(p.eventClass()<2){

        int i( static_cast<int>(log(p.energy()/E6)/log(2.35)) );
        int level = i+6;
        if( level>maxlevel) level= maxlevel;
        if( level<minlevel) return 0;

        double diff = 1e9;
        astro::SkyDir sd(0,0);

        double alpha(0.);
        //figure out the associated source
        for(int it(0);it<m_skydir.size();++it) {
            double dot = p.difference(m_skydir[it]);
            if(dot<diff) {
                diff=dot;
                sd = m_skydir[it];
                alpha = m_alphas[it][i];
            }
        }
        if(alpha<0) alpha=0;
        double p0,p1;
        //From IRF parameters of 68% containment


        if(p.eventClass()==0) {
            p0 = 0.058;
            p1 = 0.000377;
        } else {
            p0 = 0.096;
            p1 = 0.0013;
        }  
        double sigmasq = (p0*p0*pow(p.energy()/100,-1.6))+p1*p1;
        double utest = diff*diff/sigmasq/pointlike::PointSourceLikelihood::sigma_level(level)/pointlike::PointSourceLikelihood::sigma_level(level)/2;
        //if scaled deviation is within the cone and enough signal photons
        if(utest<s_umax&&alpha>0.15) {
            added=1;
            HepRotation glast = p.Rot();
            //rotate skydirs into glast centric coordinates
            Hep3Vector meas = glast.inverse()*p.dir();
            Hep3Vector tru = glast.inverse()*sd();
            //accumulate likelihood statistics
            m_roti.acc(tru,meas,sigmasq,alpha,level);
            ++m_photons[level-minlevel];
        }
    }
    return added;
}

AlignProc::LikeSurface AlignProc::fitsurface() {
    //Least squares method described in Numerical Recipes 15.4 - General Least Squares fit
    std::vector<double> params;

    double norm = m_roti.likelihood(0,0,0);
    // Rows of A are [ xi**2  xi  yi**2  yi  zi**2  zi   1 ]
    int n = RotationInfo::points();
    TMatrixD A((2*n+1)*(2*n+1)*(2*n+1),10);
    TMatrixD b((2*n+1)*(2*n+1)*(2*n+1),1);
    for(int x=-n;x<=n;++x) {
        for(int y=-n;y<=n;++y) {
            for(int z=-n;z<=n;++z) {
                int row = (2*n+1)*(2*n+1)*(x+n)+(2*n+1)*(y+n)+(z+n);
                A[row][0]=m_arcsec*m_arcsec*x*x;
                A[row][1]=m_arcsec*x;
                A[row][2]=m_arcsec*m_arcsec*y*y;
                A[row][3]=m_arcsec*y;
                A[row][4]=m_arcsec*m_arcsec*z*z;
                A[row][5]=m_arcsec*z;
                A[row][6]=1;
                A[row][7]=m_arcsec*m_arcsec*x*y;
                A[row][8]=m_arcsec*m_arcsec*x*z;
                A[row][9]=m_arcsec*m_arcsec*y*z;
                double element = m_roti.likelihood(x,y,z)-norm;
                b[row][0]=element;
            }
        }
    }
#ifdef PRINT
    for(int z=-n;z<=n;++z) {
        for(int y=-n;y<=n;++y) {
            for(int x=-n;x<=n;++x) {
                //int row = (2*n+1)*(2*n+1)*(x+n)+(2*n+1)*(y+n)+(z+n);
                surfacexyz << m_roti.likelihood(x,y,z)-norm << "\t";
            }
            surfacexyz << std::endl;
        }
    }
#endif
    Double_t* err(0);
    //Least squares equation: (At A)*x = At*b, x = (At A)**(-1)*At*b
    TMatrixD pmatrix(10,(2*n+1)*(2*n+1)*(2*n+1));
    TMatrixD p2(10,10);
    pmatrix.Transpose(A);
    p2 = pmatrix*(A);
    p2 = p2.Invert(err);
#ifdef PRINT
    for(int i=0;i<=9;++i) {
        for(int j=0;j<=9;++j) {
            outfile << sqrt(p2[i][j]) << "\t";
        }
        outfile << std::endl;
    }
#endif
    if(err) std::cout << "BAD INVERSE - Covariance Matrix is singular" << std::endl;
    TMatrixD p3(10,1);
    p3 = p2*A.T()*b;
    for(int x=0;x<10;++x) {
        params.push_back(p3[x][0]);
    }
    return AlignProc::LikeSurface(params);
}

void AlignProc::addRot(double x, double y, double z) { 
    HepRotationX mx(x*M_PI/648000);
    HepRotationY my(y*M_PI/648000);
    HepRotationZ mz(z*M_PI/648000);
    s_hr = mx*my*mz;
}

std::vector<double> AlignProc::alignment() {
    LikeSurface l = fitsurface();
    TMatrixD A(3,3);
    TMatrixD b(3,1);
    //Least squares solution to the second order likelihood surface maxmimum
    A[0][0]=2*l(0);
    A[0][1]=l(7);
    A[0][2]=l(8);
    A[1][0]=A[0][1];
    A[1][1]=2*l(2);
    A[1][2]=l(9);
    A[2][0]=A[0][2];
    A[2][1]=A[1][2];
    A[2][2]=2*l(4);
    b[0][0]=-l(1);
    b[1][0]=-l(3);
    b[2][0]=-l(5);
    Double_t* err(0);
    A = A.Invert(err);
    b = A*b;
    if(err) std::cout << "BAD INVERSE - Matrix is nearly singular" << std::endl;
    std::vector<double> alignment;
    for(int i=0;i<3;++i) {
        alignment.push_back(b[i][0]);
    }
    return alignment;
}

void  AlignProc::set_rot(const CLHEP::HepRotation& R)
{
    s_hr=R;
}
