/** 
Data Processing file, operates on a given Photon

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/AlignProc.cxx,v 1.2 2007/11/21 16:36:21 burnett Exp $

*/

#include "pointlike/AlignProc.h"
#include "CLHEP/Matrix/Matrix.h"
#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TSystem.h"

#include <float.h>
#include <stdexcept>

using namespace pointlike;
using namespace CLHEP;
//#define PRINT

namespace{

    //extra factor (kluge for now) for correcting scaled deviation
    double minE = 555;

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
        "PtDecx"//, "CTBClassLevel"
    };

    double scale[] = {1.,1.86,1.,1.};
}

//setup initally with no rotation
HepRotation AlignProc::s_hr(0,0,0);
int AlignProc::s_umax = 50;
double AlignProc::s_scale = 0.46;

AlignProc::AlignProc(std::vector<astro::SkyDir>& sources,std::vector<std::string>& files,double arcmins, double offx, double offy, double offz, int start, int stop): 
m_photons(0),
m_start(start),
m_stop(stop),
m_skydir(sources),
m_arcmin(arcmins),
m_roti(arcmins*M_PI/10800,offx*M_PI/10800,offy*M_PI/10800,offz*M_PI/10800)
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
    tt->GetEvent(0);
    int starttime = static_cast<int>(tt->GetLeaf("EvtElapsedTime")->GetValue());
    bool flag(true);
    //for each entry  
    for(int i(0);i<entries&&(flag||m_start==-1);++i) {

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
        pointlike::AlignProc::Photona p = event(row);
        if(row[3]-starttime>=m_start) {
            add(p);
            ShowPercent(i,entries,i);
        }
        if(row[3]-starttime>m_stop) {
            flag=false;
        }
        row.clear();
    }
}

void AlignProc::loadfits(const std::string& /*file */) {
    //todo : process fits file in same manner as ROOT
}

pointlike::AlignProc::Photona AlignProc::event(std::vector<float>& row) {
    float ra(0), dec(0), energy(0); // photon info
    float raz(0), decz(0), rax(90), decx(0); // sc orientation: default orthogonal
    double time(0);
    int event_class(99);
    int flag =1;
    for(unsigned int i = 0;i<row.size();++i) {
        if(row[i]<-1e7) flag=0;
    }
    if(flag) {
        time = row[3];
        event_class = static_cast<int>(row[4]);
        event_class = event_class>4? 0 : 1;  // front/back map to event class 0/1
        energy = row[2];
        if(event_class<99){
            ra = row[0];
            dec = row[1];
            raz = row[5];
            decz = row[6];
            rax = row[7];
            decx = row[8];
            energy/=scale[event_class];
        }
    }
    return pointlike::AlignProc::Photona(astro::SkyDir(ra, dec), energy, time, event_class ,
        astro::SkyDir(raz,decz),astro::SkyDir(rax,decx));
}

int AlignProc::add(pointlike::AlignProc::Photona& p){
    int added(0);
    //unused int cl = p.eventClass();
    //above healpix level 8 (mostly signal photons)
    if(p.energy()>=minE&&p.eventClass()<2){
        double diff = 1e9;
        astro::SkyDir sd(0,0);
        //figure out the associated source
        for(std::vector<astro::SkyDir>::const_iterator it = m_skydir.begin();it!=m_skydir.end();++it) {
            double dot = p.difference(*it);
            if(dot<diff) {
                diff=dot;
                sd = *it;
            }
        }
        double p0,p1;
        //From IRF parameters of 68% containment

        int i( static_cast<int>(log(p.energy()/42)/log(2.35)) );
        if( i>9-1) i= 9-1;
        int level = i+5;

        if(p.eventClass()==0) {
            p0 = 0.058;
            p1 = 0.000377;
        } else {
            p0 = 0.096;
            p1 = 0.0013;
        }  
        double sigmasq = (p0*p0*pow(p.energy()/100,-1.6))+p1*p1;
        double utest = diff*diff/sigmasq/pointlike::PointSourceLikelihood::sigma_level[level]/pointlike::PointSourceLikelihood::sigma_level[level]/2;
        //if scaled deviation is within the cone
        if(utest<s_umax) {
            added=1;
            HepRotation glast = p.Rot();
            //rotate skydirs into glast centric coordinates
            Hep3Vector meas = glast.inverse()*p.dir();
            Hep3Vector tru = glast.inverse()*sd();
            //accumulate likelihood statistics
            m_roti.acc(tru,s_hr.inverse()*meas,sigmasq,level);
            ++m_photons;
        }
    }
    return static_cast<int>(p.time());
}

std::vector<double> AlignProc::fitparam() {
    //Least squares method described in Numerical Recipes 15.4 - General Least Squares fit
    std::vector<double> params;
    HepMatrix A(125,10,0);
    HepMatrix b(125,1,0);
    // Rows of A are [ xi**2  xi  yi**2  yi  zi**2  zi   1 ] 
    for(int x=-2;x<=2;++x) {
        for(int y=-2;y<=2;++y) {
            for(int z=-2;z<=2;++z) {
                int row = 25*(x+2)+5*(y+2)+z+2;
                A[row][0]=m_arcmin*m_arcmin*x*x;
                A[row][1]=m_arcmin*x;
                A[row][2]=m_arcmin*m_arcmin*y*y;
                A[row][3]=m_arcmin*y;
                A[row][4]=m_arcmin*m_arcmin*z*z;
                A[row][5]=m_arcmin*z;
                A[row][6]=1;
                A[row][7]=m_arcmin*m_arcmin*x*y;
                A[row][8]=m_arcmin*m_arcmin*x*z;
                A[row][9]=m_arcmin*m_arcmin*y*z;
                b[row][0]=m_roti.likelihood(x,y,z);
            }
        }
    }
    int err;
#ifdef PRINT
    std::cout << A << std::endl;
#endif
    //Least squares equation: (At A)*x = At*b, x = (At A)**(-1)*At*b
    HepMatrix pmatrix = (A.T()*A).inverse(err);
    if(err) std::cout << "BAD INVERSE - Covariance Matrix is singular" << std::endl;
    pmatrix = pmatrix*A.T()*b;
    for(int x=0;x<10;++x) {
        params.push_back(pmatrix[x][0]);
    }
    return params;
}

