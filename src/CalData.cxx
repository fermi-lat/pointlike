/** @file CalData.cxx
@brief implementation of CalData

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/CalData.cxx,v 1.7 2008/11/11 01:31:16 mar0 Exp $

*/


#include "pointlike/CalData.h"
#include "skymaps/BinnedPhotonData.h"

// --ROOT --
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TSystem.h"
#include "TRandom.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include <ctime>


using astro::SkyDir;
using skymaps::BinnedPhotonData;

using namespace pointlike;


int CalData::s_class_level=3; 
int CalData::class_level(){return s_class_level;}

#ifdef WIN32
#include <float.h> // used to check for NaN
#else
#include <cmath>
#endif

namespace {

    std::string root_names[] = {"CTBBestEnergy",
        "Tkr1FirstLayer"/*"FT1ConvLayer"*/, "CTBClassLevel" ,"McDirErr","McEnergy","CTBCORE","CTBBestEnergyProb","CTBBestEnergyRatio"//,"FT1Theta"
    };

    bool isFinite(double val) {
        using namespace std; // should allow either std::isfinite or ::isfinite
#ifdef WIN32 
        return (_finite(val)!=0);  // Win32 call available in float.h
#else
        return (isfinite(val)!=0); // gcc call available in math.h 
#endif
    }

    inline static void do_load (void)
    {
        static bool first = true;
        if( first) {
            gSystem->Load("libTree");
            first=false;
        }
    }

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

    time_t t = time(NULL);
    int tseed = t%65536;

    TRandom *rand = new TRandom(tseed);

        // default binner to use
    skymaps::PhotonBinner* binner(double emin,double emax) {
        return new skymaps::PhotonBinner(4);
    }

    skymaps::PhotonBinner pb = *binner(0,0);

    //ROOT event extraction
	//Changed to use skymaps::Photon. Shouldn't change any behavior. -EEW
    skymaps::Photon events(std::vector<float>& row,int event_type, double emin, double emax) {
        float ra(0), dec(0), energy(0); // photon info

        double time(0);
        int event_class(99);
        int source_id(0);
        int class_level(0);
        int flag =1;
        for(unsigned int i = 0;i<row.size();++i) {
            if(row[i]<-1e7) {
                flag=0;
                std::cerr << "Bad data: time="<< std::fixed<< row[3]<< ", index, value: " << i << ", " << row[i]<< std::endl;
            }
        }
        if(flag) {
            event_class = static_cast<int>(row[1]);
            event_class = event_class>5? 0 : 1;  // front/back map to event class 0/1
            energy = row[0];
            //unused double mcenergy = row[4];
            class_level = static_cast<int>(row[2]);
            // NB. Selecting wired-in class definition (transient, source, diffuse
            
            if( (class_level < CalData::class_level())||(energy<emin)||(energy>emax)||(event_class!=event_type)) event_class=99;
            else{
                double r = 2*rand->Uniform()*M_PI;
                double diff = 180/M_PI*row[3];
                dec = diff*sin(r);
                ra = diff*cos(r)/cos(dec*M_PI/180);
                skymaps::Band bb = pb(skymaps::Photon(astro::SkyDir(0,0),energy,0));
                double ratio = (bb.emin()/energy);
                double rseed = rand->Uniform();
                if(ratio<rseed) event_class=99;
                if(row[5]<=0.35) event_class=99;
                if(row[6]<=0.1) event_class=99;
                if(row[7]>=5) event_class=99;
            }

        }
        return skymaps::Photon(astro::SkyDir(ra,dec),energy,time,event_class, source_id);
    }


} // anon namespace

CalData::CalData(const std::string& inputFile, int event_type, double emin, double emax, int /*source_id*/)
: m_data(new BinnedPhotonData(*binner(emin,emax)))
{
    lroot(inputFile,event_type,emin,emax);
}


CalData::~CalData()
{
    delete m_data;
}

void CalData::lroot(const std::string& inputFile,int event_type,double emin,double emax) {
    TFile *tf = new TFile(inputFile.c_str(),"READ");
    TTree *tt = static_cast<TTree*>(tf->Get("MeritTuple"));
    tt->SetBranchStatus("*", 0); // turn off all branches
    //turn on appropriate branches
    for(unsigned int j(0); j< sizeof(root_names)/sizeof(std::string); j++){
        tt->SetBranchStatus(root_names[j].c_str(), 1);
    }
    int entries = static_cast<int>(tt->GetEntries());
    std::vector<float> row;
    tt->GetEvent(0);
    //int starttime = static_cast<int>(tt->GetLeaf("EvtElapsedTime")->GetValue());
    //unused bool flag(true);
    //for each entry  
    for(int i(0);i<entries;++i) {
        tt->GetEvent(i);
        //for each
        for( size_t j(0); j< sizeof(root_names)/sizeof(std::string); j++){
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
        skymaps::Photon p = events(row,event_type,emin,emax);
#if 0 // does not compile with gcc
        p.eventClass()<99?m_data->addPhoton(p):0;
#else
        if( p.eventClass()<99 ) m_data->addPhoton(p);
#endif
        row.clear();
        //if(m_data->size()>0) {
        ShowPercent(i,entries,i);
        //if(m_data->photonCount()>10000) break;
        //}
    }
    delete tf; 
}
