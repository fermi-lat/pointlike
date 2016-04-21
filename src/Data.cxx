/** @file Data.cxx
@brief implementation of Data

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Data.cxx,v 1.74 2012/04/22 23:50:17 wallacee Exp $

*/


#include "pointlike/Data.h"
#include "pointlike/Alignment.h"
#include "EventList.h"

#include "astro/SkyDir.h"
#include "astro/GPS.h"
#include "astro/PointingTransform.h"
#include "astro/PointingHistory.h"
#include "astro/SolarSystem.h"

#include "skymaps/SkyImage.h"
#include "skymaps/BinnedPhotonData.h"
#include "skymaps/IParams.h"
#include "skymaps/Gti.h"
#include "skymaps/Photon.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"


// --ROOT --
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TSystem.h"


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
using skymaps::Gti;

using namespace pointlike;

std::string Data::s_ft2file = std::string("");
astro::PointingHistory* Data::s_history(0);
namespace {
    // file-scope pointer to gps instance
    astro::GPS* gps (astro::GPS::instance()); 

    int use_earth(0);  // flag to convert to Earth coordinates

}
void Data::useEarthCoordinates()
{
    use_earth=1;
}

//! @brief define FT2 file to use for rotation
void Data::setHistoryFile(const std::string& history)
{
    s_ft2file = history;
    if( s_history==0 ){
        s_history = new astro::PointingHistory(history);
    }else{
        s_history->readFitsData(history);
    }
    // temporary put same file into GPS
#if 0 // this screws up multiple files: do we need it?
    gps->setPointingHistoryFile(history);
#endif 
}

int Data::s_class_level=3; // diffuse selection
int Data::class_level(){return s_class_level;}
void Data::set_class_level(int level){s_class_level=level;}

double Data::s_zenith_angle_cut(105.); // standard cut
double Data::zenith_angle_cut(){ return s_zenith_angle_cut;}
void   Data::set_zenith_angle_cut(double cut){s_zenith_angle_cut=cut;}

double Data::s_theta_min_cut(0.);
double Data::theta_min_cut() {return s_theta_min_cut;}
void   Data::set_theta_min_cut(double cut){s_theta_min_cut=cut;}

double Data::s_theta_cut( 66.4 ); // standard cut (acos(0.4))
double Data::theta_cut(){ return s_theta_cut;}
void   Data::set_theta_cut(double cut){s_theta_cut=cut;}

bool Data::s_use_mc_energy=false; //use measured energy by default
bool Data::use_mc_energy(){return s_use_mc_energy;}
void Data::set_use_mc_energy(bool use_it){s_use_mc_energy = use_it;}

skymaps::Gti Data::s_gti_mask;   // default constructor is empty
void Data::set_Gti_mask(skymaps::Gti gti){
    s_gti_mask = gti;
}

// default alignment object
pointlike::Alignment* Data::s_alignment=new Alignment();

void Data::set_alignment(const std::string & filename)
{
    delete s_alignment;
    s_alignment = new Alignment(filename);
}

#ifdef WIN32
#include <float.h> // used to check for NaN
#else
#include <cmath>
#endif

namespace {

    std::string root_names[] = {"FT1Ra", "FT1Dec", 
        "CTBBestEnergy", "EvtElapsedTime",
        "FT1ConvLayer","PtRaz", 
        "PtDecz","PtRax",
        "PtDecx", "CTBClassLevel","FT1ZenithTheta","FT1Theta","CTBCORE",
        "CTBBestEnergyProb","CTBBestEnergyRatio","FT1Phi","GltTower","PtSCzenith","PtPos[3]"//THB, "McSourceId"//,"FT1ZenithTheta","CTBBestEnergyRatio","CTBCORE","CTBGAM"
    };

    std::string root_names2[] = {"FT1Ra", "FT1Dec", 
        "CTBBestEnergy", "EvtElapsedTime",
        "FT1ConvLayer","PtRaz", 
        "PtDecz","PtRax",
        "PtDecx", "CTBClassLevel","FT1ZenithTheta","FT1Theta","CTBCORE",
        "CTBBestEnergyProb","CTBBestEnergyRatio","FT1Phi","GltTower","PtSCzenith","PtPos"
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

    void ShowPercent(int sofar, int total, int found, double speed)
    {
        static int toskip(50), skipped(0);
        if(++skipped<toskip) return; skipped=0;
        static int lastpercent(-1);
        int percent( static_cast<int>(100 * sofar / total +0.5) );
        if( percent==lastpercent) return;
        lastpercent=percent;
        char   s[70];
        if ((total-sofar)/speed<60) {
            sprintf(s, "%d%%, %d found. %d events/sec about %d sec left", percent, found, (int)(speed), (int)((total-sofar)/speed));
        } else {
            sprintf(s, "%d%%, %d found. %d events/sec about %d:%d min left", percent, found, (int)(speed), (int)((total-sofar)/(speed*60)), (int)((total-sofar)/(speed))%60);
        }
        std::cout << s;
        if (sofar < total)
        {
            for (unsigned int j = 0; j < strlen(s); ++j)
                std::cout << "\b";
        }
        else
            std::cout << std::endl;
    }

    // default binner to use
    skymaps::PhotonBinner* binner( new skymaps::PhotonBinner(4) );

} // anon namespace

Data::Data(const embed_python::Module& setup)
: m_log(0)
{
    std::string pixelfile("")
        , output_pixelfile("");

    static std::string prefix("Data.");
    setup.getValue(prefix+"pixelfile", pixelfile, "");

    if(!pixelfile.empty()){
        try {
            m_data = new BinnedPhotonData(pixelfile, "BANDS");
        } catch( const std::exception& ){
            m_data = new BinnedPhotonData(pixelfile, "PHOTONMAP");
        }
        return;
    }

    int  event_class, source_id, class_level;
    std::vector<std::string> filelist;
    std::string history_filename;
    std::vector<double>alignment;

    // special flag to convert to Earth coordinates
    setup.getValue(prefix+"use_earth", use_earth, false);

    setup.getList(prefix+"files", filelist);
    setup.getValue(prefix+"event_class", event_class);
    setup.getValue(prefix+"source_id",  source_id);
    setup.getValue(prefix+"start_time", m_start, -1);
    setup.getValue(prefix+"stop_time" , m_stop, -1);
    setup.getValue(prefix+"output_pixelfile", output_pixelfile, "");
    setup.getValue(prefix+"history",   history_filename, "");
    setup.getValue(prefix+"class_level", class_level,3);
    set_class_level(class_level);
    if( !history_filename.empty()){
        log() << "loading history file: " << history_filename << std::endl;
        Data::setHistoryFile(history_filename);
    }
    setup.getList(prefix+"LATalignment", alignment);
    if( alignment.size()>0) {
        Data::set_rot(alignment);
    }
    std::string alignment_data;
    setup.getValue(prefix+"alignment_data", alignment_data, "");
    if( !alignment_data.empty() ){
        delete s_alignment;
        s_alignment = new Alignment(alignment_data);
    }

    double bins_per_decade(0); 
    setup.getValue(prefix+"bins_per_decade", bins_per_decade, bins_per_decade);

    std::vector<double>energy_bins;
    setup.getList(prefix+"energy_bins", energy_bins);

    if( energy_bins.size()> 0 ){
        m_data = new BinnedPhotonData(*(new skymaps::PhotonBinner(energy_bins)));
    }else if (bins_per_decade>0) {
        m_data = new BinnedPhotonData(bins_per_decade);
    } else { 
        m_data = new BinnedPhotonData(*binner);
    };

    if( filelist.empty()){
        throw std::invalid_argument("Data: no data specified: expected either a pixel file or a list of data files");
    }
    load_filelist(filelist, event_class, source_id);

    if( !output_pixelfile.empty() ) {
        log() << "writing output pixel file :" << output_pixelfile ;
        m_data->write(output_pixelfile);
        log() << " with GTI range " << int(m_data->gti().minValue())<<"-"<<int(m_data->gti().maxValue())
            << std::endl;
    }
}

void Data::load_filelist(const std::vector<std::string>& filelist, int event_class, int source_id)
{
    log() << "Reading from FT1 file list. " << std::endl;
    if( m_start>0 ) {
        log() << "\tselecting time range: " << static_cast<int>(m_start) << " to " << 
            static_cast<int>(m_stop) << std::endl;
    }
    log() << "\tCuts:  zenith theta and instrument theta (min,max): " << s_zenith_angle_cut << ", " 
        << s_theta_min_cut << ", " << s_theta_cut<< std::endl;
    if( use_earth) { 
        log() << "\tConverting to zenith coordinates, not applying zenith cut" << std::endl;
    }
    for( std::vector<std::string>::const_iterator it = filelist.begin(); 
        it !=filelist.end(); ++it)
    {
        const std::string& inputFile(*it);
        log() << "\t" << inputFile <<": " ;
        if (addgti(inputFile)){
            add(inputFile, event_class, source_id);
        }
    }

}


void Data::add(const std::string& inputFile, int event_type, int source_id)
{

    int photoncount(m_data->photonCount()), pixelcount(m_data->pixelCount());

    if( inputFile.find(".root") != std::string::npos) {
        lroot(inputFile, event_type);
        std::cout 
            << "photons kept: "  << (m_data->photonCount() -photoncount) << " (total: " << m_data->photonCount() <<") "
            << std::endl;


    }else {

        EventList photons(inputFile, source_id>-1, s_use_mc_energy);
        AddPhoton adder(*m_data, event_type, m_start, m_stop, source_id, s_gti_mask, photons.data_pass());
        //AddPhoton adder(*m_data, event_type, m_start, m_stop, source_id, s_gti_mask, photons.pass7());
        //AddPhoton adder(*m_data, event_type, m_start, m_stop, source_id, m_data->gti(), photons.pass7());

        AddPhoton added =std::for_each(photons.begin(), photons.end(), adder );
        log() 
            << "\t\t photons found: " << added.found() <<", kept: "  << (m_data->photonCount() -photoncount) 
            << " (total: " << m_data->photonCount() <<") "
            << std::endl;
    }

}

bool Data::addgti(const std::string& inputFile)
{
    bool ok = true;
    try
    {

        if(inputFile.find(".root") == std::string::npos) {
            // a FITS file: check for overlap
            skymaps::Gti tnew(inputFile);
            // apply gti mask
            if( s_gti_mask.getNumIntervals()>0){
                tnew = tnew & s_gti_mask;  // intersection
            }
            double tmin(tnew.minValue()), tmax(tnew.maxValue());
            ok = (m_start==0 || tmax > m_start)
                && (m_stop==0 || tmin < m_stop ) ;
            if( ok ) {
                double max(m_stop==0 ? tmax : std::min(m_stop,tmax));
                skymaps::Gti timerange(tnew.applyTimeRangeCut(m_start, max));
                m_data->addgti(timerange); 
                log() << " found interval " 
                    << int(tmin)<<"-"<< int(tmax)<<  std::endl;
            }else{
                log() <<  " not in selected range" << std::endl;
            }

        }else {

            // extract begin and end times from the ROOT file
            TFile *tf = new TFile(inputFile.c_str(),"READONLY");
            TTree *tt = static_cast<TTree*>(tf->Get("MeritTuple"));
            TLeaf *tl = tt->GetLeaf(root_names[3].c_str());
            tt->SetBranchStatus("*", 0); // turn off all branches
            // but the EvtElapsedTime
            tt->SetBranchStatus(root_names[3].c_str(), 1);
            int entries (static_cast<int>(tt->GetEntries()) );
            tt->GetEvent(0);
            double start( tl->GetValue());
            tt->GetEvent(entries-1);
            double stop( tl->GetValue() );
            delete tf;

            skymaps::Gti tnew;
            tnew.insertInterval( m_start>start? m_start:start, 
                m_stop<stop&&m_stop>0? m_stop:stop);
            m_data->addgti(tnew);
            std::cout << " found interval " 
                << int(tnew.minValue())<<"-"<< int(tnew.maxValue())
                << ", total: " << (m_data->gti().computeOntime()<0?0:m_data->gti().computeOntime()) << " s." 
                <<  std::endl;
        }


    }
    catch(const std::exception& e)
    {
        std::cerr << "\nCaught exception " << typeid(e).name() 
            << " \"" << e.what() << "\"" << std::endl;
        std::cerr << "Unable to access gti info from file " << inputFile << std::endl;
        ok = false;
    }
    return ok;
}
Data::Data(const std::string& inputFile, int event_type, double tstart, double tstop, int source_id, std::string irf)
: m_data(new BinnedPhotonData(*binner))
, m_start(tstart), m_stop(tstop)
, m_log(0)
{
    if(!irf.empty()) skymaps::IParams::init(irf);
    if( addgti(inputFile)){
    	add(inputFile, event_type, source_id);
    }
}

Data::Data(std::vector<std::string> inputfiles, int event_type, double tstart, double tstop, int source_id, std::string ft2file, std::string irf)
: m_data(new BinnedPhotonData(*binner))
, m_start(tstart), m_stop(tstop)
, m_log(0)
{
    if( !ft2file.empty() ) setHistoryFile(ft2file); // this is actually a global
    if( !irf.empty() ) skymaps::IParams::init(irf);

    load_filelist(inputfiles,event_type,source_id);
}

Data::Data(const std::string & inputFile, const std::string & tablename)
: m_data(new BinnedPhotonData(inputFile, tablename))
, m_start(0), m_stop(0)
, m_log(0)
{
    m_data->updateIrfs();
    //not needed? addgti(inputFile);
}

Data::~Data()
{
    delete m_data;
}

//ROOT event extraction, only used by Data::lroot
skymaps::Photon events(std::vector<double>& row) {

    // extract stuff from row, assuming order in root_names
    float ra(row[0]), dec(row[1]); 
    double energy(row[2])
        ,  time( row[3]) ;

    int event_class( row[4]>4? 0 : 1 );

    // for transformation
    float raz(row[5]), decz(row[6])
        , rax(row[7]), decx(row[8]); 

    // these not actually relevant here
    int class_level( static_cast<int>(row[9]) );
    double zenith_angle(row[10]); 
    int source_id(-1); // not set now?

    // create the local Photon object from these data, have it return a transformed skymaps::Photon
    Photon p(astro::SkyDir(ra, dec), energy, time, event_class , 
        source_id, astro::SkyDir(raz,decz),astro::SkyDir(rax,decx),
        zenith_angle,class_level);
#if 0
    return p.transform(Data::get_rot().inverse());
#else
    return p.transform(Data::get_rot(time));
#endif
}

void Data::lroot(const std::string& inputFile, int event_class) {
    TFile *tf = new TFile(inputFile.c_str(),"READ");
    TTree *tt = static_cast<TTree*>(tf->Get("MeritTuple"));
    double time(0);
    tt->SetBranchStatus("*", 0); // turn off all branches
    //turn on appropriate branches
    for(unsigned int j(0); j< sizeof(root_names)/sizeof(std::string); j++){
        tt->SetBranchStatus(root_names[j].c_str(), 1);
    }
    int entries = static_cast<int>(tt->GetEntries());
    int bpt=entries/2;
    int step=entries/4;
    int bstart=0;
    int bstop=entries;
    if(m_start>0) {
        for(;step>2;step/=2) {
            tt->GetEvent(bpt);
            double ctime = tt->GetLeaf(root_names[3].c_str())->GetValue();
            bpt+=(m_start>ctime?step:-step);
        }
        bstart = bpt;
        bpt=entries/2;
        step=entries/4;
        for(;step>2;step/=2) {
            tt->GetEvent(bpt);
            double ctime = tt->GetLeaf(root_names[3].c_str())->GetValue();
            bpt+=(m_stop>ctime?step:-step);
        }
        bstop=bpt;
    }
    tt->GetEvent(entries-1);
    int endtime = static_cast<int>(tt->GetLeaf("EvtElapsedTime")->GetValue());
    if(endtime<m_start&&m_start>0) return;
    std::vector<double> row;
    bool flag(true);
    //for each entry
    double stime(clock()/CLOCKS_PER_SEC);
    for(int i(bstart);i<bstop&&(flag||m_start<0);++i) {
        tt->GetEvent(i);
        //for each
        for( size_t j(0); j< sizeof(root_names)/sizeof(std::string); j++){
            const std::string name = root_names[j].c_str();
            TLeaf * tl = tt->GetLeaf(root_names[j].c_str());
            if(0==tl) {
                tl = tt->GetLeaf(("_" + root_names[j]).c_str());
                if(0==tl) {
                    //tt->Print();
                    throw std::invalid_argument(std::string("Tuple: could not find leaf ")+root_names[j]);
                }
            }
            float v(0);
                v = tl->GetValue();
                row.push_back(isFinite(v)?v:-1e8);
        }
        double time(row[3])
            ,thetazenith(row[10]), ctbclasslevel(row[9]),theta(row[11]),core(row[12]),eprob(row[13]),ratio(row[14]) ;
        // perform event selection -- note uses local function events to transform
        double curtime(clock()/CLOCKS_PER_SEC);
        double speed((i-bstart)/(curtime-stime));
        if( 
            thetazenith < Data::zenith_angle_cut() &&
            theta < Data::theta_cut() && theta>= Data::theta_min_cut()

            && ctbclasslevel>= Data::class_level() 
            && (Data::minTime()==0 || time>Data::minTime() )
            && (Data::maxTime()==0 || time<Data::maxTime() ) ) 
        {

            //set_rot(rx[tower],ry[tower],rz[tower]);
            if(i-bstart) ShowPercent(i-bstart+1,bstop-bstart+1,i-bstart+1,speed);
            skymaps::Photon p =  events(row); 
            if ((p.eventClass()==event_class||event_class==-1) 
                ) {
                    m_data->addPhoton(p);

            }
        }
        if(endtime<Data::minTime()||time>Data::maxTime())  flag = false;
        row.clear();
    }
    delete tf;
}

void Data::set_rot(double arcsecx,double arcsecy,double arcsecz) {
    delete s_alignment;
    s_alignment = new Alignment(arcsecx, arcsecy, arcsecz);
}
void Data::set_rot(std::vector<double> align) {
    assert( align.size()==3); // otherwise bad interface setup
    set_rot(align[0], align[1], align[2]);
}

const CLHEP::HepRotation& Data::get_rot(double time) {
    return s_alignment->rotation(time);
}

const astro::PointingInfo& Data::get_pointing(double time) {
    return (*s_history)(time);
}


const std::string& Data::historyfile() {
    return s_ft2file;
}

bool Data::inTimeRange(double time) {
    return time>(*s_history).startTime()&&time<(*s_history).endTime();
}

void Data::info(std::ostream& out)
{
    m_data->info(out);
}


///@brief change default binning: must be done before loading data files
void Data::setEnergyBins(const std::vector<double>& bins)
{
    delete binner;
    binner = new skymaps::PhotonBinner(bins);
}

///@brief change default binner: must be done before loading data files
void Data::setPhotonBinner(skymaps::PhotonBinner* b)
{
    //THB confused by SWIG?   delete binner;
    binner = b;
}


void Data::combine_bands()
{
    using skymaps::BinnedPhotonData;
    using skymaps::Band;
    int nside(0);
    double sigma(0);
    Band* prev (0);
    for( BinnedPhotonData::iterator it(m_data->begin()); it!=m_data->end(); ++it){
        Band& current = *it;
        if( current.nside() == nside && current.sigma() == sigma){
            // combine with previous and remove it
            current.add(*prev);
            m_data->remove(*prev);
            nside=0;
        }else {
            prev = &current;
            nside = current.nside();
            sigma = current.sigma();
        }
    }

}

void Data::write(const std::string & outputFile, bool clobber ) const
{
    m_data->write(outputFile, clobber);
}
