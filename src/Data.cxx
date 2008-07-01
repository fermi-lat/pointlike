/** @file Data.cxx
@brief implementation of Data

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/Data.cxx,v 1.42 2008/06/30 23:38:22 burnett Exp $

*/


#include "pointlike/Data.h"

#include "astro/SkyDir.h"
#include "astro/Photon.h"
#include "astro/PointingTransform.h"
#include "astro/PointingHistory.h"

#include "skymaps/SkyImage.h"
#include "skymaps/BinnedPhotonData.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include "skymaps/Gti.h"

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


using astro::SkyDir;
using skymaps::BinnedPhotonData;
using skymaps::Gti;

using namespace pointlike;

std::string Data::s_ft2file = std::string("");
astro::PointingHistory* Data::s_history;

int Data::s_class_level=2; 
int Data::class_level(){return s_class_level;}

double Data::s_zenith_angle_cut(100.);
double Data::zenith_angle_cut(){
    return s_zenith_angle_cut;
}

CLHEP::HepRotation Data::s_rot = CLHEP::HepRotationX(0)*CLHEP::HepRotationY(0)*CLHEP::HepRotationZ(0);

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
        "PtDecx", "CTBClassLevel"//THB, "McSourceId"//,"FT1ZenithTheta","CTBBestEnergyRatio","CTBCORE","CTBGAM"
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



    /** @class Photon
    @brief derive from astro::Photon to allow transformation

    */

    class Photon : public astro::Photon {
    public:
        /// ctor sets photon, also rotation info.
        Photon(const astro::SkyDir& dir, double energy, 
            double time, int event_class, int source
            , const astro::SkyDir&scz, const astro::SkyDir& scx
            , double zenith_angle
            , int ctbclasslevel=1  // note default 
            )
            : astro::Photon(dir, energy, time, event_class, source)
            , m_zenith_angle(zenith_angle)
            , m_ctbclasslevel(ctbclasslevel)
        {
            Hep3Vector scy (scz().cross(scx()));
            m_rot = HepRotation(scx(), scy, scz());
        }
        /// make transformation in GLAST frame
        /// @param corr matrix that corrects direction in the GLAST frame
        astro::Photon transform(const HepRotation& corr)const
        {
            Hep3Vector 
                local( m_rot.inverse()*dir()),
                transformed( m_rot * corr * local );

            int evtclass(eventClass());
            if( evtclass<0 || evtclass>1){ // just a check: should be 0 or 1 for front/back
                //std::stringstream err;
                std::cerr << "Data::Photon::transform--invalid eventclass value, expect 0 or 1, got: "<< evtclass<< std::endl;
                //throw std::runtime_error(err.str());
            }
            double   emeas(energy());                    // measured energy
            return astro::Photon(SkyDir(transformed), emeas,time(),evtclass, source());
        }

        double zenith_angle()const{return m_zenith_angle;}
        int class_level()const{return m_ctbclasslevel;}
    private:
        HepRotation m_rot;
        double m_zenith_angle;
        int m_ctbclasslevel;

    };

    class AddPhoton: public std::unary_function<astro::Photon, void> {
    public:
        AddPhoton (BinnedPhotonData& map, int select, double start, double stop, int source )
            : m_map(map), m_select(select), m_start(start), m_stop(stop), m_source(source)
        {}
        void operator()(const Photon& gamma)
        {
            int event_class = gamma.eventClass();
            int sourceid = gamma.source();

            if( m_select>-1 && event_class!= m_select) return;
            if( m_start>0   && gamma.time()<m_start ||  m_stop>m_start && gamma.time()>m_stop) return;
            if( m_source>-1 && sourceid != m_source)return;

            // zenith angle cut
            if( ! isFinite(gamma.zenith_angle()) ) return; // catch NaN?
            if( gamma.zenith_angle()> Data::zenith_angle_cut()) return;
            int class_level( gamma.class_level() );
            if( class_level<2 ) return; // force source for now
            { 

                double energy(gamma.energy());

                astro::Photon gcopy(gamma.dir(), energy, gamma.time(), event_class, sourceid); 
                m_map.addPhoton(gcopy);
            }
        }
        BinnedPhotonData& m_map;
        int m_select;
        double m_start, m_stop;
        int m_source;
    };
    /**
    @class EventList
    @brief Manage a list of individual photons, defined by direction and energy, by 
    adapating a file 

    */
    class EventList
    {
    public:

        /** @brief ctor sets up container
        @param infile name of the input FT1 or ROOT file
        @param selectid True if will be selecting on eventid
        @param table_name must be "EVENTS" for FT1, or "MeritTuple"
        */
        EventList( std::string infile, bool selectid=false,
            std::string table_name="EVENTS");

        ~EventList();

        // make it a container by implementing a forward iterator
        class Iterator {
        public:
            Iterator(tip::Table::ConstIterator it, bool fits, bool selectid=false)
                : m_it(it)
                , m_fits(fits)
                , m_selectid(selectid)
            {}
            Photon operator*()const;             ///< dereference
            tip::Table::ConstIterator operator++(){return ++m_it;} ///< increment operator
            bool operator!=(const Iterator& other)const{return other.m_it!=m_it;}
        private:
            tip::Table::ConstIterator m_it;
            bool m_fits, m_selectid;
        };

        /// return iterator to access 
        Iterator begin();
        Iterator end();

        static std::string root_names[];


    private:
        tip::Table::ConstIterator m_itbegin, m_itend;
        const tip::Table * m_table;
        bool m_fits;
        bool m_selectid;
    };

    EventList::EventList(const std::string infile, bool selectid,
        std::string table_name)
        : m_fits(true)
        , m_selectid(selectid)
    {
        if( infile.find(".root") != std::string::npos) {
            table_name = "MeritTuple"; 
            m_fits=false;
        }
        // connect to input data
        const tip::Table * m_table = tip::IFileSvc::instance().readTable(infile, table_name, "");

        // save the iterators
        m_itbegin= m_table->begin();
        m_itend = m_table->end();

    }
    EventList::~EventList()
    {
        // seems to create crash
        // delete m_table;
    }

    Photon EventList::Iterator::operator*()const
    {
        float ra, dec, energy; // photon info
        float raz(0), decz(0), rax(90), decx(0); // sc orientation: default orthogonal
        double time;
        double zenith_angle;
        int event_class, ctbclasslevel(1);
        int source(-1);

        // FT1 names
        static std::string fits_names[]={"RA", "DEC", 
            "ENERGY", "TIME", "EVENT_CLASS", "ZENITH_ANGLE", "CTBCLASSLEVEL", "MC_SRC_ID"};
        // corresponging names in the ROOT MeritTuple
        static std::string root_names[]={"FT1Ra", "FT1Dec", 
            "CTBBestEnergy", "EvtElapsedTime"
            , "FT1ConvLayer", "FT1ZenithAngle", "McSourceId"};

        std::string *names = m_fits?  fits_names : root_names;

        (*m_it)[*names++].get(ra);
        (*m_it)[*names++].get(dec);
        (*m_it)[*names++].get(energy);
        (*m_it)[*names++].get(time);
        (*m_it)[*names++].get(event_class);
        (*m_it)[*names++].get(zenith_angle);
        if( ! isFinite(zenith_angle) || zenith_angle<1e-10 ){ // latter seems to be what fits gives?
            zenith_angle=180.; // will be cut
        }
        try{
        (*m_it)[*names++].get(ctbclasslevel);
        }catch(const std::exception&){}

        if( m_selectid) { // check for source id only if requested
            (*m_it)[*names++].get(source);
        }

        if( !isFinite(energy) || !isFinite(dec) || !isFinite(ra) ){
            std::stringstream s;
            s << "Bad data: time = " << std::setprecision(10)<< time;
            s << "\n\tenergy, ra, dec: " << energy <<", " << ra<<", " << dec;
            throw std::runtime_error(s.str());
        }
        if(! m_fits){
            // A root file: apply standard cuts
            if( energy==0) event_class=99; // reject if not well measured
            else {
                event_class = event_class>4? 0 : 1;  // front/back map to event class 0/1
                double class_level; 
                (*m_it)["CTBClassLevel"].get(class_level);
                if( class_level<Data::class_level()) event_class=99;
                else{
                    // gets S/C pointing info for this event
                    (*m_it)["PtRaz" ].get(raz);
                    (*m_it)["PtDecz"].get(decz);
                    (*m_it)["PtRax" ].get(rax);
                    (*m_it)["PtDecx"].get(decx);                   
                }
            }
        }else{
            ///@todo: process FT2 to extact pointing information as above.
            if(!Data::historyfile().empty()&&Data::inTimeRange(time)) {
                astro::PointingInfo pi = Data::get_pointing(time);
                raz = pi.zAxis().ra();
                decz = pi.zAxis().dec();
                rax = pi.xAxis().ra();
                decx = pi.xAxis().dec();
            }
        }

        return Photon(astro::SkyDir(ra, dec), energy, time, event_class , source, 
            SkyDir(raz,decz),SkyDir(rax,decx), zenith_angle, ctbclasslevel);
    }


    EventList::Iterator EventList::begin()
    { 
        return Iterator(m_itbegin, m_fits, m_selectid);
    }

    EventList::Iterator EventList::end()
    {
        return Iterator(m_itend, m_fits);
    }

    //ROOT event extraction
    Photon events(std::vector<float>& row) {
        float ra(0), dec(0), energy(0); // photon info
        float raz(0), decz(0), rax(90), decx(0); // sc orientation: default orthogonal
        double time(0);
        double zenith_angle(0); //TODO: set this?
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
            time = row[3];
            event_class = static_cast<int>(row[4]);
            event_class = event_class>4? 0 : 1;  // front/back map to event class 0/1
            energy = row[2];
            class_level = static_cast<int>(row[9]);
            source_id = static_cast<int>(row[10]);
            // NB. Selecting wired-in class definition (transient, source, diffuse
            if( class_level < Data::class_level()) event_class=99;
            else{
                ra = row[0];
                dec = row[1];
                raz = row[5];
                decz = row[6];
                rax = row[7];
                decx = row[8];
            }
        }
        Photon p(astro::SkyDir(ra, dec), energy, time, event_class , source_id, astro::SkyDir(raz,decz),astro::SkyDir(rax,decx),zenith_angle,class_level);
        astro::Photon ap = p.transform(Data::get_rot().inverse());
        return Photon(ap.dir(),ap.energy(),ap.time(),ap.eventClass(),source_id,astro::SkyDir(raz,decz),astro::SkyDir(rax,decx),zenith_angle, class_level);
    }

    // default binner to use
    skymaps::PhotonBinner* binner( new skymaps::PhotonBinner() );
} // anon namespace

Data::Data(const embed_python::Module& setup)
{
    std::string pixelfile("")
        , output_pixelfile("");

    static std::string prefix("Data.");
    setup.getValue(prefix+"pixelfile", pixelfile, "");

    if(!pixelfile.empty()){
        try {
            m_data = new BinnedPhotonData(pixelfile, "PHOTONMAP" );
        } catch( const std::exception& ){
            m_data = new BinnedPhotonData(pixelfile, "BANDS" );
        }
        return;
    }

    int  event_class, source_id;
    std::vector<std::string> filelist;
    std::string history_filename;
    std::vector<double>alignment;

    setup.getList(prefix+"files", filelist);
    setup.getValue(prefix+"event_class", event_class);
    setup.getValue(prefix+"source_id",  source_id);
    setup.getValue(prefix+"start_time", m_start, -1);
    setup.getValue(prefix+"stop_time" , m_stop, -1);
    setup.getValue(prefix+"output_pixelfile", output_pixelfile, "");
    setup.getValue(prefix+"history",   history_filename, "");
    if( !history_filename.empty()){
        std::cout << "loading history file: " << history_filename << std::endl;
        Data::setHistoryFile(history_filename);
    }
    setup.getList(prefix+"LATalignment", alignment);
    if( alignment.size()>0) Data::set_rot(alignment);

    double bins_per_decade(0); 
    setup.getValue(prefix+"bins_per_decade", bins_per_decade, bins_per_decade);

    std::vector<double>energy_bins;
    setup.getList(prefix+"energy_bins", energy_bins);

    if( energy_bins.size()> 0 ){
        m_data = new BinnedPhotonData(skymaps::PhotonBinner(energy_bins));
    }else{
        m_data = new BinnedPhotonData(bins_per_decade);
    }

    if( filelist.empty()){
        throw std::invalid_argument("Data: no data specified: expected either a pixel file or a list of data files");
    }
    for( std::vector<std::string>::const_iterator it = filelist.begin(); 
        it !=filelist.end(); ++it)
    {
        const std::string& inputFile(*it);
        add(inputFile, event_class, source_id);
        addgti(inputFile);
    }
    if( !output_pixelfile.empty() ) {
        std::cout << "writing output pixel file :" << output_pixelfile ;
        m_data->write(output_pixelfile);
        std::cout << " with GTI range " << int(m_data->gti().minValue())<<"-"<<int(m_data->gti().maxValue())
            << std::endl;
    }
}

void Data::add(const std::string& inputFile, int event_type, int source_id)
{
    std::cout << "Loading data from file " << inputFile 
        << ", selecting event type " << event_type ;
    if( source_id>-1) {
        std::cout << " and source id " << source_id;
    }
    if( m_start>0 || m_stop>0 ) std::cout << " time range: (" << int(m_start) << " to " << int(m_stop) << ")";
    std::cout  << std::endl;

    int photoncount(m_data->photonCount()), pixelcount(m_data->pixelCount());
    if( inputFile.find(".root") != std::string::npos) {
        lroot(inputFile);
    }else {
        EventList photons(inputFile, source_id>-1);
        AddPhoton adder(*m_data, event_type, m_start, m_stop, source_id);

        std::for_each(photons.begin(), photons.end(), adder );
    }
    std::cout 
        << "photons found: "  << (m_data->photonCount() -photoncount) << " (total: " << m_data->photonCount() <<") "
        << "  pixels created: " << (m_data->pixelCount() -pixelcount) << " (total: " << m_data->pixelCount() << ") "
        << std::endl;

}

void Data::addgti(const std::string& inputFile)
{
    try
    {

        std::cout << "Loading gti info from file " << inputFile << "...";
        if(inputFile.find(".root") == std::string::npos) {

            skymaps::Gti tnew(inputFile); 
            m_data->addgti(tnew); 
            std::cout << " found interval " 
                << int(tnew.minValue())<<"-"<< int(tnew.maxValue())<<  std::endl;
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
                << ", total: " << m_data->gti().computeOntime()<< " s." 
                <<  std::endl;
        }

 
    }
    catch(const std::exception& e)
    {
        std::cerr << "\nCaught exception " << typeid(e).name() 
            << " \"" << e.what() << "\"" << std::endl;
        std::cerr << "Unable to access gti info from file " << inputFile << std::endl;
    }
}
Data::Data(const std::string& inputFile, int event_type, double tstart, double tstop, int source_id)
: m_data(new BinnedPhotonData(*binner))
, m_start(tstart), m_stop(tstop)
{
    add(inputFile, event_type, source_id);
    addgti(inputFile);
}

Data::Data(std::vector<std::string> inputFiles, int event_type, double tstart, double tstop, int source_id, std::string ft2file)
: m_data(new BinnedPhotonData(*binner))
, m_start(tstart), m_stop(tstop)
{

    for( std::vector<std::string>::const_iterator it = inputFiles.begin(); 
        it !=inputFiles.end(); ++it)
    {
        const std::string& inputFile(*it);
        add(inputFile, event_type, source_id);
        addgti(inputFile);
    }
}

Data::Data(const std::string & inputFile, const std::string & tablename)
: m_data(new BinnedPhotonData(inputFile, tablename))
{
    addgti(inputFile);
}

Data::~Data()
{
    delete m_data;
}



void Data::lroot(const std::string& inputFile) {
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
    bool flag(true);
    //for each entry  
    for(int i(0);i<entries&&(flag||m_start==-1);++i) {
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
        double time(row[3]), energy(row[2]), convlayer(row[4]);
        if( m_start<=0 || time >= m_start){
            if( m_stop<=0 || time < m_stop){
                m_data->addPhoton( events(row) );
            }else{
                flag=false; // here if off the end
            }
        }
        row.clear();
    }
    delete tf;
}

CLHEP::HepRotation Data::set_rot(double arcsecx,double arcsecy,double arcsecz) {
    CLHEP::HepRotation current = s_rot;
    s_rot = CLHEP::HepRotationX(arcsecx*M_PI/648000)*CLHEP::HepRotationY(arcsecy*M_PI/648000)*CLHEP::HepRotationZ(arcsecz*M_PI/648000);
    return current;
}
void Data::set_rot(std::vector<double> align) {
    using CLHEP::HepRotation;
    assert( align.size()==3); // otherwise bad interface setup
    static double torad(M_PI/(180*60*60));
    s_rot = HepRotationX(align[0]*torad) 
        * HepRotationY(align[1]*torad) 
        * HepRotationZ(align[2]*torad);
}

const CLHEP::HepRotation& Data::get_rot() {
    return s_rot;
}

const astro::PointingInfo& Data::get_pointing(double time) {
    return (*s_history)(time);
}

//! @brief define FT2 file to use for rotation
void Data::setHistoryFile(const std::string& history)
{
    s_ft2file = history;
    s_history = new astro::PointingHistory(history);
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