/** @file EventList.cxx 
@brief declaration of the EventList wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/EventList.cxx,v 1.7 2009/02/26 23:54:18 kerrm Exp $
*/

#include "EventList.h"

#include "skymaps/Gti.h"

#include "astro/PointingInfo.h"
#include "astro/GPS.h"

#include <iomanip>
using namespace pointlike;
using skymaps::BinnedPhotonData;




using astro::SkyDir;

namespace{

    // file-scope pointer to gps instance
    astro::GPS* gps (astro::GPS::instance()); 

    int use_earth(0);  // flag to convert to Earth coordinates FIXME!

#ifdef WIN32
#include <float.h> // used to check for NaN
#else
#include <cmath>
#endif



    bool isFinite(double val) {
        using namespace std; // should allow either std::isfinite or ::isfinite
#ifdef WIN32 
        return (_finite(val)!=0);  // Win32 call available in float.h
#else
        return (isfinite(val)!=0); // gcc call available in math.h 
#endif
    }
    bool apply_correction(false); // set true to apply correction. disabled now


}// anon namespace

void AddPhoton::operator()(const Photon& gamma)
    {
        m_found++;
        int event_class = gamma.eventClass();
        int sourceid = gamma.source();

        if( m_select>-1 && event_class!= m_select) return;
        if( m_start>0   && gamma.time()<m_start ||  m_stop>m_start && gamma.time()>m_stop) return;
        if( m_source>-1 && sourceid != m_source)return;

        // theta cut: define FOV
        if( gamma.theta()>Data::theta_cut() ) return;

        // zenith angle cut, unless in Earth coordinates
        if( ! isFinite(gamma.zenith_angle()) ) return; // catch NaN?
        if( use_earth==0 && gamma.zenith_angle()> Data::zenith_angle_cut()) return;
        int class_level( gamma.class_level() );
        if( class_level< pointlike::Data::class_level() ) return; // select class level
        m_kept++;


        if( apply_correction ){
            // using GPS to make the correction, including aberration
            // 10Nov08: this does not seem to work, but now moot. leave code for testing with debugger
            gps->enableAberration();
            gps->setAlignmentRotation(Data::get_rot(gamma.time()));
            SkyDir fixed(gps->correct(gamma.dir(), gamma.time()));


#if 1 // to look at a few values with the debugger
            double ra(gamma.ra()), dec(gamma.dec());
            double raf(fixed.ra()), decf(fixed.dec());
#endif
#if 1 // alternate code, from old version
            SkyDir oldfix( gamma.transform(Data::get_rot(gamma.time())) );
            double raf2(oldfix.ra()), decf2(oldfix.dec());
            fixed=oldfix; // replace!!!

#endif


            m_map.addPhoton(astro::Photon(fixed, gamma.energy(),gamma.time(),gamma.eventClass()));
        }else{
            // no correction: just add the photon
            m_map.addPhoton(gamma);
        }
    }


EventList::EventList(const std::string infile, bool selectid, bool use_mc_energy,
                     std::string table_name)
                     : m_fits(true)
                     , m_selectid(selectid)
                     , m_use_mc_energy(use_mc_energy)
{
    if( infile.find(".root") != std::string::npos) {
        table_name = "MeritTuple"; 
        m_fits=false;
    }
    // connect to input data
    m_table = tip::IFileSvc::instance().readTable(infile, table_name, "");

    // save the iterators
    m_itbegin= m_table->begin();
    m_itend = m_table->end();

}
EventList::EventList()
{
    std::cout << "Default ctor!" << std::endl;
}
EventList::~EventList()
{
    // seems to create crash
    std::cout << "deleting table" << std::endl;
    delete m_table;
}

Photon EventList::Iterator::operator*()const
{
    float ra, dec, energy, mc_energy; // photon info
    float raz(0), decz(0), rax(90), decx(0); // sc orientation: default orthogonal
    double time;
    double zenith_angle;
    double theta;
    int event_class, ctbclasslevel(1);
    int source(-1);

    // FT1 names
    static std::string fits_names[]={"RA", "DEC", 
        "ENERGY", "TIME", "EVENT_CLASS", "ZENITH_ANGLE","THETA", "CTBCLASSLEVEL", "MC_SRC_ID"};
    // corresponging names in the ROOT MeritTuple
    static std::string root_names[]={"FT1Ra", "FT1Dec", 
        "CTBBestEnergy", "EvtElapsedTime"
        , "FT1ConvLayer", "FT1ZenithAngle", "FT1Theta", "McSourceId"};

    std::string *names = m_fits?  fits_names : root_names;

    (*m_it)[*names++].get(ra);
    (*m_it)[*names++].get(dec);

    if( use_earth>0) {
        double az, zen;
        // convert to zenith-centered coordinate system: z axis to zenith, x to East
        (*m_it)["EARTH_AZIMUTH_ANGLE"].get(az);
        (*m_it)["ZENITH_ANGLE"].get(zen);
        dec = 90-zen; // latitude
        ra  = az; // 0 is North, 90 East, etc. (not Petry's convention)
    }

    (*m_it)[*names++].get(energy);
    (*m_it)[*names++].get(time);
    (*m_it)[*names++].get(event_class);
    (*m_it)[*names++].get(zenith_angle);
    if( ! isFinite(zenith_angle) || zenith_angle<1e-10 ){ // latter seems to be what fits gives?
        zenith_angle=180.; // will be cut
    }
    (*m_it)[*names++].get(theta);
    try{
        (*m_it)[*names++].get(ctbclasslevel);
    }catch(const std::exception&){
        ctbclasslevel=3;
    }
    if( m_selectid) { // check for source id only if requested
        (*m_it)[*names++].get(source);
    }
	if( m_use_mc_energy ) { //get MC_ENERGY and replace ENERGY, if using it
		(*m_it)["MCENERGY"].get(mc_energy);
		energy = mc_energy;
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
        SkyDir(raz,decz),SkyDir(rax,decx), zenith_angle, theta, ctbclasslevel);
}


EventList::Iterator EventList::begin()
{ 
    return Iterator(m_itbegin, m_fits, m_selectid, m_use_mc_energy);
}

EventList::Iterator EventList::end()
{
    return Iterator(m_itend, m_fits);
}

