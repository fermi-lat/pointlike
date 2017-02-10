/** @file EventList.cxx 
@brief declaration of the EventList wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/EventList.cxx,v 1.28 2016/04/21 00:22:22 wallacee Exp $
*/

#include "EventList.h"

#include "skymaps/Gti.h"

#include "astro/PointingInfo.h"
#include "astro/GPS.h"
#include <iostream>
#include <iomanip>
#include <cstdlib>
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
int count(100);

void AddPhoton::operator()(const Photon& gamma)
{
    m_found++;
    int event_class = gamma.eventClass();
    int event_type = m_binner.event_type_index(gamma.event_type()); // convert to 0-6 value
    int sourceid = gamma.source();

#if 0 // disable event type or class or timing selection for now
    if( m_select>-1){
        if( m_data_pass < 8){
            if( event_class!= m_select) return;
        }else{
            if( (event_type & (1<<m_select))==0 ) return;
        }
    }
   
    if( m_source>-1 && sourceid != m_source){
        if(count>0){std::cout<< "wrong source" << std::endl; }
        return;
    }
#else
    if(count==100){
        std::cout << "AddPhoton: not checking event type or sourceid" << std::endl;
    }
#endif
    if (count-- >0){
            std::cout << "AddPhoton: time, type, energy, theta, zenith: " << gamma.time() <<" "
            << event_type    <<" " 
            << gamma.energy() << " "
            << gamma.theta() << " " << gamma.zenith_angle()<<"..." ;
        }

    // timing: either start/stop interval, or a Gti object
    if( m_use_gti ){
        if( !m_gti.accept(gamma.time()) ) return;}
    else{
            if( m_start>0   && gamma.time()<m_start ||  m_stop>m_start && gamma.time()>m_stop){
                if(count>0){ std::cout << "Failed time check" << std::endl;}
                 return;
            }
    }

    // zenith angle cut, unless in Earth coordinates
    double zenith_angle = gamma.zenith_angle();
    double energy(gamma.energy());
    if( use_earth==0 ){
        if( event_type>1 ){
            // Make zenith angle cut: wire in the values for PSF selection
            if(count>0){std::cout << "PSF zenith chk...";}
            if ( energy< 100 && zenith_angle > 80 
            || energy  < 300 && zenith_angle > 90 
            || energy  <1000 && zenith_angle > 100 
            || zenith_angle >105 ){
                if (count>0){std::cout << "fail zcut" << std::endl; }
                return;
                }     
        } else { 
            // Doing FB: apply theta cut here, not iwn PSF
            if( gamma.theta()>Data::theta_cut() || gamma.theta()<=Data::theta_min_cut() ){
                if(count>0){std::cout << "fail theta cut" << std::endl;}
                return;
            }
            if ( zenith_angle > 100 ){
                if (count>0){std::cout << "fail zcut" << std::endl; }
                return;
            } 
        }
    }
    m_kept++;
    if (count>0){ std::cout << "ok" << std::endl;  }
    
    m_map.addPhoton(gamma);
}


EventList::EventList(const std::string infile, bool selectid, bool use_mc_energy,
                     std::string table_name)
                     : m_fits(true)
                     , m_selectid(selectid)
                     , m_use_mc_energy(use_mc_energy)
                     //, m_pass7(true), m_evclass_bitarray(true)
					 , m_data_pass(8)
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
    if(m_fits) {
        try{
            double dif;
            // (*m_itbegin)["DIFRSP1"].get(dif);
            (*m_itbegin)["CTBCLASSLEVEL"].get(dif);
            //m_pass7=false;
			m_data_pass=6;
        } catch (const std::exception& ){}
    }

    int evclsver(0); //Is this ever used?
    
    const tip::Header & header(m_table->getHeader());
    try {
      header["EVCLSVER"].get(evclsver);
    } catch(tip::TipException) {
      // keyword missing so use default value
    }
    std::string pass_ver;
    try {
      header["PASS_VER"].get(pass_ver);
    } catch(tip::TipException) {
      // keyword missing so use default value
      pass_ver = "NONE";
    }

    //if (pass_ver == "NONE" || pass_ver.substr(0, 2) == "P7") {
    //  m_evclass_bitarray = false;
    //}
	
	//NB: This assumes that all Pass 7 and higher FT1s will
	//have the pass_ver header key. -EEW
	m_data_pass = pass_ver=="NONE"?6:atoi(&pass_ver[1]);

}
EventList::EventList()
{
    std::cout << "Default ctor!" << std::endl;
}
EventList::~EventList()
{
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
    int event_type(0); // NEW

    // NB: the internal variables (event_class, ctbclasslevel) no longer match FT1 names.
    // event_class == CONVERSION_TYPE, ctbclasslevel = "EVENT_CLASS"

    // FT1 names
    static std::string fits_names[]={"RA", "DEC", 
        "ENERGY", "TIME", "CONVERSION_TYPE", "ZENITH_ANGLE","THETA", "EVENT_CLASS", "MC_SRC_ID"};
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

    if (m_data_pass>=8) {
       //std::cout << "Testing ..." << std::endl;
#if 1
      tip::BitStruct tip_event_class;
      (*m_it)[*names++].get(tip_event_class);
      ctbclasslevel = static_cast<int>(tip_event_class);
#endif 
#if 1 // seems to cause exception??
      // NEW STUFF for EVENT_TYPE -- convert bit array to int, like EVENT_CLASS
      tip::BitStruct tip_event_type;
      (*m_it)["EVENT_TYPE"].get(tip_event_type);
      event_type = static_cast<int>(tip_event_type);
#endif
     
       //std::cout << "End Testing: event type, class: " << event_type << ", " << ctbclasslevel << std::endl;
      
    } else if(m_data_pass==7)  {
        (*m_it)[*names++].get(ctbclasslevel);
		event_type = event_class;  //
    } else  {
      (*m_it)["CTBCLASSLEVEL"].get(ctbclasslevel);
	  names++; // Still need to keep the iterator lined up...
	  event_type = event_class;
    }

    if( m_selectid) { // check for source id only if requested
      (*m_it)[*names++].get(source); // source id for monte carlo testing
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
    
    return Photon(astro::SkyDir(ra, dec), energy, time, event_class , source, // note sets event_class in astro::Photon
        SkyDir(raz,decz),SkyDir(rax,decx), zenith_angle, theta, ctbclasslevel
          , event_type //NEW
            ); 
}


EventList::Iterator EventList::begin()
{ 
  //return Iterator(m_itbegin, m_fits, m_selectid, m_use_mc_energy, m_pass7, m_evclass_bitarray);
  return Iterator(m_itbegin, m_fits, m_selectid, m_use_mc_energy, m_data_pass);
}

EventList::Iterator EventList::end()
{
    return Iterator(m_itend, m_fits);//Does it matter if data_pass is given to this one?
}

