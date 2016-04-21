/** @file EventList.h 
@brief declaration of the EventList wrapper class

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/EventList.h,v 1.14 2015/06/26 14:34:06 burnett Exp $
*/


#ifndef pointlike_EventList_h
#define pointlike_EventList_h

#include "pointlike/Data.h"
#include "skymaps/Photon.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"
#include "astro/SkyDir.h"

#include "skymaps/BinnedPhotonData.h"
#include <string>

/** @class Photon
@brief derive from skymaps::Photon to allow transformation

*/

class Photon : public skymaps::Photon {
public:
    /// ctor sets photon, also rotation info.
    Photon(const astro::SkyDir& dir, double energy, 
        double time, int event_class, int source
        , const astro::SkyDir&scz, const astro::SkyDir& scx
        , double zenith_angle, double theta
        , int ctbclasslevel=pointlike::Data::class_level()  // note default 
        , int event_type=-1 //default invalid: expect to be filled from EVENT_TYPE field
        )
        : skymaps::Photon(dir, energy, time, event_class, source, event_type)
        , m_zenith_angle(zenith_angle)
        , m_theta(theta)
        , m_ctbclasslevel(ctbclasslevel)
    {
        CLHEP::Hep3Vector scy (scz().cross(scx()));
        m_rot = CLHEP::HepRotation(scx(), scy, scz());
    }
    /// make transformation in GLAST frame
    /// @param corr matrix that corrects direction in the GLAST frame
    skymaps::Photon transform(const CLHEP::HepRotation& corr)const
    {
        CLHEP::Hep3Vector 
            local( m_rot.inverse()*dir()),
            transformed( m_rot * corr * local );

        int evtclass(eventClass());
        if( evtclass<0 || evtclass>1){ // just a check: should be 0 or 1 for front/back
            //std::stringstream err;
            std::cerr << "Data::Photon::transform--invalid eventclass value, expect 0 or 1, got: "<< evtclass<< std::endl;
            //throw std::runtime_error(err.str());
        }
        double   emeas(energy());                    // measured energy
        return skymaps::Photon(SkyDir(transformed), emeas,time(),evtclass, source(), event_type());
    }
    double zenith_angle()const{return m_zenith_angle;}
    double theta()const{return m_theta;}
    int class_level()const{return m_ctbclasslevel;}
    int event_type_index(int event_type);
private:
    CLHEP::HepRotation m_rot;
    double m_zenith_angle;
    double m_theta;
    int m_ctbclasslevel;

};


class AddPhoton: public std::unary_function<skymaps::Photon, void> {
public:
    AddPhoton (skymaps::BinnedPhotonData& map, int select, 
        double start, double stop, int source,  skymaps::Gti gti=skymaps::Gti(),
        //bool pass7=true, )
		int data_pass=8)
        : m_map(map), m_select(select)
        , m_start(start), m_stop(stop), m_source(source)
        , m_found(0), m_kept(0)
        , m_gti(gti)
        , m_use_gti(gti.getNumIntervals()>0)
        //, m_pass7(pass7)
		, m_data_pass(data_pass)
    {}
    void operator()(const Photon& gamma);

    int found()const{return m_found;}
    int kept()const{return m_kept;}
    skymaps::BinnedPhotonData& m_map;
    int m_select;
    double m_start, m_stop;
    int m_source;
    int m_found, m_kept;
    skymaps::Gti m_gti;
    bool m_use_gti;//,m_pass7;
	int m_data_pass;
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
	@param use_mc_energy True if will be using MC_ENERGY
    @param table_name must be "EVENTS" for FT1, or "MeritTuple"
    */
    EventList( std::string infile, bool selectid=false, bool use_mc_energy=false,
        std::string table_name="EVENTS");

    EventList(); 
    ~EventList();

    // make it a container by implementing a forward iterator
    class Iterator {
    public:
  	// these traits needed  for STL functions like accumulate
		typedef const Photon& reference;
		typedef const Photon* pointer;
		typedef Photon value_type;
		typedef std::forward_iterator_tag iterator_category;
		typedef int difference_type; //??? needed, hope this kluge does not break anythin

    //Iterator(tip::Table::ConstIterator it, bool fits, bool selectid=false, bool use_mc_energy=false, bool pass7=false, bool evclass_bitarray=true)
    Iterator(tip::Table::ConstIterator it, bool fits, bool selectid=false, bool use_mc_energy=false, int data_pass=8)
            : m_it(it)
            , m_fits(fits)
            , m_selectid(selectid)
            , m_use_mc_energy(use_mc_energy)
			// Just specify the pass directly
            //, m_pass7(pass7)						
            //, m_evclass_bitarray(evclass_bitarray)
			, m_data_pass(data_pass) //Note default change to 8, need to verify that it doesn't break anything - EEW
        {}
        Photon operator*()const;             ///< dereference
        tip::Table::ConstIterator operator++(){return ++m_it;} ///< increment operator
        bool operator!=(const Iterator& other)const{return other.m_it!=m_it;}
    private:
        tip::Table::ConstIterator m_it;
        bool m_fits, m_selectid, m_use_mc_energy; //,m_pass7, m_evclass_bitarray;
		int m_data_pass;
    };

    /// return iterator to access 
    Iterator begin();
    Iterator end();

    void close();
    static std::string root_names[];
    //bool pass7() {return m_pass7;}
	int data_pass()const{return m_data_pass;}

    
private:
    tip::Table::ConstIterator m_itbegin, m_itend;
    const tip::Table * m_table;
    bool m_fits;
    bool m_selectid;
    bool m_use_mc_energy;
    //bool m_pass7;
    //bool m_evclass_bitarray;
	int m_data_pass;
};
#endif

