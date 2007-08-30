/** @file CalData.cxx
@brief implementation of CalData

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/CalData.cxx,v 1.1 2007/08/27 23:24:00 mar0 Exp $

*/


#include "pointlike/CalData.h"

#include "astro/SkyDir.h"
#include "astro/PointingTransform.h"

#include "map_tools/PhotonMap.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <cassert>


using namespace astro;
using namespace pointlike;
using namespace CLHEP;
double CalData::s_scale[4]={1.0, 1.86, 1.0, 1.0}; // wired in for front, back !!



namespace {
    /** @class Photon
    @brief derive from astro::Photon to allow transformation

    */

    class Photon : public astro::Photon {
    public:
        /// ctor sets photon, also rotation info.
        Photon(const astro::SkyDir& dir, double energy, 
            double time, int event_class, int source
            , const astro::SkyDir&scz, const astro::SkyDir& scx
            )
            : astro::Photon(dir, energy, time, event_class, source)
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

            // account for worse energy resolution of back events by simply scaling up the energy.
            int evtclass(eventClass());
            assert( evtclass>=0 && evtclass<3); // just a check: should be 0 or 1 if not DC2
            double 
                emeas(energy()),                    // measured energy
                eff_energy( emeas/CalData::s_scale[evtclass] ); // "effective" energy
            return astro::Photon(astro::SkyDir(transformed), eff_energy,time(),evtclass, source());
        }
        HepRotation getRotation() const {return m_rot;}

    private:
        HepRotation m_rot;

    };

    class AddPhoton: public std::unary_function<astro::Photon, void> {
    public:
        AddPhoton (map_tools::PhotonMap& map, HepRotation& hr, std::vector<astro::SkyDir>& sources, int select, int source )
            : m_map(map), m_rot(hr), m_select(select), m_source(source), m_sources(sources)
        {}
        double rescale(double energy, int eventclass)
        {
            if( eventclass==0)     return energy;   // front events: pass through
            //    back events
            if( energy< 6500) return energy/1.85; // below 6.5 GeV: rescale
            if( energy<21000) return 5000;        // below 21 GeV: 5 GeV is level 11
            return 10000;                         // above 21 Gev: 10 GeV 1s level 12
        }
        void operator()(const Photon& gamma)
        {

            int event_class = gamma.eventClass();
            int sourceid = gamma.source();

            if( m_select>-1 && event_class!= m_select) return;
            if( m_source>-1 && sourceid != m_source)return;

            astro::SkyDir sd;
            double min = 1e9;
            for(std::vector<astro::SkyDir>::const_iterator it = m_sources.begin();it!=m_sources.end();++it) {
                double newdiff = it->difference(gamma());
                if(newdiff<min) {
                    sd = *it;
                    min = newdiff;
                }
            }
            if(min<7*M_PI/180&&gamma.energy()>=100) {
                double energy(gamma.energy());
                // rescale according to event class
                HepRotation glast_rot = gamma.getRotation();
                Hep3Vector gammavec = gamma();
                Hep3Vector sdvec = sd();
                Hep3Vector transformed = glast_rot.inverse()*(sdvec-gammavec);
                m_map.addPhoton(astro::Photon(astro::SkyDir(Hep3Vector(1,transformed.y(),transformed.x())),rescale(energy,event_class),gamma.time(),event_class,sourceid));
            }
        }
        map_tools::PhotonMap& m_map;
        const HepRotation m_rot; 
        int m_select;
        int m_source;
        const std::vector<astro::SkyDir> m_sources;
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
        @param table_name must be "EVENTS" for FT1, or "MeritTuple"
        */
        EventList( std::string infile, 
            std::string table_name="EVENTS");

        ~EventList();

        // make it a container by implementing a forward iterator
        class Iterator {
        public:
            Iterator(tip::Table::ConstIterator it, bool fits):m_it(it), m_fits(fits){}
            Photon operator*()const;             ///< dereference
            tip::Table::ConstIterator operator++(){return ++m_it;} ///< increment operator
            bool operator!=(const Iterator& other)const{return other.m_it!=m_it;}
        private:
            tip::Table::ConstIterator m_it;
            bool m_fits;
        };

        /// return iterator to access 
        Iterator begin();
        Iterator end();


    private:
        tip::Table::ConstIterator m_itbegin, m_itend;
        const tip::Table * m_table;
        bool m_fits;
    };

    EventList::EventList(const std::string infile, 
        std::string table_name)
        :m_fits(true)
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
        int event_class;
        int source;

        // FT1 names
        static std::string fits_names[]={"RA", "DEC", 
            "ENERGY", "TIME", "EVENT_CLASS","MC_SRC_ID"};
        // corresponging names in the ROOT MeritTuple
        static std::string root_names[]={"FT1Ra", "FT1Dec", 
            "CTBBestEnergy", "EvtElapsedTime"
            , "FT1ConvLayer", "McSourceId"};

        std::string *names = m_fits?  fits_names : root_names;

        (*m_it)[*names++].get(ra);
        (*m_it)[*names++].get(dec);
        (*m_it)[*names++].get(energy);
        (*m_it)[*names++].get(time);
        (*m_it)[*names++].get(event_class);
        (*m_it)[*names++].get(source);

        if(! m_fits){
            // A root file: apply standard cuts
            event_class = event_class>4? 0 : 1;  // front/back map to event class 0/1
            if( energy==0) event_class=99; // reject if not well measured
            else {
                double gammaprob; 
                (*m_it)["CTBGAM"].get(gammaprob);
                if( gammaprob<0.1) event_class=99;
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
        }

        return Photon(astro::SkyDir(ra, dec), energy, time, event_class , source, 
            SkyDir(raz,decz),SkyDir(rax,decx));
    }


    EventList::Iterator EventList::begin()
    { 
        return Iterator(m_itbegin, m_fits);
    }

    EventList::Iterator EventList::end()
    {
        return Iterator(m_itend, m_fits);
    }


} // anon namespace


void CalData::add(const std::string& inputFile, int event_type, int source_id)
{
    std::cout << "Loading data from file " << inputFile 
        << ", selecting event type " << event_type ;
    if( source_id>0 ) {
        std::cout << " and source id " << source_id;
    }
    std::cout  << std::endl;

    int photoncount(m_data->photonCount()), pixelcount(m_data->pixelCount());

    EventList photons(inputFile);
    HepRotation hr(0,0,0);
    AddPhoton adder(*m_data, hr, m_sources, event_type, source_id);
    EventList::Iterator it = photons.begin();
    double starttime = (*(photons.begin())).time();
    for(EventList::Iterator it = photons.begin();it!=photons.end()&&((((*it).time()-starttime<=m_stop)&&((*it).time()-starttime>=m_start))||m_start==-1);++it) {
        adder(*it);
    }

    std::cout << "photons found: "<< (m_data->photonCount() -photoncount) << " (total: " << m_data->photonCount() <<") "
        <<"  pixels created: " << (m_data->pixelCount() -pixelcount) << " (total: " << m_data->pixelCount() << ") "
        <<std::endl;

}
CalData::CalData(const std::string& inputFile, std::vector<astro::SkyDir> sources, int start, int stop, int event_type, int source_id)
: m_data(new map_tools::PhotonMap())
, m_ft2file("")
, m_sources(sources)
, m_start(start)
, m_stop(stop)
{
    add(inputFile, event_type, source_id);
}

CalData::CalData(std::vector<std::string> inputFiles, std::vector<astro::SkyDir> sources, int start, int stop, int event_type, int source_id, std::string ft2file )
: m_data(new map_tools::PhotonMap())
, m_ft2file(ft2file)
, m_sources(sources)
, m_start(start)
, m_stop(stop)
{

    for( std::vector<std::string>::const_iterator it = inputFiles.begin(); 
        it !=inputFiles.end(); ++it)
    {
        const std::string& inputFile(*it);
        add(inputFile, event_type, source_id);
    }
}


