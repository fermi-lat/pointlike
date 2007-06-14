/** @file Data.cxx
@brief implementation of Data

$Header: /nfs/slac/g/glast/ground/cvs/users/burnett/pointlike/src/Data.cxx,v 1.1.1.1 2007/06/10 01:05:26 burnett Exp $

*/


#include "pointlike/Data.h"

#include "astro/SkyDir.h"

#include "map_tools/PhotonMap.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>


using namespace astro;
using namespace pointlike;

namespace {

    class AddPhoton: public std::unary_function<astro::Photon, void> {
    public:
        AddPhoton (map_tools::PhotonMap& map, int select, int source=-1)
            : m_map(map), m_select(select), m_source(source)
        {}
        void operator()(astro::Photon& gamma)
        {
            int event_class = gamma.eventClass();
            int sourceid = gamma.source();

            if( m_select>-1 && event_class!= m_select) return;
            if( m_source>-1 && sourceid != m_source)return;
            m_map.addPhoton(gamma);
        }
        map_tools::PhotonMap& m_map;
        int m_select;
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
        @param infile name of the input ROOT file, expect to be FT1 
        */
        EventList(const std::string infile, 
            const std::string table_name="EVENTS");

        ~EventList();

        // make it a container by implementing a forward iterator
        class Iterator {
        public:
            Iterator(tip::Table::ConstIterator it):m_it(it){}
            astro::Photon operator*()const;             ///< dereference
            tip::Table::ConstIterator operator++(){return ++m_it;} ///< increment operator
            bool operator!=(const Iterator& other)const{return other.m_it!=m_it;}
        private:
            tip::Table::ConstIterator m_it;
        };

        /// return iterator to access 
        Iterator begin();
        Iterator end();


    private:
        tip::Table::ConstIterator m_itbegin, m_itend;
        const tip::Table * m_table;
    };

    EventList::EventList(const std::string infile, 
        const std::string table_name)
    {
        // connect to  input data
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

    astro::Photon EventList::Iterator::operator*()const
    {
        float ra, dec, energy;
        double time;
        short layer = 0;
        int event_class;
        int source;

        // FT1 names
        static std::string names[]={"RA", "DEC", 
            "ENERGY", "TIME", "EVENT_CLASS","MC_SRC_ID"};

        int n=0;
        (*m_it)[names[n++]].get(ra);
        (*m_it)[names[n++]].get(dec);
        (*m_it)[names[n++]].get(energy);
        (*m_it)[names[n++]].get(time);
        (*m_it)[names[n++]].get(event_class);
        (*m_it)[names[n++]].get(source);

        return astro::Photon(astro::SkyDir(ra, dec), energy, time, event_class , source);
    }


    EventList::Iterator EventList::begin()
    { 
        return Iterator(m_itbegin);
    }

    EventList::Iterator EventList::end()
    {
        return Iterator(m_itend);
    }


} // anon namespace


void Data::add(const std::string& inputFile, int event_type, int source_id)
{
    std::cout << "Loading data from file " << inputFile 
        << ", selecting event type " << event_type ;
    if( source_id>0 ) {
        std::cout << " and source id " << source_id;
    }
    std::cout  << std::endl;

    int photoncount(m_data->photonCount()), pixelcount(m_data->pixelCount());

    EventList photons(inputFile);
    std::for_each(photons.begin(), photons.end(), AddPhoton(*m_data, event_type, source_id));
    std::cout << "photons found: "<< (m_data->photonCount() -photoncount) << " (total: " << m_data->photonCount() <<") "
        <<"  pixels created: " << (m_data->pixelCount() -pixelcount) << " (total: " << m_data->pixelCount() << ") "
        <<std::endl;

}
Data::Data(const std::string& inputFile, int event_type, int source_id)
: m_data(new map_tools::PhotonMap())
{
    add(inputFile, event_type, source_id);
}

Data::Data(std::vector<std::string> inputFiles, int event_type, int source_id)
: m_data(new map_tools::PhotonMap())
{
    for( std::vector<std::string>::const_iterator it = inputFiles.begin(); 
        it !=inputFiles.end(); ++it)
    {
        const std::string& inputFile(*it);
        add(inputFile, event_type, source_id);
    }
}


