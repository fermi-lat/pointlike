/** @file PointSourceLikelihood.h
@brief declaration of classes Source and SourceList

$Header$
*/


#include "pointlike/SourceList.h"
#include "pointlike/PointSourceLikelihood.h"
#include "skymaps/BinnedPhotonData.h"
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace pointlike;
using astro::SkyDir;

namespace{
    // refit parameters: TODO, provide access
    double group_radius(2.0);


}
const skymaps::BinnedPhotonData * SourceList::s_data(0);
const skymaps::BinnedPhotonData * SourceList::data(){return s_data;}
void SourceList::set_data(const skymaps::BinnedPhotonData * data){s_data=data;}

Source::Source(const std::string& name, const astro::SkyDir& sdir, double TS)
: m_name(name)
, m_dir(sdir)
, m_fit(0)
, m_TS(TS)
, m_sigma(0)
, m_neighbor(0)
{
    if( SourceList::data()==0){
        throw std::invalid_argument("Source::Source: no data set");
    }
    if(m_fit==0){
        m_fit = new PointSourceLikelihood(*SourceList::data(), name, m_dir);
    }
    // inital maximize unless TS already set.
    if( TS=0 )  m_TS = m_fit->maximize();
}

double Source::localize(){
    m_fit->maximize();   // may not be needed
    m_sigma = m_fit->localize() ;
    m_dir = m_fit->dir();
    m_TS = m_fit->TS();
    return m_sigma;
}

Source::~Source()
{
    //delete m_fit; //TODO: copy problem with this
}

SourceList::SourceList(std::map<std::string, std::vector<double> > input_list)
{
    std::map<std::string, std::vector<double> >::const_iterator it(input_list.begin());
    for( ; it!= input_list.end(); ++it ){

        std::string name(it->first);
        const std::vector<double>& v(it->second);
        double TS( v.size()>2? v[2] : 0);
        push_back( Source(name, SkyDir(v.at(0), v.at(1)),TS));
    }

}
SourceList::SourceList(const std::string& filename)
{
    std::ifstream input_file;
    input_file.open(filename.c_str());

    if(!input_file.is_open()) {
        std::cerr << "ERROR:  Unable to open:  " << filename << std::endl;
        throw std::invalid_argument("SourceList: could not open file");
    }
    while (!input_file.eof()){
        std::string line; std::getline(input_file, line);
        if( line[0]=='#' ) continue;
        double start; 
        std::stringstream buf(line); 
        std::string name; buf >> name;
        double ra, dec, TS(0);
        buf >> ra >> dec;
        if( !buf.eof()) buf >> TS;
        push_back( Source(name, SkyDir(ra, dec), TS));
    }
}

class GreaterTS{
public:
    GreaterTS(){}
    bool operator()(const Source& a, const Source& b){
        return a.TS() > b.TS();
    }
private:

};

void SourceList::sort_TS()
{
    std::sort( begin(), end(), GreaterTS() );
}

void SourceList::dump(std::ostream& out)const
{
    using std::setw;
    using std::setprecision;
    out << std::left << std::setw(20) <<"name" << "     TS     error    ra     dec\n";


    for( const_iterator it= begin(); it!= end(); ++it){
        const Source& s(*it);
        out << std::left << std::setw(20) << s.name() 
                << std::setprecision(2) << std::setw(8) << std::fixed << std::right
                << s.TS() 
                << std::setprecision(4) 
                << std::setw(10) << s.sigma()
                << std::setw(10) << s.dir().ra() 
                << std::setw(10) << s.dir().dec() 
                ;
        if( s.neighbor()!=0){
            
            out << std::setw(20) << (s.neighbor()->name());
        }
        out << std::endl;
    }
}

void SourceList::refit()
{
    for( iterator it= begin(); it!= end(); ++it){
        Source& cand(*it);
        SkyDir currentpos( cand.dir() );

        // See if we have already found a candidate near this location
        iterator neighbor;
        double max_value(0.);
        cand.fit()->clearBackgroundPointSource(); // clear any existing
        cand.set_neighbor(0);

        if (group_radius > 0.0) {
            for( iterator it2 = begin();  it2 != end(); ++it2) {
                Source& other( *it2);
                double diff = currentpos.difference(other.dir())*180/M_PI;
                if( diff <= group_radius &&  other.TS() > cand.TS() ) {
                    if( other.TS() > max_value) {
                        max_value = other.TS();
                        neighbor = it2;
                    }
                }
            }
        }
        if (max_value>0.) {
            // Add strong neighbor to this candidate's background
            cand.fit()->addBackgroundPointSource(neighbor->fit());
            cand.set_neighbor( &*neighbor );
        }
        // adjust position to maximize likelhood
        cand.localize();
    }
}

void SourceList::createRegFile(std::string filename, std::string color, double tsmin)const
{
    std::ofstream out(filename.c_str());
    out << "global color="<< color 
        << " font=\"helvetica 10 normal\" select=1 edit=1 move=0 delete=1 include=1 fixed=0 width=2;fk5;"
        << std::fixed << std::setprecision(4) << std::endl;
    int n(0);
    for( const_iterator it = begin(); it != end();  ++it)  {
        const Source& cand( * it);
        if(cand.TS()< tsmin) continue;
        out << "cross point("<< cand.dir().ra()<< ","<<cand.dir().dec() <<") # text={TS=" 
            << int(cand.TS()+0.5) << "};\n";
        ++n;
    }
    out.close();
    std::cout << "Wrote "<< n << " entries to  reg file "<< filename << ", with TS >= " <<tsmin << std::endl;
}
