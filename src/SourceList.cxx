/** @file PointSourceLikelihood.h
@brief declaration of classes Source and SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SourceList.cxx,v 1.5 2008/07/06 06:41:34 burnett Exp $
*/


#include "pointlike/SourceList.h"
#include "pointlike/PointSourceLikelihood.h"
#include "skymaps/BinnedPhotonData.h"
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>

using namespace pointlike;
using astro::SkyDir;

namespace{
    // refit parameters: TODO, provide access
    double group_radius(2.0);


}
const skymaps::BinnedPhotonData * SourceList::s_data(0);
const skymaps::BinnedPhotonData * SourceList::data(){return s_data;}
void SourceList::set_data(const skymaps::BinnedPhotonData * data){s_data=data;}

Source::Source(const std::string& name, const astro::SkyDir& seed_dir, double TS)
: m_name(name)
, m_dir(seed_dir)
, m_seed_dir(seed_dir)
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
    if( TS==0 )  m_TS = m_fit->maximize();
}

double Source::localize(){
    m_fit->maximize();   // may not be needed
    m_sigma = m_fit->localize() ;
    m_dir = m_fit->dir();
    m_TS = m_fit->TS();
    return m_sigma;
}

double Source::moved()const{
    if(m_sigma>1 || m_sigma==0) return 0; // failed, or not localized
    return dir().difference(seed_dir())*180/M_PI;
}

void Source::header(std::ostream& out){
    out << std::left << std::setw(20) 
        <<"# name           ra        dec          TS    localization moved  [neighbor]\n";
//         vela             128.8465  -45.1863    248.20    0.0287    0.0000
}

void Source::info(std::ostream& out)const{
        out << std::left << std::fixed 
            << std::setw(15) << name() 
            << std::setprecision(4) << std::right
            << std::setw(10) << dir().ra() 
            << std::setw(10) << dir().dec() 
            << std::setprecision(2) 
            << std::setw(10) << TS() 
            << std::setprecision(4) 
            << std::setw(10) << sigma()
            << std::setw(10) << moved()
                ;
        if( neighbor()!=0){
            out << "  " << std::setw(20)<< std::left << (neighbor()->name());
        }
        out << std::endl;
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
        if( line.size()<5 || line[0]=='#' ) continue; // comment or empty
        double start; 
        std::stringstream buf(line); 
        std::string name; buf >> name;
        double ra, dec, TS(0);
        buf >> ra >> dec;
        if( !buf.eof()) buf >> TS;
        push_back( Source(name, SkyDir(ra, dec), TS));
    }
}

class GreaterTS{ // simple functor to sort the list by TS
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
    Source::header(out);
    for( const_iterator it= begin(); it!= end(); ++it){it->info(out); }
}
void SourceList::dump(const std::string& outfilename)const
{
    std::ofstream out(outfilename.c_str());
    dump( out );
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
    static int size(30); // 11 is default: make them stand out
    out << "global color="<< color 
        << " font=\"helvetica 10 normal\" select=1 edit=1 move=0 delete=1 include=1 fixed=1 width=2;fk5;"
        << std::fixed << std::setprecision(4) << std::endl;
    int n(0);
    for( const_iterator it = begin(); it != end();  ++it)  {
        const Source& cand( * it);
        if(cand.TS()< tsmin) continue;
        out << "point("<< cand.dir().ra()<< ","<<cand.dir().dec() <<") point=cross "
            << size << " # text={TS=" 
            << int(cand.TS()+0.5) << "};\n";
        ++n;
    }
    out.close();
    std::cout << "Wrote "<< n << " entries to  reg file "<< filename << ", with TS >= " <<tsmin << std::endl;
}
