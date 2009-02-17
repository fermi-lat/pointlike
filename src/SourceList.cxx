/** @file SourceList.cxx
@brief implementation of classes Source and SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SourceList.cxx,v 1.17 2008/11/09 23:51:37 burnett Exp $
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
#include <time.h>
#include <cstring>

using namespace pointlike;
using astro::SkyDir;

namespace{
    // refit parameters: TODO, provide access
    double group_radius(2.0); // how far to consider for correlated fit
    double too_close(0.25); // too close for an independent source

    double examine_radius(180.), prune_radius(0.25);

    double ts_min(5.0);
    double emin(500);
    int    nside(256);
    double pixel_fraction(1.0);
    astro::SkyDir examine_dir;



    void timer(std::string label="",std::ostream& out = std::cout)
    {
        static bool first=true;
        static time_t start;
        if(first){ first=false; ::time(&start);}
        time_t aclock;
        ::time( &aclock );   
        char tbuf[25]; ::strncpy(tbuf, asctime( localtime( &aclock ) ),24);
        tbuf[24]=0;
        if( !label.empty()) out << label<<std::endl;;
        out<<  "Current time: " << tbuf
            << " ( "<< ::difftime( aclock, start) <<" s elapsed)" << std::endl;
    }

}
double SourceList::set_group_radius(double value)
{
    double old(group_radius);
    group_radius= value;
    return old;
}
const skymaps::BinnedPhotonData * SourceList::s_data(0);
const skymaps::BinnedPhotonData * SourceList::data(){return s_data;}
void SourceList::set_data(const skymaps::BinnedPhotonData * data){s_data=data;}

Source::Source(const std::string& name, const astro::SkyDir& seed_dir, double TS)
: m_name(name)
, m_dir(seed_dir)
, m_seed_dir(seed_dir)
, m_fit(0)
, m_seedTS(TS)
, m_TS(TS)
, m_sigma(0)
, m_neighbor(0)
{
    setup();
}
void Source::setup()
{
    if( SourceList::data()==0){
        throw std::invalid_argument("Source::Source: no data set");
    }
    if(m_fit==0){
        m_fit = new PointSourceLikelihood(*SourceList::data(), m_name, m_dir);
    }
    // inital maximize at seed position, if TS < 10 to start (0 to start)
    // otherwise will sort on the input
    if( m_TS<10.) {
        m_seedTS = m_TS = m_fit->maximize();
    }
}
Source::Source(const std::string& name, double ra, double dec, double TS)
: m_name(name)
, m_dir(SkyDir(ra,dec))
, m_seed_dir(m_dir)
, m_seedTS(TS)
, m_fit(0)
, m_TS(TS)
, m_sigma(0)
, m_neighbor(0)
{
    setup();
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
        <<"name                  ra        dec        TS    localization moved  neighbor\n";
//         B0833-45              128.8339  -45.1629   5698.49    0.0069    0.0136
}

void Source::info(std::ostream& out)const{
        out << std::left << std::fixed 
            << std::setw(20) << name() 
            << std::setprecision(4) << std::right
            << std::setw(10) << dir().ra() 
            << std::setw(10) << dir().dec() 
            << std::setprecision(2) 
            << std::setw(10) << TS() 
            << std::setprecision(4) 
            << std::setw(10) << (TS()<=0? 99: sigma() )
            << std::setprecision(1)
            << std::setw(10) << (TS()<=0? 0. : moved()/sigma() )
                ;
        if( neighbor()!=0){
            out << "  " << std::setw(20)<< std::left << (neighbor()->name());
        }else{
            out << std::setw(10) << " -"; 
        }

        out << std::endl;
}

Source::~Source()
{
    //delete m_fit; //TODO: copy problem with this
}

SourceList::SourceList(std::map<std::string, std::vector<double> > input_list)
: m_verbose(false)
{
    std::map<std::string, std::vector<double> >::const_iterator it(input_list.begin());
    for( ; it!= input_list.end(); ++it ){

        std::string name(it->first);
        const std::vector<double>& v(it->second);
        double TS( v.size()>2? v[2] : 0);
        push_back( Source(name, SkyDir(v.at(0), v.at(1)),TS));
    }

}
SourceList::SourceList(const std::string& filename, bool verbose)
: m_verbose(verbose)
{
    std::ifstream input_file;
    input_file.open(filename.c_str());

    if(!input_file.is_open()) {
        std::cerr << "ERROR:  Unable to open:  " << filename << std::endl;
        throw std::invalid_argument("SourceList: could not open file");
    }
    while (!input_file.eof()){
        std::string line; std::getline(input_file, line);
        if( line.substr(0,5)=="name ") continue; // skip header line
        if( line.size()<5 || line[0]=='#' ) continue; // comment or empty
        // check for a blank in the name
        size_t i(line.find_first_of(" "));
        if( i<5 && line[i+1]!=' '){
            std::cout << "replacing blank in name: " << line << std::endl;
            line[i]='_';
        }
        std::stringstream buf(line); 
        std::string name; buf >> name;
        double ra, dec, TS(0);
        buf >> ra >> dec;
        if( !buf.eof()) buf >> TS;
        push_back( Source(name, SkyDir(ra, dec), TS));
    }
}

class GreaterTS{ // simple functor to sort the list by decreasing TS
public:
    GreaterTS(){}
    bool operator()(const Source& a, const Source& b){
        return a.TS() > b.TS();
    }
};
class Less_ra{ // functor to sort by increasing ra
public:
    Less_ra(){}
    bool operator()(const Source& a, const Source&b){
        return a.dir().ra() < b.dir().ra();
    }
};

void SourceList::sort_TS()
{
    this->sort(  GreaterTS() );
}
void SourceList::sort_ra()
{
    this->sort(  Less_ra() );
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
        if( verbose() ){
            std::cout << ".." << cand.name() << std::endl;
        }
        SkyDir currentpos( cand.dir() );
        double ra(currentpos.ra()), dec(currentpos.dec()); // debug only

        // See if we have already found a candidate near this location
        iterator neighbor;
        double max_value(0.), min_dist(99);
        cand.fit()->clearBackgroundPointSource(); // clear any existing
        cand.set_neighbor(0);

        for( iterator it2 = begin();  it2 != end(); ++it2) {
            Source& other( *it2);
            double diff = currentpos.difference(other.dir())*180/M_PI;
            if( diff <= group_radius &&  other.TS() > cand.TS() ) {
                if( other.TS() > max_value) {
                    max_value = other.TS();
                    neighbor = it2;
                    min_dist=diff;
                }
            }
        }

        if (max_value>0.) {
            // Add strong neighbor to this candidate's background
            cand.fit()->addBackgroundPointSource(neighbor->fit());
            cand.set_neighbor( &*neighbor );
        }
        if( min_dist > too_close){
            // adjust position to maximize likelhood
            cand.localize();
        }else{
            // no, too close to another source, just flag as such
            cand.TS()=-1; 
        }
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
        out << "point("<< cand.dir().ra()<< ","<<cand.dir().dec() <<") # point=cross "
            << size << " text={TS=" 
            << int(cand.TS()+0.5) << "};\n";
        ++n;
    }
    out.close();
    std::cout << "Wrote "<< n << " entries to  reg file "<< filename << ", with TS >= " <<tsmin << std::endl;
}

void SourceList::filter_TS(double threshold)
{

    iterator it(begin());
    int oldsize(size());
    while( it != end()){
        iterator next(it); ++next;
        double TS((*it).TS());
        if( TS < threshold){
            this->erase(it); // it not valid
        }
        it = next;
    }
    int newsize(size());
    std::cout << "applied TS>"<< threshold<<": reduced from: " << oldsize << " to " << newsize << std::endl;
}

#if 0 // under development
void SourceList::examineRegion(void) 
{  
    //
    // ---------------- phase 1: make list of seed positions ---------------------
    //
    double  radius(examine_radius);
    if( radius>=180){
        radius = 179.99999; // bug in gcc version of healpix code
        out() << "Examining full sky"<< std::endl;
    }else{
        out() << "Examining cone of radius "<< radius<<
            " about (ra,dec)=" 
            << std::fixed << std::setprecision(2) 
            << examine_dir.ra()<<", " << examine_dir.dec() 
            << "; (l,b)= "<< examine_dir.l()<<", " << examine_dir.b() 
            << std::resetiosflags(std::ios::fixed)<< std::setprecision(6) 
            << std::endl;
    }
    double emin,emax;
    PointSourceLikelihood::get_energy_range(emin, emax);
    out() << "minimum energy used by PointSourceLikelihood: " 
           << emin << std::endl;


    typedef  std::vector< std::pair<int, int> > PixelVector;
    typedef std::map<int, int> PixelMap;
    PixelMap m;
    skymaps::BinnedPhotonData::const_iterator bpit1 = m_data.begin();
    skymaps::BinnedPhotonData::const_iterator bpit;
    for(; bpit1 != m_data.end(); ++bpit1)
    {
        int tmp(bpit1->nside());
        if(nside != bpit1->nside()) continue;
        bpit = bpit1;
        // load pixels
        PixelVector v;
        bpit->query_disk(examine_dir, radius*M_PI/180, v);
        // add to the map
        for( PixelVector::const_iterator it(v.begin()); it!=v.end(); ++it){
            m[it->first] += it->second;
        }

    }
}
#endif
