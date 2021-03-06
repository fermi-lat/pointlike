/** @file SourceList.cxx
@brief implementation of classes Source and SourceList

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SourceList.cxx,v 1.23 2009/06/03 22:18:09 burnett Exp $
*/


#include "pointlike/SourceList.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/LeastSquaresFitter.h"
#include "skymaps/BinnedPhotonData.h"
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <time.h>
#include <cstring>

#define USELSQ

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

    std::ostream* log_stream(&std::cout); // log steam to write to
    std::ostream& logger(){return *log_stream;} 

    void timer(std::string label="")
    {
        static bool first=true;
        static time_t start;
        if(first){ first=false; ::time(&start);}
        time_t aclock;
        ::time( &aclock );   
        char tbuf[25]; ::strncpy(tbuf, asctime( localtime( &aclock ) ),24);
        tbuf[24]=0;
        if( !label.empty()) logger() << label<<std::endl;;
        logger()<<  "Current time: " << tbuf
            << " ( "<< ::difftime( aclock, start) <<" s elapsed)" << std::endl;
    }

}
void SourceList::set_log(std::ostream* logstream){log_stream=logstream;}
double SourceList::set_group_radius(double value)
{
    double old(group_radius);
    group_radius= value;
    return old;
}
const skymaps::BinnedPhotonData * SourceList::s_data(0);
const skymaps::BinnedPhotonData * SourceList::data(){return s_data;}
void SourceList::set_data(const skymaps::BinnedPhotonData * data){s_data=data;}
bool SourceList::s_uselsq(false);
bool SourceList::set_uselsq(bool q){
    bool t(q); s_uselsq=q; return t;
}
bool SourceList::uselsq(){return s_uselsq;}

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
    m_fit->set_ostream(logger());
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
    m_sigma = m_fit->errorCircle(); // estimate of the sigma
    if( m_sigma>1){
        m_sigma=100;
        return m_sigma; // flag no data.
    }
    SkyDir idir( m_fit->dir() );
    double ra(idir.ra()), dec(idir.dec());
    // standard
    double sigma = m_fit->localize() ; // standard gradient fit
    if( sigma<1){ 
        m_sigma=sigma; // keep if OK.
        m_TS = m_fit->TS();
        m_dir = m_fit->dir();

    }else{
        // failed: restore dir in case localize did not
        m_fit->setDir(idir);
    }
    if( SourceList::uselsq()){ // use lsq here
        double iTS = m_fit->TSmap(idir);
        // try this
        LeastSquaresFitter lsq(*m_fit, m_sigma);
        sigma = lsq.err();
        m_fitparams = lsq.ellipse();
    }else{
        //std::cout << "Note: using radius of 1 sigma" << std::endl;
        m_fitparams = TScircle(3.0); // just dump the values after gradient fit
    }
    m_fit->printSpectrum();
    return m_sigma;
}

std::vector<double> Source::TScircle(double radius)const
{
    static double d(1/sqrt(2.) );
    static double x[]= { 1,  d,  0, -d, -1, -d,  0,  d};
    static double y[]= { 0,  d,  1,  d,  0, -d, -1, -d};
    double ra( m_fit->dir().ra()), dec(m_fit->dir().dec());
    double decsig( radius*m_sigma ); // how far to go in sigma units
    double rasig( decsig/cos(dec*M_PI/180) ); // scale ra accordingly (fails if very near a pole.)

    std::vector<double> ret;
    double zTS = m_fit->TSmap(SkyDir(ra,dec)); // value at center
    for(int i(0); i< 8; ++i){
        double dx(x[i]*rasig), dy(y[i]*decsig);
        double ts(m_fit->TSmap( SkyDir( ra+dx, dec+dy) ));
        //std::cout << (ra+dx) << ", " << (dec+dy) << ", " << ts << ", " << (zTS-ts) << std::endl;
        ret.push_back(zTS- ts);
    }
    return ret;
}

double Source::moved()const{
    if(m_sigma>1 || m_sigma==0) return 0; // failed, or not localized
    return dir().difference(seed_dir())*180/M_PI;
}

void Source::header(std::ostream& out){
    out << std::left << std::setw(20) 
        <<"#name                  ra        dec        TS    localization   moved  neighbor"; 
        
//         B0833-45              128.8339  -45.1629   5698.49    0.0069    0.0136
    if( SourceList::uselsq() ){
        out << LeastSquaresFitter::header();
    }else{
        out << "\n";
    }
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
            << std::setw(10) << (TS()<=0? 0. : moved()/sigma() 
            )
                ;
        if( neighbor()!=0){
            out << "  " << std::setw(20)<< std::left << (neighbor()->name());
        }else{
            out << std::setw(10) << " -"; 
        }

        if( SourceList::uselsq() ){
            if(! m_fitparams.empty()){
                out << std::setprecision(4);
                for(std::vector<double>::const_iterator it = m_fitparams.begin(); it!= m_fitparams.end(); ++it){
                    out << std::setw(12) << (*it);
                }
            }else{
                for(int i(0); i<5; ++i){ // length wired in :-(
                    out << std::setw(12) << "0";
                }
                out << std::setw(12) << "100"; //last is quality 
            }
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
SourceList::SourceList(const std::string& filename, bool verbose,std::ostream* log)
: m_verbose(verbose)
{
    timer("SourceList start:");
    std::ifstream input_file;
    input_file.open(filename.c_str());

    if(!input_file.is_open()) {
        logger() << "ERROR:  Unable to open:  " << filename << std::endl;
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
            logger() << "replacing blank in name: " << line << std::endl;
            line[i]='_';
        }
        std::stringstream buf(line); 
        std::string name; buf >> name;
        double ra, dec, TS(0);
        buf >> ra >> dec;
        if( !buf.eof()) buf >> TS;
        push_back( Source(name, SkyDir(ra, dec), TS));
        logger() << "---- " << name << ", initial TS: " << back().TS() << std::endl;
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
    if( outfilename.find(".xml") == std::string::npos ) {
        dump( out );
    }else{
        dump_xml(out, outfilename);
    }
}
void SourceList::refit()
{
    timer("SourceList::refit begin");
    for( iterator it= begin(); it!= end(); ++it){
        Source& cand(*it);
        if( verbose() ){
            logger() << "---------" << cand.name() << std::endl;
        }
        SkyDir currentpos( cand.dir() );
        //double ra(currentpos.ra()), dec(currentpos.dec()); // debug only

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
            double sig = cand.localize();
            if( sig>1){
                logger() << "---Failed localization" << std::endl;
            }
        }else{
            // no, too close to another source, just flag as such
            cand.TS()=-1; 
        }
    }
    timer("SourceList::refit end");;
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
    logger()<< "Wrote "<< n << " entries to  reg file "<< filename << ", with TS >= " <<tsmin << std::endl;
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
    logger() << "applied TS>"<< threshold<<": reduced from: " << oldsize << " to " << newsize << std::endl;
}

void SourceList::dump_xml(std::ostream& out, std::string name)const
{
    out <<"<?xml version='1.0'?> \n"
"<VOTABLE version=\"1.1\"\n"
" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
" xsi:schemaLocation=\"http://www.ivoa.net/xml/VOTable/v1.1\n" 
" http://www.ivoa.net/xml/VOTable/v1.1\" \n"
" xmlns=\"http://www.ivoa.net/xml/VOTable/v1.1\"> \n"
" <RESOURCE>\n"
" <TABLE name=\"" << name << "\" nrows=\"" << size() << "\">\n"
" <FIELD arraysize=\"15\" datatype=\"char\" name=\"name\"/> \n"
" <FIELD datatype=\"float\" name=\"ra\" ucd=\"pos.eq.ra;meta.main\"/>\n"
" <FIELD datatype=\"float\" name=\"dec\" ucd=\"pos.eq.dec;meta.main\"/>\n"
" <FIELD datatype=\"float\" name=\"TS\"/>\n"
" <FIELD datatype=\"float\" name=\"localization\"/>\n"
" <FIELD datatype=\"float\" name=\"moved\"/>\n"
" <FIELD arraysize=\"15\" datatype=\"char\" name=\"neighbor\"/>\n"
" <DATA>\n"
"<TABLEDATA>\n";
    for( const_iterator it = begin(); it != end();  ++it)  {
        const Source& cand( * it);
        out <<"<TR>\n <TD>" << cand.name() << "</TD>\n"
            << " <TD>" << cand.dir().ra() << "</TD>\n"
            << " <TD>" << cand.dir().dec() << "</TD>\n"
            << " <TD>" << cand.TS() << "</TD>\n"
            << " <TD>" << cand.sigma() << "</TD>\n"
            << " <TD>" << (cand.TS()<0? 0 : cand.moved()/cand.sigma()) << "</TD>\n"
            << " <TD>" << (cand.neighbor()!=0?cand.neighbor()->name(): "-") << "</TD>\n"
            << "</TR>\n";
    }
    out << "</TABLEDATA>\n</DATA>\n</TABLE>\n</RESOURCE>\n</VOTABLE>" << std::endl;

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
