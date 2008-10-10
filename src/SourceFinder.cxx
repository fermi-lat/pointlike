/** @file SourceFinder.cxx
@brief implementation of SourceFinder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SourceFinder.cxx,v 1.45 2008/07/28 21:49:15 burnett Exp $
*/

#include "pointlike/SourceFinder.h"
#include "pointlike/PointSourceLikelihood.h"
#include "healpix/HealPixel.h"
#include "healpix/Healpix.h"

#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <cassert>
#include <errno.h>
#include <stdexcept>

namespace {


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

    // Show percent complete
    void ShowPercent(int sofar, int total, int found)
    {
        static int toskip(50), skipped(0);
        if(++skipped<toskip) return; skipped=0;
        static int lastpercent(-1);
        int percent( static_cast<int>(100 * sofar / total +0.5) );
        if( percent==lastpercent) return;
        lastpercent=percent;
        char   s[50];
        sprintf(s, "%d%%, %d found.", percent, found);
        std::cout << s;
        if (sofar < total)
        {
            for (size_t j = 0; j < strlen(s); ++j)
                std::cout << "\b";
        }
        else
            std::cout << std::endl;
    }


    // source finding parameters
    int skip1(2), skip2(3);
    double sigma_max(0.25); // maximum allowable sigma


    double examine_radius(180.), group_radius(2.0), prune_radius(0.25);

    double  ts_min(5.0);
    double emin(500);
    int nside(256);
    double pixel_fraction(1.0);
    astro::SkyDir examine_dir;
    std::string outfile;
    std::string regfile;
    double regtsmin(0);
    std::string regcolor("white");
    std::string fitsfile;
    std::string imagefile;
    std::string logfile;
    double imageresolution(0.1);
    int smooth(1);

    static std::string prefix("SourceFinder.");
} // anon namespace

using namespace pointlike;
using healpix::HealPixel;
using astro::SkyDir;

void SourceFinder::setParameters(const embed_python::Module & module)
{
    module.getValue(prefix+"TSmin", ts_min, ts_min);
    module.getValue(prefix+"emin", emin, emin);
    module.getValue(prefix+"pass1_nside", nside, nside);
    module.getValue(prefix+"pixel_fraction", pixel_fraction, pixel_fraction);

    module.getValue(prefix+"examine_radius", examine_radius);
    module.getValue(prefix+"group_radius", group_radius, group_radius);
    module.getValue(prefix+"prune_radius", prune_radius, prune_radius);

    module.getValue(prefix+"outfile", outfile, "");
    module.getValue(prefix+"regfile", regfile, "");
    module.getValue(prefix+"regtsmin", regtsmin, regtsmin);
    module.getValue(prefix+"regcolor", regcolor, regcolor);
    module.getValue(prefix+"fitsfile", fitsfile, "");
    module.getValue(prefix+"imagefile", imagefile, "");
    module.getValue(prefix+"imageresolution", imageresolution, imageresolution);
    module.getValue(prefix+"smoothed", smooth, smooth);

    module.getValue(prefix+"logfile", logfile, "");
    double l,b,ra,dec;
    module.getValue(prefix+"l", l, 999.);
    module.getValue(prefix+"b", b, 999.);
    module.getValue(prefix+"ra", ra, 999.);
    module.getValue(prefix+"dec", dec, 999.);
    if( l <999 && b < 999) {
        examine_dir = astro::SkyDir(l,b, astro::SkyDir::GALACTIC);
    }else if( ra<999 && dec<999) {
        examine_dir = astro::SkyDir(ra,dec);
    }

}



SourceFinder::SourceFinder(const pointlike::Data& map)
: m_pmap(map)
{
    if( ! logfile.empty() ){
        m_log = new std::ofstream(logfile.c_str());
    }else{
        m_log = & std::cout;
    }
    timer("---------------SourceFinder----------------", out());  

    map.map().info(out());
}
SourceFinder::~SourceFinder()
{
    if( m_log != &std::cout){
        dynamic_cast<std::ofstream*>(m_log)->close();
    }
}

/** @brief
Analyze range of likelihood significance values for all pixels at a particular level  
*/
void SourceFinder::examineRegion(void) 
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
    out() << "minimum energy used by PointSourceLikelihood: " << emin << std::endl;


    typedef  std::vector< std::pair<int, int> > PixelVector;
    typedef std::map<int, int> PixelMap;
    PixelMap m;
#if 0
    skymaps::BinnedPhotonData::const_iterator bpit1 = m_pmap.begin();
    skymaps::BinnedPhotonData::const_iterator bpit;
    for(; bpit1 != m_pmap.end(); ++bpit1)
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
#else // just do all the pixels 
    healpix::Healpix hpx(nside, healpix::Healpix::RING, SkyDir::GALACTIC);
    std::vector<int> v;
    hpx.query_disc(examine_dir, radius*M_PI/180, v);
    for( std::vector<int>::const_iterator it(v.begin()); it!=v.end(); ++it){
        m[*it] = 0; // just set
    }

#endif
    if( m.empty() )
    {
        throw std::invalid_argument("SourceFinder: did not find a Band with requested nside");
    }

    // ------------------- phase 2: examine all pixels with data ---------------------------
    //
    skymaps::SkySpectrum* saved_diffuse = PointSourceLikelihood::set_diffuse(0); // clear the diffuse for this

    out() <<  "First pass will examine " << m.size() 
              << " pixels with nside = " << nside << std::endl;
    Prelim can; // for list of candidates indexed by TS
    can.clear();
    size_t found(0), total(m.size()), count(0);


    // examine un-localized likelihood for each returned pixel
    for(PixelMap::const_iterator it = m.begin(); it != m.end(); ++it, ++count)
    {
        ShowPercent(count, total, found);
#if 0
        SkyDir sd(bpit->dir(it->first));
#else
        healpix::Healpix::Pixel pix(it->first, hpx);

        SkyDir sd( pix() );
#endif
        PointSourceLikelihood ps(m_pmap, "test", sd );
        double ts = ps.maximize();

        if (ts > ts_min)
        {
            can.insert(std::make_pair(ts, CanInfo(ts, 0, sd) ) );
            ++found;
        }
    }
    //
    // --------------------- phase 3: refit the found guys ----------------------------------
    //
    PointSourceLikelihood::set_diffuse(saved_diffuse);  // need full diffuse for this

    size_t fract = static_cast<size_t>(found * pixel_fraction); 

    out() <<  "Found " << found 
        << " pixels with TS >= " << ts_min
        << "; will examine " << fract<< std::endl;

    timer("start refit", out());

    m_can.clear();
    size_t i(0);
    int nbr_skipped(0);


    // Examine full likelihood fit for most significant preliminary candidates
    for (Prelim::reverse_iterator it = can.rbegin(); it != can.rend() && i < fract; ++ it, ++i)
    {
        CanInfo& candidate(it->second);
        ShowPercent(i, fract, m_can.size());

        const SkyDir& currentpos = candidate.dir();
        
        // See if we have already found a candidate near this location
        int neighbor_index(-1);
        double max_value = 0.0;
        bool neighbor_found(false);
        if (group_radius > 0.0)
        {
            for (Candidates::iterator it2 = m_can.begin();  it2 != m_can.end(); ++it2)
            {
                CanInfo& other( it2->second);
                double diff = currentpos.difference(other.dir())*180/M_PI;
                if (diff <= group_radius && other.value() > it->first)
                {
                    if (it2->second.value() > max_value)
                    {
                        max_value = other.value();
                        neighbor_found = true;
                        neighbor_index = it2->first;
                    }
                }
            }
        }
        
        // calculate initial likelihood at current position
        PointSourceLikelihood ps(m_pmap, "test", currentpos);
        // If we found a previous candidate near this location, add it to background
        if (neighbor_found)
        {
            // Recalculate likelihood for strong neighbor
            CanInfo& strong ( m_can[neighbor_index] );
            SkyDir neighbor_dir ( strong.dir() );
            if( strong.fit()==0){
                strong.set_fit(new PointSourceLikelihood(m_pmap, "Strong", neighbor_dir));
                strong.fit()->maximize();
            }

          //  double strong_ts = strong.fit().maximize(); // not used

            // Add strong neighbor to this candidate's background
            ps.addBackgroundPointSource(strong.fit());
            candidate.setStrongNeighbor(neighbor_index);
        }
        double ts = ps.maximize();
        if (ts < ts_min)
        {
            ++ nbr_skipped;
            continue;  // apply initial threshold
        }

        // adjust position to maximize likelhood
        double error = ps.localize(skip1, skip2);
        if (error >= sigma_max) // quit if fail to find maximum or > 1 degree
        {
            ++ nbr_skipped;
            continue;  // apply initial threshold
        }
        ts = ps.maximize(); // readjust likelihood at current position
        if (ts <ts_min)
        {
            ++ nbr_skipped;
            continue;  // apply initial threshold
        }


        // found candidate
        // add to the final list, indexed according to level 13 healpixel
        SkyDir sd = ps.dir();
        int healpix_index (healpix::HealPixel(ps.dir(), 13).index() );
        Candidates::iterator newit = m_can.find(healpix_index);
        if( newit != m_can.end() ) 
        {
            //std::cout << "--duplicate!" << std::endl;
            continue;
        }

        m_can[healpix_index] = CanInfo(ts, error, sd, neighbor_index);


    }
    out() << m_can.size() << " sources found before pruning neighbors.\n";
    timer("", out());
}   


// Eliminate weaker neighbors
// criterion is closer than tolerance, or 3-sigma circles overlap
void SourceFinder::prune_neighbors(void)
{
    double radius(prune_radius);
    if(radius==0) return;
    int pruned_by_distance(0), pruned_by_sigma(0);

    out() << "Eliminating weaker neighbors using radius of " << radius << " degrees...\n";

    // Mark weaker neighbors: loop over all pairs
    for (Candidates::iterator it1 = m_can.begin(); it1 != m_can.end(); ++it1) 
    {
        CanInfo & cand1 ( it1->second);

        double sigma1 ( cand1.sigma());
        for (Candidates::iterator it2 =it1;  it2 != m_can.end(); ++it2)
        {
            if (it1 == it2)   continue;  // Don't compare to yourself.  
            double sigma2(it2->second.sigma());

            double diff = (it1->second.dir().difference(it2->second.dir()))*180/M_PI;
            if (diff <= radius || diff < 3*(sigma1+sigma2) )
            {
                if (diff <= radius)
                    ++ pruned_by_distance;
                if (diff < 3*(sigma1+sigma2))
                    ++ pruned_by_sigma;
                if (sigma2 < sigma1 )it1->second.setDelete();
                else                 it2->second.setDelete();
            }
        }
    }

    // Delete marked entries
    for (Candidates::iterator it = m_can.begin(); it != m_can.end();) 
    {
        Candidates::iterator it2 = it;
        ++ it;
        if (it2->second.is2bdeleted())
            m_can.erase(it2);
    }

    out() << pruned_by_distance << " pruned by distance, " << radius << "degrees\n";
    out() << pruned_by_sigma << " pruned by sigma\n. ";
    out() << m_can.size() << " source Candidates remain.\n";
    timer("", out());
}

void SourceFinder::createReg(const std::string& fileName, double radius, const std::string& color)
{
    std::cout << "Writing results to the reg file " << fileName << std::endl;

    std::ofstream reg;
    reg.open(fileName.c_str());
    reg << "global color="
        << color.c_str()
        << " font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0 width=2;fk5;";
    for (Candidates::const_iterator it = m_can.begin(); it != m_can.end(); ++it)
    {
        int value = static_cast<int>(it->second.value() + 0.5);
        if (radius < -1.0) // Use a cross instead of a circle
        {
            reg << "point("
                << it->second.dir().ra() << ", "
                << it->second.dir().dec();
        }
        else
        {
            reg << "circle("
                << it->second.dir().ra() << ", "
                << it->second.dir().dec() << ", ";
            if (radius > 0.0)
                reg << radius;
            else
                reg << it->second.sigma() * 100.0;
        }
        reg << ") # point=cross " 
            << " text = {" << value << "};";
    }
    reg.close();
}

static float precision(double x, double s) // s is power of 10, like 0.01 
{
    return float(int(x/s+0.5))*s;
}
void SourceFinder::createTable(const std::string& fileName)
{
    std::cout << "Writing results to the table " << fileName << std::endl;

    std::ofstream table;
    table.open(fileName.c_str());
    std::string delim("\t");
    table << "id\tra\tdec\tTS\tsigma\tneighbor ";

    table << std::endl;


    for (Candidates::iterator it = m_can.begin(); it != m_can.end(); ++it) {
        CanInfo & cand ( it->second);
        SkyDir dir = cand.dir();

        table
            << it->first << delim
            << std::fixed << std::setprecision(4) 
            << dir.ra()  << delim
            << dir.dec() << delim
            << cand.value() << delim  // TS
            << cand.sigma() << delim // error
            << std::resetiosflags(std::ios::fixed)<< std::setprecision(6) 
            << cand.strongNeighbor() 
            << std::endl;
            ;
    }
    table.close();
    timer("",out());
}

std::vector<CanInfo> SourceFinder::candidateList()const
{
    std::vector<CanInfo> ret;
    for (Candidates::const_iterator it = m_can.begin(); it != m_can.end(); ++it) {
        if( !it->second.is2bdeleted() ){
            ret.push_back(it->second);
        }
    }
    return ret;
}

void SourceFinder::createFitsFile(const std::string & outputFile,
                         const std::string & tablename,
                         bool clobber) const
{
    std::cout << "Writing results to the FITS file, table " 
        << outputFile <<", " << tablename << std::endl;

    if (clobber)
    {
        int rc = std::remove(outputFile.c_str());
        if( rc == -1 && errno == EACCES ) 
            throw std::runtime_error(std::string(" Cannot remove file " + outputFile));
    }

    // now add a table to the file
    tip::IFileSvc::instance().appendTable(outputFile, tablename);
    tip::Table & table = *tip::IFileSvc::instance().editTable( outputFile, tablename);
    table.appendField("NAME", "24A");
    table.appendField("TYPE", "30A");
    table.appendField("RA", "1E");
    table.appendField("DEC", "1E");
    table.appendField("L", "1E");
    table.appendField("B", "1E");
    table.appendField("F100", "1E");
    table.appendField("EMIN", "1E");
    table.appendField("EMAX", "1E");
    table.appendField("G1", "1E");
    table.appendField("G2", "1E");
    table.appendField("EB", "1E");
    table.appendField("SRCID", "1J");
    table.appendField("Z", "1E");
    table.appendField("NEVT", "1J");
    table.appendField("TO", "1E");
    table.appendField("CLOSEST", "1J");
    table.appendField("DISTIN", "1E");
    table.appendField("CONFUSED", "1I");
    table.appendField("TS", "1E");
    table.appendField("ERRX", "1E");
    table.appendField("NPRED", "1E");
    table.setNumRecords(m_can.size());

    // get iterator for the Table 
    tip::Table::Iterator itor = table.begin();

    // now just copy
    for (Candidates::const_iterator it = m_can.begin(); it != m_can.end(); ++itor, ++it)
    {
        (*itor)["NAME"].set(" ");
        (*itor)["TYPE"].set(" ");
        (*itor)["RA"].set(it->second.dir().ra());
        (*itor)["DEC"].set(it->second.dir().dec());
        (*itor)["L"].set(it->second.dir().l());
        (*itor)["B"].set(it->second.dir().b());
        (*itor)["F100"].set(0.0);
#ifdef OLD
        (*itor)["EMIN"].set(m_pmap.energyBins().front());
        (*itor)["EMAX"].set(m_pmap.energyBins().back());
#endif
        (*itor)["G1"].set(0.0);
        (*itor)["G2"].set(0.0);
        (*itor)["EB"].set(0.0);
        (*itor)["SRCID"].set(0);
        (*itor)["Z"].set(0.0);
        (*itor)["NEVT"].set(0.0);
        (*itor)["TO"].set(0.0);
        (*itor)["CLOSEST"].set(0.0);
        (*itor)["DISTIN"].set(0.0);
        (*itor)["CONFUSED"].set(0);
        (*itor)["TS"].set(it->second.value());
        (*itor)["ERRX"].set(it->second.sigma());
        (*itor)["NPRED"].set(0.0);
    }

    // set the headers (TODO: do the comments, too)
    tip::Header& hdr = table.getHeader();
    hdr["NAXIS1"].set((20 * sizeof(long)) + 50);

    // close it?
    delete &table;
}
void SourceFinder::createRegFile(std::string filename, std::string color, double tsmin)const
{
    std::ofstream out(filename.c_str());
    out << "global color="<< color 
        << " font=\"helvetica 10 normal\" select=1 edit=1 move=0 delete=1 include=1 fixed=0 width=1;fk5;"
        << std::fixed << std::setprecision(4) << std::endl;
    int n(0);
    for( Candidates::const_iterator it = m_can.begin(); it != m_can.end();  ++it)  {
        const CanInfo& cand = it->second;
        if(cand.value()< tsmin) continue;
        out << "point("<< cand.ra()<< ","<<cand.dec() <<") # point=cross 20 text={TS=" 
            << int(cand.value()+0.5) << "};\n";
        ++n;
    }
    out.close();
    std::cout << "Wrote "<< n << " entries to  reg file "<< filename << ", with TS >= " <<tsmin << std::endl;
}


void SourceFinder::run()
{

    // draw the data region if requested
    if( !imagefile.empty()) {
        Draw drawer(m_pmap);
        double fov(examine_radius);
        drawer.region(examine_dir, imagefile, imageresolution, fov, smooth);
    }
    
    // do the work
    examineRegion();

    //  prune 
    prune_neighbors();   


    // and write out the ascii table, reg or fits files
    if( ! outfile.empty() )  createTable(outfile);   

    if( ! regfile.empty() ) createRegFile(regfile, regcolor, regtsmin);

    if( ! fitsfile.empty() ) createFitsFile(fitsfile);
}


