/** @file SourceFinder.cxx
@brief implementation of SourceFinder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SourceFinder.cxx,v 1.31 2008/03/22 17:43:55 burnett Exp $
*/

#include "pointlike/SourceFinder.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/PowerLawFilter.h"
#include "tip/IFileSvc.h"
#include "tip/Table.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <cassert>
#include <errno.h>
#include <stdexcept>

namespace {


    void timer(std::string label="",std::ostream& out = std::cout)
    {
        time_t t1;
        time(&t1);
        out << "\n"<< ctime(&t1);
        if( !label.empty()) out << label<<std::endl;;
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


    double examine_radius, group_radius, prune_radius;

    double  ts_min;
    int  pix_level;
    int  skip_TS_levels;
    double pixel_fraction;
    astro::SkyDir examine_dir;
    std::string outfile;
    std::string regfile;
    std::string fitsfile;
    int refit(1);

    static std::string prefix("SourceFinder.");

    void getParameters(const embed_python::Module & module)
    {
        module.getValue(prefix+"TSmin", ts_min, 10);
        module.getValue(prefix+"skipTSlevels", skip_TS_levels, 0);
        module.getValue(prefix+"pixLevel", pix_level, 8);
        module.getValue(prefix+"pixel_fraction", pixel_fraction, 0.5);

        module.getValue(prefix+"examine_radius", examine_radius);
        module.getValue(prefix+"group_radius", group_radius, 1.0);
        module.getValue(prefix+"prune_radius", prune_radius, 0.25);
        module.getValue(prefix+"refit", refit, 1);

        module.getValue(prefix+"outfile", outfile, "");
        module.getValue(prefix+"regfile", regfile, "");
        module.getValue(prefix+"fitsfile", fitsfile, "");

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

        std::cout << "\nSourceFinder parameters:\n"
            << "  Likelihood minimum: " << ts_min<< "\n" 
            << "  skip when localizing: " << skip1<< " to " << skip2  <<"\n"
            << "  sigma_max: " << sigma_max << "\n"
            << "  skip TS levels: " << skip_TS_levels << "\n"
            << std::endl;


    }

} // anon namespace

using namespace pointlike;
using namespace embed_python;
using healpix::HealPixel;
using astro::SkyDir;

SourceFinder::SourceFinder(const pointlike::Data& map, const Module & py_module)
: m_pmap(map)
, m_counts(0)
, m_module(py_module)
{
    getParameters(py_module);
}


/** @brief
Analyze range of likelihood significance values for all pixels at a particular level  
*/
void SourceFinder::examineRegion(void) 
{  
    timer("---------------SourceFinder::examineRegion----------------");  
    getParameters(m_module); // load all parameters

    double  radius(examine_radius);
    if( radius>=180){
        radius = 179.99999; // bug in gcc version of healpix code
        std::cout << "Examining full sky"<< std::endl;
    }else{
        std::cout << "Examining cone of radius "<< radius<<
            " about (ra,dec)=" 
            << std::setprecision(5) << examine_dir.ra()<<", " << examine_dir.dec() << std::endl;
    }


    std::vector<std::pair<healpix::HealPixel, int> > v;

    // Extract the pixels to be examined
    m_pmap.extract_level(examine_dir, radius, v, pix_level, true);
    Prelim can;
    can.clear();
    size_t num(0);


    // examine all non-zero pixels
    for(std::vector<std::pair<healpix::HealPixel, int> >::const_iterator it = v.begin();
        it != v.end(); ++it)
    {
        astro::SkyDir sd;
        int count = static_cast<int>(m_pmap.photonCount(it->first, sd));

        if (count > 0)
        {
            can.insert(std::pair<int, CanInfo>(count, CanInfo(count, 0, sd)) );
            ++num;
        }
    }

    size_t fract = static_cast<size_t>(num*pixel_fraction); 

    std::cout <<  "Found " << num << " non-empty pixels at level " << pix_level 
        << "; will examine " << fract<< std::endl;

    m_can.clear();
    size_t i(0);
    int nbr_skipped(0);


    // Examine likelihood fit for each preliminary candidate
    for (Prelim::reverse_iterator it = can.rbegin(); it != can.rend() && i<fract; ++ it, ++i)
    {
        ShowPercent(i, can.size(), m_can.size());
        int check = it->first;

        const astro::SkyDir& currentpos = it->second.dir();

        // perform likelihood analysis at the current candidate position  
        PointSourceLikelihood ps(m_pmap, "test", currentpos );

        double ts = ps.maximize(skip_TS_levels);
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
        ts = ps.maximize(skip_TS_levels); // readjust likelihood at current position
        if (ts <ts_min)
        {
            ++ nbr_skipped;
            continue;  // apply initial threshold
        }

        // found candidate
        // add to the final list, indexed according to level 13 location
        HealPixel px(ps.dir(), 13); 
        //std::cout << "Found candidate: pixel, ts:" << px.index() <<", " << ts << std::endl;
        Candidates::iterator newit = m_can.find(px);
        if( newit != m_can.end() ) {
            //std::cout << "--duplicate!" << std::endl;
            continue;
        }
        m_can[px] = CanInfo(ts, error, ps.dir());
        for(int id =ps.minlevel() ;id<=ps.maxlevel();++id)
        {
            double levTS = ps.levelTS(id);
            if( levTS<-1 ) {
                std::cout << "   Warning: negative TS for level:" << id << ", " << levTS << std::endl; 
            }
            m_can[px].setValue(id,ps.levelTS(id));
            m_can[px].setPhotons(id,(ps[id]->photons()) * (ps[id]->alpha()));
            m_can[px].setSigalph(id,ps[id]->sigma_alpha());
        }

        // Calculate and store power law fit values.  New as of 6/5/07
        std::vector<std::pair<double, double> > values;
        values.clear();
        std::vector<double> energyBins = m_pmap.energyBins();
        // ignore first two levels
        for( int i = 2, lvl = m_pmap.minLevel() + 2; lvl < m_pmap.minLevel() + m_pmap.levels(); ++i, ++lvl)
        {
            double count = m_can[px].photons(lvl);
            values.push_back(std::make_pair(energyBins[i], count));
        }
        pointlike::PowerLawFilter pl(values);
        m_can[px].set_pl_slope(pl.slope());
        m_can[px].set_pl_constant(pl.constant());
        m_can[px].set_pl_confidence(pl.metric());
        m_can[px].set_weighted_count(check);
    }
    std::cout << m_can.size() << " sources found before pruning neighbors.\n";
    timer();
}       
void SourceFinder::reExamine(void) 
{  
    timer("---------------SourceFinder::reExamine----------------");  
    int i(0), nbr_to_examine(m_can.size()), nbr_purged(0);

    // Re-examine likelihood fit for each candidate that has a strong neighbor.
    // Add the strong neighbor to the background first
    for (Candidates::iterator it = m_can.begin(); it != m_can.end(); ) 
    {
        // 2nd iterator makes it possible to delete and still iterate with the other one
        Candidates::iterator it2 = it;  
        ++ it;

        ShowPercent(i++, nbr_to_examine, nbr_purged);
        if (! (it2->second.hasStrongNeighbor())) continue;  // Only care about ones with strong neighbors

        CanInfo& cand = it2->second;
        double oldts(cand.value());

        //not used? const astro::SkyDir& currentpos = it2->second.dir();
        PointSourceLikelihood::clearBackgroundPointSource();

        // Recalculate likelihood for strong neighbor
        PointSourceLikelihood strong(m_pmap, "test", m_can[it2->second.strongNeighbor()].dir());
        double strong_ts = strong.maximize(skip_TS_levels); // not used

        // Add strong neighbor to this candidate's background
        PointSourceLikelihood::addBackgroundPointSource(& strong);

        // Recalculate likelihood for this candidate
        PointSourceLikelihood ps(m_pmap, "test", it2->second.dir());
        double ts = ps.maximize(skip_TS_levels);

        // perform likelihood analysis at the current candidate position 
        //already done ts = ps.maximize(skip_TS_levels);

        // eliminate candidate if now below threshold
        if (ts < ts_min)
        {
            ++ nbr_purged;
            m_can.erase(it2);
            continue; 
        }

        // adjust position to maximize likelhood
        double error = ps.localize(skip1, skip2);
        if (error >= sigma_max) // quit if fail to find maximum or > 1 degree
        {
            ++ nbr_purged;
            m_can.erase(it2);
            continue;  
        }

        ts = ps.maximize(skip_TS_levels); // readjust likelihood at current position
        if (ts < ts_min)
        {
            ++ nbr_purged;
            m_can.erase(it2);
            continue;  // apply initial threshold
        }

        // candidate is still good!  update CanInfo values
        it2->second.set_total_value(ts);
        it2->second.set_sigma(error);
        it2->second.set_dir(ps.dir());
        for(int id =ps.minlevel(); id <= ps.maxlevel(); ++id)
        {
            it2->second.setValue(id,ps.levelTS(id));
            it2->second.setPhotons(id,(ps[id]->photons()) * (ps[id]->alpha()));
            it2->second.setSigalph(id,ps[id]->sigma_alpha());
        }

        // Calculate and store power law fit values.  
        std::vector<std::pair<double, double> > values;
        values.clear();
        std::vector<double> energyBins = m_pmap.energyBins();
        // ignore first two levels
        for( int i = 2, lvl = m_pmap.minLevel() + 2; lvl < m_pmap.minLevel() + m_pmap.levels(); ++i, ++lvl)
        {
            double count = it2->second.photons(lvl);
            values.push_back(std::make_pair(energyBins[i], count));
        }
        pointlike::PowerLawFilter pl(values);
        it2->second.set_pl_slope(pl.slope());
        it2->second.set_pl_constant(pl.constant());
        it2->second.set_pl_confidence(pl.metric());
    }
    std::cout << nbr_purged << " candidates purged,  " << m_can.size() << " candidates left.\n";
    timer();
}     
void SourceFinder::checkDir(astro::SkyDir & sd,
                            double eq_TS_min,
                            double mid_TS_min,
                            double polar_TS_min,
                            int    pix_level, 
                            int count_threshold,
                            bool   /*background_filter*/,
                            int	skip_TS_levels,
                            double equator_boundary,
                            double polar_boundary) 
{
    int skip1(2), skip2(3),
        photon_count_check(2);
    double sigma_max(0.25); // maximum allowable sigma

    healpix::HealPixel hp(sd, pix_level);
    astro::SkyDir new_dir;
    double count = m_pmap.photonCount(hp, new_dir);
    std::cout << "\nChecking for a source at (ra, dec) = ("
        << sd.ra() << ", " << sd.dec() <<
        ") (l, b) = (" << sd.l() << ", " << sd.b() <<
        ") pixel level: " << pix_level << std::endl;
    std::cout << "HealPixel direction = ("
        << hp().ra() << ", " << hp().dec() <<
        ") (" << hp().l() << ", " << hp().b() <<
        ") pixel index: " << hp.index() << std::endl;
    std::cout << "Weighted count: " << count << std::endl;
    std::cout << "Weighted direction = ("
        << new_dir.ra() << ", " << new_dir.dec() <<
        ") ("  << new_dir.l() << ", " << new_dir.b() << ")\n";
    if (count < count_threshold)
        std::cout << "  ** count is less than threshold of " << count_threshold << std::endl;

    double ts_min, abs_b = fabs(new_dir.b());
    if (abs_b < equator_boundary)
        ts_min = eq_TS_min;
    else if (abs_b > polar_boundary)
        ts_min = polar_TS_min;
    else
        ts_min = mid_TS_min;
    PointSourceLikelihood ps(m_pmap, "test", new_dir);
    ps.set_verbose(true);

    double ts = ps.maximize(skip_TS_levels);
    ps.printSpectrum();
    std::cout << "Initial ts: " << ts << std::endl;
    if (ts < ts_min)
        std::cout << "  ** ts is less than minimum of " << ts_min << std::endl;

    double error = ps.localize(skip1, skip2);
    if (error >= sigma_max)
    {
        std::cout << "  ** No max found. " << error << " returned from localize.\n";
        return;
    }
    ts = ps.maximize(skip_TS_levels); // readjust likelihood at current position
    ps.printSpectrum();
    std::cout << "Final ts: " << ts << std::endl;
    if (ts < ts_min)
        std::cout << "  ** ts is less than threshold of " << ts_min << std::endl;
    std::cout << "Maximized direction = ("
        << ps.dir().ra() << ", " << ps.dir().dec() <<
        ") (" << ps.dir().l() << ", " << ps.dir().b() << ")\n";

    HealPixel px_check(ps.dir(), 8);
    count = m_pmap.photonCount(px_check, true, false);
    if (count < photon_count_check)
        std::cout << "  ** Parent level 8 pixel photon count < " << photon_count_check << std::endl;
    std::cout << "sigma = " << ps.errorCircle() << std::endl;

    HealPixel px_final(ps.dir(), 13);
    std::cout << "Final level 13 pixel direction = ("
        << px_final().ra() << ", " << px_final().dec() <<
        ") (" << px_final().l() << ", " << px_final().b() << ")\n";
    double distance = sd.difference(px_final()) * (180 / M_PI);
    std::cout << "Distance from starting direction: " << distance << " degrees.  distance/sigma = " << distance/ps.errorCircle() << std::endl;
}





// List selected pixels
void SourceFinder::list_pixels()
{
    std::cout << "\nl \t b \t ra \t dec \t level \t index \t count \n";
    for (Candidates::const_iterator it = m_can.begin(); it != m_can.end(); ++it)
    {
        std::cout << it->first().l() << "\t" 
            << it->first().b() << "\t"
            << it->first().ra() << "\t"
            << it->first().dec() << "\t"
            << it->first.level() << "\t"
            << it->first.index() << "\t"
            << it->second.value() << "\n";
    }
}
// Eliminate candidates that don't meet power law criteria
void SourceFinder::prune_power_law(void)
{
    // Get parameters 
    double  plEqTSmin;
    m_module.getValue(prefix+"plEqTSmin", plEqTSmin, 38.0);

    double  plMidTSmin;
    m_module.getValue(prefix+"plMidTSmin", plMidTSmin, 19.0);

    double  plPolarTSmin;
    m_module.getValue(prefix+"plPolarTSmin", plPolarTSmin, 24.0);

    double  plSlopeCutoff;
    m_module.getValue(prefix+"plSlopeCutoff", plSlopeCutoff, -1.5);

    double  plFitCutoff;
    m_module.getValue(prefix+"plFitCutoff", plFitCutoff, 0.9);

    double  eqBoundary;
    m_module.getValue(prefix+"eqBoundary", eqBoundary, 6.0);

    double  polarBoundary;
    m_module.getValue(prefix+"polarBoundary", polarBoundary, 6.0);

    int minlevel = m_pmap.minLevel() + 2;  // Ignore lowest 2 TS levels

    std::cout << "Eliminating candidates by power law test:\n" 
        << "  Slope cutoff: " << plSlopeCutoff << std::endl
        << "  Confidence cutoff: " << plFitCutoff << std::endl
        << "  Ignore power law fit if total TS for levels >= " << minlevel << std::endl
        << "    >= " << plEqTSmin << " in equatorial region" << std::endl
        << "    >= " << plMidTSmin << " in middle region" << std::endl
        << "    >= " << plPolarTSmin << " in polar region" << std::endl;

    for (Candidates::iterator it1 = m_can.begin(); it1 != m_can.end();) 
    {
        // 2nd iterator makes it possible to delete and still iterate with the other one
        Candidates::iterator it2 = it1;  
        ++ it1;

        double totalTS = 0;
        for (int i = minlevel; i < m_pmap.minLevel() + m_pmap.levels(); ++i)
            totalTS += it2->second.values(i);

        double abs_b = fabs(it2->first().b());


        // Skip power law test for candidates with big enough TS
        if (abs_b < eqBoundary)
        {
            if (totalTS >= plEqTSmin) continue;
        }
        else if (abs_b > polarBoundary)
        {
            if (totalTS >= plPolarTSmin) continue;
        }
        else if (totalTS >= plMidTSmin) continue;

        if (it2->second.pl_slope() > plSlopeCutoff 
            || it2->second.pl_confidence() < plFitCutoff)
            m_can.erase(it2);
    }

    std::cout << m_can.size() << " source Candidates remain.\n";
    timer();
}

// Eliminate weaker neighbors
// criterion is closer than tolerance, or 3-sigma circles overlap
void SourceFinder::prune_neighbors(void)
{
    double radius(prune_radius);
    if(radius==0) return;
    int pruned_by_distance(0), pruned_by_sigma(0);

    std::cout << "Eliminating weaker neighbors using radius of " << radius << " degrees...\n";

    // Mark weaker neighbors: loop over all pairs
    for (Candidates::iterator it1 = m_can.begin(); it1 != m_can.end(); ++it1) 
    {
        //        if( it1->second.is2bdeleted() ) continue; // if this one aleady deleted, ignore
        double sigma1 ( it1->second.sigma());
        for (Candidates::iterator it2 =it1;  it2 != m_can.end(); ++it2)
        {
            if (it1 == it2)   continue;  // Don't compare to yourself.  
            //            if( it2->second.is2bdeleted() ) continue; // already gone
            double sigma2(it2->second.sigma());

            double diff = (it1->first)().difference((it2->first)())*180/M_PI;
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

    std::cout << pruned_by_distance << " pruned by distance.\n";
    std::cout << pruned_by_sigma << " pruned by sigma.\n";
    std::cout << m_can.size() << " source Candidates remain.\n";
    timer();
}

// Eliminate weaker neighbors
void SourceFinder::prune_adjacent_neighbors()
{
    std::cout << "Eliminating weaker neighbors...";

    // Mark weaker neighbors.
    for (Candidates::iterator it = m_can.begin(); it != m_can.end(); ++it) 
    {
        std::vector<healpix::HealPixel> hv = it->first.neighbors();
        for (std::vector<healpix::HealPixel>::const_iterator n = hv.begin();
            n != hv.end(); ++n)
        {
            Candidates::const_iterator p = m_can.find(*n);
            if (p != m_can.end() && it->second.value() < p->second.value())
            {
                it->second.setDelete();
                break;
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

    std::cout << m_can.size() << " source Candidates remain.\n";
    timer();
}

// Group nearby candidates to facilitate further examination
// Any candidate that is within "group_radius" of another is matched with its strongest neighbor
void SourceFinder::group_neighbors(void)
{
    double radius(group_radius); // use radius from file-scope varialbe

    int nbr_found = 0;
    if(radius==0) return;

    std::cout << "Grouping neighbors using radius of " << radius << " degrees...";

    // Mark candidates with stronger neighbors: loop over all pairs
    for (Candidates::iterator it1 = m_can.begin(); it1 != m_can.end(); ++it1) 
    {
        double max_value = 0.0;
        for (Candidates::iterator it2 = m_can.begin();  it2 != m_can.end(); ++it2)
        {
            if (it1 == it2)   continue;  // Don't compare to yourself.  

            double diff = (it1->first)().difference((it2->first)())*180/M_PI;
            if (diff <= radius && it2->second.value() > it1->second.value())
            {
                if (it2->second.value() > max_value)
                {
                    max_value = it2->second.value();
                    it1->second.setHasStrongNeighbor(true);
                    it1->second.setStrongNeighbor(it2->first);
                }
            }
        }
        if (max_value > 0.0)
            ++ nbr_found;
    }

    // Follow chains.  If A points to B and B points to C, have A point to C.
    for (Candidates::iterator it1 = m_can.begin(); it1 != m_can.end(); ++it1)
    {
        if (! (it1->second.hasStrongNeighbor())) continue; // Only look at ones with stronger neighbor.

        healpix::HealPixel px = it1->second.strongNeighbor();
        while (m_can[px].hasStrongNeighbor())
        {
            px = m_can[px].strongNeighbor();
            it1->second.setStrongNeighbor(px);
        }
    }

    std::cout << nbr_found << " candidates have stronger neighbors.\n";
    timer();
}

void SourceFinder::createReg(const std::string& fileName, double radius, const std::string& color)
{
    std::cout << "Writing results to the reg file " << fileName << std::endl;

    std::ofstream reg;
    reg.open(fileName.c_str());
    reg << "global color="
        << color.c_str()
        << " font=\"helvetica 10 normal\" select=1 edit=1 move=1 delete=1 include=1 fixed=0 width=2;";
    for (Candidates::const_iterator it = m_can.begin(); it != m_can.end(); ++it)
    {
        int value = static_cast<int>(it->second.value() + 0.5);
        if (radius < -1.0) // Use a cross instead of a circle
        {
            reg << "fk5; cross point("
                << (it->first)().ra() << ", "
                << (it->first)().dec();
        }
        else
        {
            reg << "fk5; circle("
                << it->second.dir().ra() << ", "
                << it->second.dir().dec() << ", ";
            if (radius > 0.0)
                reg << radius;
            else
                reg << it->second.sigma() * 100.0;
        }
        reg << ") " << "# text = {" << value << "};";
    }
    reg.close();
}

static float precision(double x, double s) // s is power of 10, like 0.01 
{
    return float(int(x/s+0.5))*s;
}
void SourceFinder::createTable(const std::string& fileName, 
                               bool/* get_background*/,   // not used?
                               int /*skip_TS*/)
{
    std::cout << "Writing results to the table " << fileName << std::endl;

    std::ofstream table;
    table.open(fileName.c_str());
    std::string delim=(" ");
    table << "Name  source_flag ra  dec  sigma  circle  l b ";
    table << "TS6 TS7 TS8 TS9 TS10 TS11 TS12 TS13 ";
    table << "n6 n7 n8 n9 n10 n11 n12 n13 ";
    table << "pl_m pl_b pl_sigma pl_chi_sq ";
    table << "weighted_count skipped";

    table << std::endl;


    for (Candidates::iterator it = m_can.begin(); it != m_can.end(); ++it) {
        // refit to get energy spectrum
        astro::SkyDir dir = it->second.dir();

        table
            << "UW_J"<< int(dir.ra()*10+0.5) 
            << ( dir.dec()>0? "p":"m") << int(fabs(dir.dec())*10+0.5) <<delim// create a name
            << it->second.isSource() << delim
            << dir.ra()  << delim
            << dir.dec() << delim
            << sqrt(it->second.value())  << delim  // sigma
            << 3.*it->second.sigma() << delim // error circle
            << dir.l() << delim
            << dir.b() << delim
            ;
        for( int i=6; i<14; ++i){
            table << precision(it->second.values(i),0.1) << delim;
        }
        // measured number of signal photons per band
        for( int i=6; i<14; ++i){
            table << precision(it->second.photons(i),0.1) << delim;
        }

        // check power law fit
        std::vector<std::pair<double, double> > values;
        values.clear();
        std::vector<double> energyBins = m_pmap.energyBins();
        for( int i=8; i<14; ++i)
        {
            double count = it->second.photons(i);
            if (count > 1e-10)
                values.push_back(std::make_pair(energyBins[i - 6], count));
        }
        pointlike::PowerLawFilter pl(values);
        table << pl.slope() << delim;
        table << pl.constant() << delim;
        table << pl.metric() << delim;

        // Calculate chi_squared
        double chi_sq = 0.0;
        for( int i=8; i<14; ++i)
        {
            double sigma =  1 - it->second.sigalph(i);
            if (fabs(sigma) < 1e-10)
                sigma = 1e-10;
            double chi = (it->second.photons(i)) - (pl.constant() * pow(energyBins[i - 6], pl.slope()))
                / sigma;
            chi_sq += chi * chi;
        }

        table << chi_sq << delim;

        // add weighted count and nbr skipped
        table << it->second.weighted_count() << delim;
        table << it->second.skipped() << delim;
        table << std::endl;

    }
    table.close();
    timer();

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
        (*itor)["EMIN"].set(m_pmap.energyBins().front());
        (*itor)["EMAX"].set(m_pmap.energyBins().back());
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
void SourceFinder::createRegFile(std::string filename, std::string color)const
{
    std::ofstream out(filename.c_str());
    out << "global color="<< color 
        << " font=\"helvetica 10 normal\" select=1 edit=1 move=0 delete=1 include=1 fixed=0 width=2;fk5;"
        << std::fixed << std::setprecision(4) << std::endl;
    for( Candidates::const_iterator it = m_can.begin(); it != m_can.end();  ++it)  {
        const CanInfo& cand = it->second;
        out << "cross point("<< cand.ra()<< ","<<cand.dec() <<") # text={TS=" 
            << int(cand.value()+0.5) << "};\n";
    }
    out.close();
}




void SourceFinder::run()
{
    examineRegion();
    // inital prune
    prune_neighbors();   

    if( refit) {    // group nearby candidates with strongest neighbor   
        group_neighbors();   
        // reexamine the groups of candidates   
        reExamine();   
        // prune the result   
    prune_neighbors();   
    }

    // and write out the ascii table, reg or fits files
    if( ! outfile.empty() )  createTable(outfile);   

    if( ! regfile.empty() ) createRegFile(regfile);

    if( ! fitsfile.empty() ) createFitsFile(fitsfile);
}


