/** @file SourceFinder.cxx
@brief implementation of SourceFinder

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/SourceFinder.cxx,v 1.11 2007/09/12 02:44:35 burnett Exp $
*/

#include "pointlike/SourceFinder.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/PowerLawFilter.h"

#include <fstream>
#include <iostream>
#include <time.h>
#include <cmath>

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

} // anon namespace

using namespace map_tools;
using namespace astro;
using namespace pointlike;
using namespace embed_python;

SourceFinder::SourceFinder(const pointlike::Data& map, Module & Mod)
: m_pmap(map)
, m_counts(0)
, m_module(Mod)
{}

SourceFinder::SourceFinder(const pointlike::CalData& map,Module & Mod)
: m_pmap(map)
, m_counts(0)
, m_module(Mod)
{}


/** @brief
Analyze range of likelihood significance values for all pixels at a particular level  
*/
void SourceFinder::examineRegion(void) 
{  
    int skip1(2), skip2(3);
    double min_alpha(0.15); // miniumum alpha for TS computation
    double sigma_max(0.25); // maximum allowable sigma

    // Get parameters
    double l, b;
    m_module.getValue("l", l, 0);
    m_module.getValue("b", b, 0);
    astro::SkyDir dir(l,b, astro::SkyDir::GALACTIC);

    double  radius;
    m_module.getValue("radius", radius, 180);

    double  eq_TS_min;
    m_module.getValue("eqTSmin", eq_TS_min, 10);

    double  mid_TS_min;
    m_module.getValue("midTSmin", mid_TS_min, 10);

    double  polar_TS_min;
    m_module.getValue("polarTSmin", polar_TS_min, 10);

    int  pix_level;
    m_module.getValue("pixLevel", pix_level, 8);

    int  count_threshold;
    m_module.getValue("countThreshold", count_threshold, 16);

    int  final_pix_lvl;
    m_module.getValue("finalPixlvl", final_pix_lvl, 8);

    int  final_count_threshold;
    m_module.getValue("finalCntThresold", final_count_threshold, 2);

    int  skip_TS_levels;
    m_module.getValue("skipTSlevels", skip_TS_levels, 2);

    double  equator_boundary;
    m_module.getValue("eqBoundary", equator_boundary, 6);

    double  polar_boundary;
    m_module.getValue("polarBoundary", polar_boundary, 40);

    timer("---------------SourceFinder::examineRegion----------------");  
    std::vector<std::pair<astro::HealPixel, int> > v;
    std::cout << "\nAnalyzing likelihood significance using:\n"
        << "  " << radius << " degree radius about gal(" << dir.l() << ", " << dir.b() << 
        ") =eq("<<dir.ra()<<", "<< dir.dec()<< ")\n"
        << "  Pixelization level: " << pix_level << "\n"
        << "  Weighted count threshold: " << count_threshold << "\n"
        << "  Likelihood minimum: (boundaries: " <<equator_boundary <<", "<< polar_boundary<<")\n"
        << "    Equatorial " << eq_TS_min << "\n"
        << "    Middle     " << mid_TS_min << "\n"
        << "    Polar      " << polar_TS_min << "\n"
        << "  skip when localizing: " << skip1<< " to " << skip2  <<"\n"
        << "  min alpha: " << min_alpha << "\n"
        << "  sigma_max: " << sigma_max << "\n"
        << "  skip TS levels: " << skip_TS_levels << "\n"
        << "  Equatorial boundary: " << equator_boundary << "\n"
        << "  Polar boundary: " << polar_boundary << "\n"
        << std::endl;

    m_pmap.extract_level(dir, radius, v, pix_level, true);
    std::cout << v.size() << " pixels will be examined.\n";
    Candidates can;
    can.clear();
    int num(0);
    for(std::vector<std::pair<astro::HealPixel, int> >::const_iterator it = v.begin();
        it != v.end(); ++it)
    {
        ShowPercent(num,v.size(),can.size());
        astro::SkyDir sd;
        //if (it->first.level() != pix_level) continue; // Only want to examine at pixelization level at this point.
        double abs_b = fabs((it->first)().b());
        int count = static_cast<int>(m_pmap.photonCount(it->first, sd));
        
        if (count >= count_threshold)
            can[it->first] = CanInfo(count, 0, sd);
        ++num;
    }

    std::cout << can.size() << " pixels at level " << pix_level << " passed weighted count test.\n";
    //        prune_neighbors(can);

    m_can.clear();
    int i(0);

    for (Candidates::const_iterator it = can.begin(); it != can.end(); ++ it)
    {
        ShowPercent(i++, can.size(), m_can.size());

        const astro::SkyDir& currentpos = it->second.dir();
        double ts_min, abs_b = fabs(currentpos.b());
        if (abs_b < equator_boundary)
            ts_min = eq_TS_min;
        else if (abs_b > polar_boundary)
            ts_min = polar_TS_min;
        else
            ts_min = mid_TS_min;


        // perform likelihood analysis at the current candidate position  
        PointSourceLikelihood ps(m_pmap, "test",currentpos );

        double ts = ps.maximize(skip_TS_levels);
        if (ts < ts_min) continue;  // apply initial threshold

        // adjust position to maximize likelhood
        double error = ps.localize(skip1, skip2);
        if (error >= sigma_max) continue; // quit if fail to find maximum or > 1 degree
        ts = ps.maximize(skip_TS_levels); // readjust likelihood at current position
        if (ts <ts_min) continue;

        // found candidate

        // also check number of photons in pixel
        HealPixel px_check(ps.dir(), final_pix_lvl);
        int count = static_cast<int>(m_pmap.photonCount(px_check, true, false));
        if (count >= final_count_threshold)
        {  

            // add to the final list, indexed according to level 13 location
            HealPixel px(ps.dir(), 13); 
            m_can[px] = CanInfo(ts, error, ps.dir());
            for(int id = 6;id<14;++id)
            {
                m_can[px].setValue(id,ps.levelTS(id));
                m_can[px].setPhotons(id,ps[id]->photons()*ps[id]->alpha());
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
        }
    }
    std::cout << m_can.size() << " sources found before pruning neighbors.\n";
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

    astro::HealPixel hp(sd, pix_level);
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
    m_module.getValue("plEqTSmin", plEqTSmin, 38.0);

    double  plMidTSmin;
    m_module.getValue("plMidTSmin", plMidTSmin, 19.0);

    double  plPolarTSmin;
    m_module.getValue("plPolarTSmin", plPolarTSmin, 24.0);

    double  plSlopeCutoff;
    m_module.getValue("plSlopeCutoff", plSlopeCutoff, -1.5);

    double  plFitCutoff;
    m_module.getValue("plFitCutoff", plFitCutoff, 0.9);

    double  eqBoundary;
    m_module.getValue("eqBoundary", eqBoundary, 6.0);

    double  polarBoundary;
    m_module.getValue("polarBoundary", polarBoundary, 6.0);

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
// criterion is closer that tolerance, or 3-sigma circles overlap
void SourceFinder::prune_neighbors(void)
{
    double  radius;
    m_module.getValue("prune_radius", radius, 0.25);
    
    std::cout << "Eliminating weaker neighbors using radius of " << radius << " degrees...";

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
            if (diff <= radius || diff < 3*(sigma1+sigma2) ) {
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
        std::vector<astro::HealPixel> hv = it->first.neighbors();
        for (std::vector<astro::HealPixel>::const_iterator n = hv.begin();
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
        reg << "fk5; circle("
            << (it->first)().ra() << ", "
            << (it->first)().dec() << ", ";
        if (radius > 0.0)
            reg << radius;
        else
            reg << it->second.sigma() * 100.0;
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
    table << "pl_m pl_b pl_sigma pl_chi_sq";

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

        // table << pl.chi_sq() << delim;
        table << chi_sq << delim;
        table << std::endl;

    }
    table.close();
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



