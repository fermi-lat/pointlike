/** @file finder_main.cxx
    @brief  Finder

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/finder/finder_main.cxx,v 1.1 2007/07/19 14:07:57 burnett Exp $

*/
#include "pointlike/SourceFinder.h"
#include "pointlike/Data.h"
#include "embed_python/Module.h"

#include <iostream>
#include <iomanip>
#include <fstream>

void help(){
    std::cout << "This program expects a single command-line parameter which is the path to a folder containing a file\n"
        "\tfinder_setup.py" 
        << std::endl;
        
}


int main(int argc, char** argv)
{
    using namespace astro;
    using namespace pointlike;
    using namespace embed_python;

    int rc(0);
    try{

        std::string python_path("../python");

        if( argc>1){
            python_path = argv[1];
        }

        Module setup(python_path , "finder_setup",  argc, argv);

        int   pix_level(8)
            , count_threshold(5)
            , skip_TS(2)
            , event_class(-1)    // selection of A/B?
            , source_id
            ;
        double
              ra(250.), dec(-47.)
            , radius(180) // 180 for all sky
            , prune_radius(0.25)//1.0)
            , eq_TS_min(18.)
            , mid_TS_min(18.)
            , polar_TS_min(18.)
            , eq_boundary(6.)
            , TSmin
            ;

        pointlike::SourceFinder::RegionSelector region =
            pointlike::SourceFinder::ALL;

        bool includeChildren (true), 
             weighted( true),
	     background_filter(false);
        int verbose(0); // set true to see fit progress

        setup.getValue("ra", ra);
        setup.getValue("dec", dec);

        setup.getValue("radius",     radius, 180);
        setup.getValue("event_class", event_class, 0);
        setup.getValue("source_id",  source_id, -1);
        setup.getValue("TSmin",      TSmin, 10);
        setup.getValue("prune_radius", prune_radius, 0.25);


        setup.getValue("verbose", verbose, 0);
    
        std::vector<std::string> filelist;
        std::string ft2file, outfile;
 
        setup.getList("files", filelist);
        setup.getValue("outfile",  outfile, "");



        // use the  Data class to create the PhotonData object
        Data healpixdata(filelist,  event_class, source_id);

        std::ostream* out = &std::cout;
        if( !outfile.empty() ) {
            out = new std::ofstream(outfile.c_str());
        }
  

        pointlike::SourceFinder finder(healpixdata);

        astro::SkyDir dir(ra, dec);


        finder.examineRegion( dir, radius, 
            eq_TS_min, mid_TS_min, polar_TS_min, 
            pix_level, 
            count_threshold, true, true, 
            background_filter, 
            skip_TS, 
            region, 
            eq_boundary);
        
        finder.prune_neighbors(prune_radius);

        finder.createTable(outfile);
        
        if( !outfile.empty()){
            delete out;
        }

    }catch(const std::exception& e){
        std::cerr << "Caught exception " << typeid(e).name() 
            << " \"" << e.what() << "\"" << std::endl;
        help();
        rc=1;
    }
     return rc;
}

