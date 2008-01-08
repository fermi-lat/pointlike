/** @file finder_main.cxx
    @brief  Finder

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfind/finder_main.cxx,v 1.11 2008/01/02 19:15:01 burnett Exp $

*/
#include "pointlike/SourceFinder.h"
#include "pointlike/Data.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/DiffuseFunction.h"

#include "embed_python/Module.h"

#include <iostream>
#include <iomanip>
#include <fstream>

void help(){
    std::cout << "This program expects a single command-line parameter which is the path to a folder containing a file\n"
        "\tpointfind_setup.py" 
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

        Module setup(python_path , "pointfind_setup",  argc, argv);

        // direction and radius to use: default all sky
        double  radius;
        setup.getValue("radius", radius, 180);
        SkyDir dir(0,0, SkyDir::GALACTIC);
        if( radius<180){
            // if radius is not all sky, look for direction
            double ra, dec;
            setup.getValue("ra", ra, -999);
            setup.getValue("dec", dec,-999);
            if( ra!=-999 && dec !=-999){
                dir=SkyDir(ra, dec);
            }else {
                // ra, dec not specified: use l,b if there
                double l,b;
                setup.getValue("l", l, 0);
                setup.getValue("b", b, 0);
                dir  = SkyDir(l, b, SkyDir::GALACTIC);
            }
        }

        // set up output summary file
        std::string  outfile;
        setup.getValue("outfile",  outfile, "");
        std::ostream* out = &std::cout;
        if( !outfile.empty() ) {
            out = new std::ofstream(outfile.c_str());
        }
  
        // create healpix database using parameters in the setup file
        Data healpixdata(setup);

        // define all parameters used by PointSourceLikelihood
        PointSourceLikelihood::setParameters(setup);

        // create the SourceFinder
        pointlike::SourceFinder finder(healpixdata, setup);

        // look for sources
        finder.examineRegion();

       // group nearby candidates with strongest neighbor
        finder.group_neighbors();

        // reexamine the groups of candidates
        finder.reExamine();

        // prune the result
        finder.prune_neighbors();

        // and write out the table
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

/** @page pointfind The pointfind application



*/
