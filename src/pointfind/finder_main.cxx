/** @file finder_main.cxx
    @brief  Finder

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfind/finder_main.cxx,v 1.12 2008/01/08 22:40:15 burnett Exp $

*/
#include "pointlike/SourceFinder.h"
#include "pointlike/Data.h"
#include "pointlike/PointSourceLikelihood.h"

#include "embed_python/Module.h"

#include <iostream>

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
  
        // create healpix database using parameters in the setup file
        Data healpixdata(setup);

        // define all parameters used by PointSourceLikelihood
        PointSourceLikelihood::setParameters(setup);

        // create and run the SourceFinder
        pointlike::SourceFinder finder(healpixdata, setup);

        finder.run();



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
