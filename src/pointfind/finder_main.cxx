/** @file finder_main.cxx
    @brief  Finder

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfind/finder_main.cxx,v 1.16 2008/06/23 14:33:44 burnett Exp $

*/
#include "pointlike/SourceFinder.h"
#include "pointlike/Data.h"
#include "pointlike/PointSourceLikelihood.h"

#include "embed_python/Module.h"

#include <iostream>

void help(){
    std::cout <<
        "Usage:  pointfind [path] [setup_file]\n"
        "\n"
        "\tpath        [../python] path to folder containing a python setup file\n"
        "\tsetup_file  [pointfind_setup] name of a python module file (without the .py)\n"
         << std::endl;
        
}


int main(int argc, char** argv)
{
    using namespace astro;
    using namespace pointlike;
    using namespace embed_python;

    int rc(0);
    try{

        std::string python_path("../python"), setup_file("pointfind_setup");

        if( argc>1) python_path = argv[1];
        if( argc>2) setup_file = argv[2];

        Module setup(python_path , setup_file,  argc, argv);
  
        // create healpix database using parameters in the setup file
        Data healpixdata(setup);

        // define all parameters used by PointSourceLikelihood
        PointSourceLikelihood::setParameters(setup);
        SourceFinder::setParameters(setup);

        // create and run the SourceFinder
        pointlike::SourceFinder finder(healpixdata);

        finder.run();

        // try to call the finish function, if defined
        PyObject* finish(setup.attribute("finish", false) );
        if(finish!=0) setup.call(finish);


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
