/** @file pointfit_main.cxx
    @brief  Main program for pointlike localization fits

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfit/pointfit_main.cxx,v 1.24 2008/07/06 06:41:34 burnett Exp $

*/
#include "pointlike/SourceList.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/Data.h"
#include "pointlike/ParamOptimization.h"

#include "embed_python/Module.h"

#include <iostream>
#include <iomanip>
#include <fstream>

void help(){
    std::cout << "This program expects a single command-line parameter which is the path to a folder containing a file\n"
        "\tpointfit_setup.py" 
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

        Module setup(python_path , "pointfit_setup",  argc, argv);
        std::string outfile;
        setup.getValue("outfile",  outfile, "");

        // use the  Data class to create the PhotonData object
        Data healpixdata(setup);

        // print out summary of the the data used?
        healpixdata.info(); 

        // define all parameters used by PointSourceLikelihood
        PointSourceLikelihood::setParameters(setup);

        std::ostream* out = &std::cout;
        if( !outfile.empty() ) {
            out = new std::ofstream(outfile.c_str());
        }
        SourceList::set_data(& healpixdata.map());

        SourceList* sl(0);
        std::map<std::string, std::vector<double> > sourcedict;
        try{
            setup.getDict("sources", sourcedict);
            sl = new SourceList(sourcedict);
        }catch( const std::exception&) {
            std::string sourcelistfile;
            setup.getValue("sourcelistfile", sourcelistfile);
            sl = new SourceList(sourcelistfile);
        }        

        sl->refit(); 
        sl->dump(*out); 

        if( !outfile.empty()){
            delete out;
        }
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

