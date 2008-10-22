/** @file pointfit_main.cxx
    @brief  Main program for pointlike localization fits

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/src/pointfit/pointfit_main.cxx,v 1.28 2008/08/28 22:42:15 burnett Exp $

*/
#include "pointlike/SourceList.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/Data.h"
#include "pointlike/ParamOptimization.h"

#include "embed_python/Module.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <typeinfo>

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
        if( setup.attribute("sources", false)!=0 ){
            setup.getDict("sources", sourcedict);
            sl = new SourceList(sourcedict);
        }else if( setup.attribute("sourcelistfile", false)!=0 ) {
            std::string sourcelistfile;
            setup.getValue("sourcelistfile", sourcelistfile);
            sl = new SourceList(sourcelistfile);
        }else {
            throw std::invalid_argument(
                "expected either 'sources' or 'sourcelistfile' to define seeds");
        }

        sl->sort_TS(); // initial sort by decreasing TS
        sl->refit(); 
        sl->filter_TS(10); // filter
        sl->sort_ra(); // now by ra
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

