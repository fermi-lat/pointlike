/** @file pointfit_main.cxx
    @brief  Main program for pointlike localization fits

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfit/pointfit_main.cxx,v 1.3 2007/06/25 20:59:25 burnett Exp $

*/
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/Data.h"
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

        double radius(7.0)
            , TSmin (5)
            ; // 
        int   event_type(0)  // 0: select front only; -1 no selection
            , source_id(-1)  // -1: all sources 
            , minlevel(8)    // minimum level to use for fits (>-6)  
            , maxlevel(13)   // maximum level for fits  (<=13)
            , skip1(2)       // inital number of layers to skip in localization fit
            , skip2(4)       // don't skip beyond this
            , itermax(2)     // maximum number of iterations
            ;
        int verbose(0); // set true to see fit progress

        setup.getValue("radius", radius, 7.0);
        setup.getValue("event_type", event_type, 0);
        setup.getValue("source_id", source_id, -1);
        setup.getValue("TSmin",   TSmin, 5);
        setup.getValue("minlevel", minlevel, 8);
        setup.getValue("maxlevel", maxlevel, 13);
        setup.getValue("skip1", skip1, 2);
        setup.getValue("skip2", skip2, 4);
        setup.getValue("itermax", itermax, 2);
        setup.getValue("verbose", verbose, 0);
    
        std::vector<std::string> filelist;
        std::string ft2file, outfile;
        setup.getList("files", filelist);
        setup.getValue("FT2file",  ft2file, "");
        setup.getValue("outfile",  outfile, "");


        // separate lists of names, ra's, and dec's
        std::vector<std::string> names;
        std::vector<double> ras, decs;
        setup.getList("name", names);
        setup.getList("ra", ras);
        setup.getList("dec", decs);

        // use the  Data class to create the PhotonData object
        Data healpixdata(filelist,  event_type, source_id, ft2file);

        std::ostream* out = &std::cout;
        if( !outfile.empty() ) {
            out = new std::ofstream(outfile.c_str());
        }
        (*out) << std::left << std::setw(15) <<"name" << "     TS   error    ra     dec\n";

        for( size_t n=0; n< names.size(); ++n){
            astro::SkyDir dir(ras[n], decs[n]);
            std::string name(names[n]);

            // fit the point: create the fitting object 

            PointSourceLikelihood like(healpixdata, name, dir, radius, minlevel, maxlevel);
            like.set_verbose(verbose>0);
            
            // initial fit to all levels at current point
            like.maximize(); 

            // now localize it, return error circle radius
            double sigma =like.localize(skip1, skip2, itermax, TSmin);

            // add entry to table with name, total TS, localizatino sigma, fit direction
            (*out) << std::left << std::setw(15) << name 
                << std::setprecision(2) << std::setw(8) << std::fixed << std::right
                << like.TS() 
                << std::setprecision(4) 
                << std::setw(10) << sigma
                << std::setw(10) << like.dir().ra() 
                << std::setw(10) << like.dir().dec() 
                << std::endl;
        }
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

