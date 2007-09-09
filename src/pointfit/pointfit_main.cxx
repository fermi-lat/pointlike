/** @file pointfit_main.cxx
    @brief  Main program for pointlike localization fits

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfit/pointfit_main.cxx,v 1.11 2007/09/03 23:32:23 burnett Exp $

*/
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/Data.h"
#include "pointlike/SigmaOptimization.h"
#include "pointlike/DiffuseFunction.h"

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


        // separate lists of names, ra's, and dec's
        std::vector<std::string> names;
        std::vector<double> ras, decs;
        std::vector<astro::SkyDir> directions;
        setup.getList("name", names);
        setup.getList("ra", ras);
        setup.getList("dec", decs);

        // flag, if present, to run sigma fitter
        int check_sigma;
        setup.getValue("check_sigma", check_sigma, 0);

        // use the  Data class to create the PhotonData object
        Data healpixdata(setup);

        // define all parameters used by PointSourceLikelihood
        PointSourceLikelihood::setParameters(setup);

        std::ostream* out = &std::cout;
        if( !outfile.empty() ) {
            out = new std::ofstream(outfile.c_str());
        }
        (*out) << std::left << std::setw(15) <<"name" << "     TS   error    ra     dec\n";

        for( size_t n=0; n< names.size(); ++n){
            astro::SkyDir dir(ras[n], decs[n]);
            std::string name(names[n]);

            // fit the point: create the fitting object 
            PointSourceLikelihood like(healpixdata, name, dir);
            // initial fit to all levels at current point
            like.maximize(); 
            // now localize it, return error circle radius
            double sigma =like.localize();

            // add entry to table with name, total TS, localizatino sigma, fit direction
            (*out) << std::left << std::setw(15) << name 
                << std::setprecision(2) << std::setw(8) << std::fixed << std::right
                << like.TS() 
                << std::setprecision(4) 
                << std::setw(10) << sigma
                << std::setw(10) << like.dir().ra() 
                << std::setw(10) << like.dir().dec() 
                << std::endl;

            directions.push_back(like.dir());
        }
        if( check_sigma){
            int minlevel(6), maxlevel(13);
            SigmaOptimization so(healpixdata,directions,out,minlevel,maxlevel);
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

