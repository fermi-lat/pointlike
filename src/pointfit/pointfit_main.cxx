/** @file pointfit_main.cxx
    @brief  Main program for pointlike localization fits

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfit/pointfit_main.cxx,v 1.32 2009/03/08 01:36:23 burnett Exp $

*/
#include "pointlike/SourceList.h"
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/Data.h"
#include "pointlike/ParamOptimization.h"

#include "embed_python/Module.h"
#include "skymaps/IParams.h"
#include "skymaps/BinnedPhotonData.h"

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

void fix_psf(pointlike::Data &data)
{
    using std::setw;
    using std::setprecision;
    using std::left;
    using std::fixed;

    int w(12);

    skymaps::BinnedPhotonData& pdata  = data.map();
    std::cout << " Setting new values for PSF parameters" << std::endl;
    std::cout << setw(w)<<"emin" << setw(w) << "class" << setw(w) << "sigma" << setw(w) << "gamma" << std::endl;
    for( skymaps::BinnedPhotonData::iterator it = pdata.begin(); it!=pdata.end(); ++it){
        skymaps::Band& b = *it;
        double ebar( sqrt(b.emin()*b.emax()) );
        int eclass(b.event_class());
        double 
            sigma( skymaps::IParams::sigma(ebar, eclass) ),
            gamma( skymaps::IParams::gamma(ebar, eclass) );
        b.setSigma(sigma);
        b.setGamma(gamma);

        std::cout << setw(w) << "  "<<setprecision(0)<< fixed<< left<< b.emin() 
            << setw(w) << eclass
            << setprecision(5) 
            << setw(w) << b.sigma() 
            << setw(w) << setprecision(2) <<b.gamma() 
            << std::endl; 
               
    }
    std::cout << setprecision(5);

}

int main(int argc, char** argv)
{
    using namespace astro;
    using namespace pointlike;
    using namespace embed_python;

    int rc(0);
    try{

        skymaps::IParams::init();

        std::string python_path("../python");

        if( argc>1){
            python_path = argv[1];
        }

        Module setup(python_path , "pointfit_setup",  argc, argv);
        std::string outfile;
        setup.getValue("outfile",  outfile, "");

        // use the  Data class to create the PhotonData object
        Data healpixdata(setup);

        fix_psf(healpixdata);

        // print out summary of the the data used?
        healpixdata.info(); 

        // define all parameters used by PointSourceLikelihood
        PointSourceLikelihood::setParameters(setup);

        std::ostream* out = &std::cout;
        if( !outfile.empty() ) {
            out = new std::ofstream(outfile.c_str());
        }
        SourceList::set_data(& healpixdata.map());
        std::string logfile;
        setup.getValue("logfile", logfile, logfile);
        std::ofstream* log(0);
        if( !logfile.empty() ){
            log = new std::ofstream(logfile.c_str());
            SourceList::set_log(log);
            std::cout << "logging to file " << logfile<< std::endl;
        }

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

        double tsmin(10);
        setup.getValue("tsmin", tsmin, tsmin);

        sl->sort_TS(); // order for refit
        sl->refit(); 
        if(tsmin>0) sl->filter_TS(tsmin); // filter
        sl->sort_ra(); // now by ra
        if( !outfile.empty() && outfile.find(".xml") != std::string::npos){
            sl->dump_xml(*out);
        }else{
            sl->dump(*out); 
        }

        if( !outfile.empty()){
            delete out;
        }
        // try to call the finish function, if defined
        PyObject* finish(setup.attribute("finish", false) );
        if(finish!=0) setup.call(finish);

        delete log; // close
        
    }catch(const std::exception& e){
        std::cerr << "Caught exception " << typeid(e).name() 
            << " \"" << e.what() << "\"" << std::endl;
        help();
        rc=1;
    }
     return rc;
}

