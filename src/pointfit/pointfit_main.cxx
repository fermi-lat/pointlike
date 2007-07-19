/** @file pointfit_main.cxx
    @brief  Main program for pointlike localization fits

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/pointfit/pointfit_main.cxx,v 1.7 2007/07/19 15:20:59 burnett Exp $

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

double goldensearch(std::vector<astro::SkyDir> directions, pointlike::Data healpixdata, int num2look, int level) {
    double tol = 1e-4;
    double C = (3-sqrt(5.))/2;
    int iter2=0;
    double R = 1-C;
    double ax = 1e-6;
    double bx = 1.;
    double cx = 1.e6;
    double x0 = ax;
    double x3 = cx;
    double x1,x2;
    double xmin,fmin;
    if (fabs(cx-bx) > fabs(bx-ax)) {
        x1 = bx;
        x2 = bx + C*(cx-bx);
    }else {
        x2 = bx;
        x1 = bx - C*(bx-ax);
    }
    double f1 = 0;
    double f2 = 0;
    for(std::vector<astro::SkyDir>::iterator it=directions.begin(); (it!=directions.end())&& (iter2<num2look);++it,++iter2) {
        pointlike::PointSourceLikelihood ps(healpixdata, "test", (*it),7.);
        ps.maximize(2);
        pointlike::PointSourceLikelihood::iterator ite = ps.find(level);
        if(ite->second->photons()==0) continue;
        f1+=ite->second->feval(x1);
        f2+=ite->second->feval(x2);
    }
    if(f1==0||f2==0) return -1;
    int k = 1;
    while(fabs(x3-x0) > tol*(fabs(x1)+fabs(x2))) {
        iter2=0;
        double a1=0;
        double a2=0;
        for(std::vector<astro::SkyDir>::iterator it=directions.begin(); (it!=directions.end())&& (iter2<num2look);++it,++iter2) {
            pointlike::PointSourceLikelihood ps(healpixdata, "test", (*it),7.);
            ps.maximize(2);
            pointlike::PointSourceLikelihood::iterator ite = ps.find(level);
            if(ite->second->photons()==0) continue;
            a1+=ite->second->feval(R*x1+C*x0);
            a2+=ite->second->feval(R*x2+C*x3);
        }
        if (f2 < f1){
            x0 = x1;
            x1 = x2;
            x2 = R*x1 + C*x3;   // x2 = x1+c*(x3-x1)
            f1 = f2;
            f2 = a2;
        }
        else {
            x3 = x2;
            x2 = x1;
            x1 = R*x2 + C*x0;   // x1 = x2+c*(x0-x2)
            f2 = f1;
            f1 = a1;
        }
        k = k+1;
    }
    if (f1 < f2) {
        xmin = x1;
        fmin = f1;
    }else {
        xmin = x2;
        fmin = f2;
    }
    return xmin;
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
        std::vector<astro::SkyDir> directions;
        setup.getList("name", names);
        setup.getList("ra", ras);
        setup.getList("dec", decs);

        // flag, if present, to run sigma fitter
        int check_sigma;
        setup.getValue("check_sigma", check_sigma, 0);

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

            directions.push_back(like.dir());
        }
        if( check_sigma){
            std::cout << "Creating Likelihood set from data" << std::endl;
            for(int i=minlevel;i<=maxlevel;i++) {
                PointSourceLikelihood::set_sigma_level(i,1.);
            }
            int num2look = 10;
            int timeout = 30;
            double tol = 1e-3;
            for(int iter=minlevel;iter<=maxlevel;++iter) {
                int whileit =0;
                double maxfactor = 0;
                double osigma=0;
                while(maxfactor>=0 && fabs(maxfactor-1.)>tol && whileit<timeout){
                    maxfactor = goldensearch(directions,healpixdata,num2look,iter);
                    if(maxfactor>0) {
                        osigma = PointSourceLikelihood::set_sigma_level(iter,pow(maxfactor,-0.5));
                        osigma*= PointSourceLikelihood::set_sigma_level(iter,osigma*pow(maxfactor,-0.5));
                    }
                    whileit++;
                }
                std::cout << "Best value of sigma for level " << iter << " was " << (maxfactor>0?osigma:-1) << " with a scale factor of: " << (maxfactor>0?maxfactor:-1) << std::endl;
            }
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

