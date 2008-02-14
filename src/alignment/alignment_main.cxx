/** @file alignment_main.cxx
@brief  LAT alignment main program

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/alignment/alignment_main.cxx,v 1.1 2007/08/27 23:24:01 mar0 Exp $

*/
#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/Data.h"
#include "pointlike/AlignProc.h"
#include "pointlike/ParamOptimization.h"
#include "skymaps/PhotonMap.h"
#include "embed_python/Module.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

//#define SELECT

using namespace pointlike;

void help(){
    std::cout << "This program expects a single command-line parameter which is the path to a folder containing a file\n"
        "\talignment_setup.py" 
        << std::endl;

}


class AlignmentSources {
public:
    AlignmentSources(const embed_python::Module& setup )
        :healpixdata(setup)
    {
        // set these in calling program from the Data class
        setup.getList("Data.files",m_filelist);  // get same list of (root for now) files as used in the Data instance

        // get list of candidate sources from the setup module 
        double separation,maxerr;
        std::vector<double> ras,decs;
        std::vector<std::string> names;
        std::vector<astro::SkyDir> sources,temp;
        setup.getList("ra_list", ras);
        setup.getList("dec_list", decs);
        setup.getList("name_list",names);
        setup.getValue("separation",separation,M_PI);
        setup.getValue("maxerr",maxerr,97);
        separation*=M_PI/180;

        std::vector<double>::iterator id = decs.begin();
        int current(0);
        for(std::vector<double>::iterator ir = ras.begin();ir!=ras.end()&&id!=decs.end();++id,++ir) {
            temp.push_back(astro::SkyDir(*ir,*id));
        }
        for(std::vector<astro::SkyDir>::iterator it1=temp.begin();it1!=temp.end();++it1) {
            double min = 1e9;
            for(std::vector<astro::SkyDir>::iterator it2=temp.begin();it2!=temp.end();++it2) {
                double diff = it1->difference(*it2);
                if(min>diff&&diff>0) {
                    min = diff;
                }
            }
            if(min>separation) sources.push_back(*it1);
        }

        int usefitdirection(1);
        setup.getValue("usefitdirection", usefitdirection, 1);
        // run likelihood on all of them, select strong ones

        PointSourceLikelihood::setParameters(setup);
        int minlevel = PointSourceLikelihood::minlevel();
        int maxlevel = PointSourceLikelihood::maxlevel();
        std::cout << std::left << std::setw(20) <<"name" << "     TS   error    ra     dec\n";
        int added = 0;
        for(int i(0);i<sources.size();++i) {
            PointSourceLikelihood ps(healpixdata,names[i],sources[i]);
            ps.maximize();
            double err = ps.localize();
            if(ps.TS()>25 && err<maxerr) { //must be 5-sigma and 1 arcmin resolution 
                ++added;
                used_sources.push_back( usefitdirection!=0? ps.dir() : sources[i] );
                std::vector<double> alphat;
                for(int j(minlevel);j<=maxlevel;++j) {
                    alphat.push_back(ps.find(j)->second->alpha());
                }
                used_alphas.push_back(alphat);
                used_names.push_back(names[i]);
                std::cout << std::left << std::setw(20) << names[i] 
                << std::setprecision(2) << std::setw(8) << std::fixed << std::right
                    << ps.TS() 
                    << std::setprecision(4) 
                    << std::setw(10) << err
                    << std::setw(10) << ps.dir().ra() 
                    << std::setw(10) << ps.dir().dec() 
                    << std::endl;
            }
        }
        std::cout << "Used " << added << " of " << sources.size() << " sources" << std::endl;
    }
    const std::vector<astro::SkyDir>& usedSourceList()const{return used_sources;}
    const std::vector<std::vector<double> >& usedAlphas() const{return used_alphas;}
    const std::vector<std::string>& usedNames() const{return used_names;}
    double tstart()const{return healpixdata.minTime();}
    double tstop()const {return healpixdata.maxTime();}
    const std::vector<std::string>& filelist() { return m_filelist;}
private:
    Data healpixdata;
    std::vector<std::string> m_filelist;
    std::vector<astro::SkyDir> used_sources;
    std::vector<std::vector<double> > used_alphas;
    std::vector<std::string> used_names;

};


int main(int argc, char** argv)
{

    int rc(0);
    try{

        std::string python_path("../python");

        if( argc>1){
            python_path = argv[1];
        }

        embed_python::Module setup(python_path , "alignment_setup",  argc, argv);


        // set up a list of sources for alignment, by fitting them.
        // get the parameters for rerunning the data from the Data object (list of files, time limits)
        AlignmentSources sources(setup);
        double start(sources.tstart()), stop(sources.tstop());

        std::vector<astro::SkyDir> used_sources = sources.usedSourceList();
        std::vector<std::vector<double> > used_alphas = sources.usedAlphas();
        std::vector<std::string> used_names = sources.usedNames();
        std::vector<std::string> filelist = sources.filelist();


        // copy the alignment matrix from the Data class
        AlignProc::set_rot(Data::get_rot());

        // analyze the data
        double resolution;
        setup.getValue("Alignment.resolution",resolution,1);
        // report results
        std::ostream* out = &std::cout;
        std::string outfile;
        setup.getValue("Alignment.outfile",outfile,"");
        if( !outfile.empty() ) {
            out = new std::ofstream(outfile.c_str());
        }
#ifdef SELECT
        for(int j(0);j<used_sources.size();++j) {
            std::vector<astro::SkyDir> used_sources2;
            std::vector<std::vector<double> > used_alphas2;
            used_sources2.push_back(used_sources[j]);
            used_alphas2.push_back(used_alphas[j]);

            AlignProc ap(used_sources2,used_alphas2
#else
        AlignProc ap(used_sources,used_alphas
#endif
            ,filelist, resolution, start, stop);




        std::vector<double> b = ap.alignment();

        //output most likely rotation
        std::cout << "Alignment corrections (arcsec):\n"
            << std::setw(10) << "x" << std::setw(10)<< "y"<< std::setw(10)<<"z" << "\n   ";
#ifdef SELECT
        (*out) << used_names[j]<< std::endl;
#endif
        for(int ir(0);ir<=2;++ir) {
            std::cout << std::setprecision(1) << std::setw(10) << (b[ir]);
            (*out) << std::setprecision(1) << std::setw(10) << (b[ir]);
        }

        std::cout << std::endl << "+/-";
        //output sigma values
        AlignProc::LikeSurface l = ap.fitsurface();
        for(int ir(0);ir<=2;++ir) {
            std::cout << std::setprecision(1) << std::setw(10) << pow(-2*l.curvature()[ir],-0.5);
            (*out) << std::setprecision(1) << std::setw(10) << pow(-2*l.curvature()[ir],-0.5);
        }
        (*out) << std::endl;
        std::cout << std::endl;

#if 0
        std::cout << "\nFit Parameters" << std::endl;
        (*out) << "\nFit Parameters" << std::endl;

        //output fit parameters
        std::cout << std::setprecision(6); // restore
        (*out) << std::setprecision(6);
        for(int k(0);k<10;++k) {
            std::cout << "a[" << k << "]= " <<l(k) << std::endl;
            (*out) << l(k) << std::endl;
        }
#endif
#ifdef SELECT
    }
#endif
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