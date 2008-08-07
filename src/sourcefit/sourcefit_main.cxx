/** @file pointfit_main.cxx
    @brief  Main program for pointlike localization fits

    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/sourcefit/sourcefit_main.cxx,v 1.3 2008/06/27 05:59:26 markusa Exp $

*/
#include "pointlike/SourceLikelihood.h"
#include "pointlike/Data.h"
#include "pointlike/ParamOptimization.h"
#include "pointlike/ResultsFile.h"

#include "embed_python/Module.h"
#include "tip/IFileSvc.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>

void help(){
  std::cout << "This program expects a single command-line parameter "
	    << " which is the path to a folder containing a file\n"
	    << "\tsourcefit_setup.py" 
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
    std::string setup_file("sourcefit_setup");
    tip::IFileSvc& fileSvc=tip::IFileSvc::instance();
    
    if( argc>1) python_path = argv[1];
    if( argc>2) setup_file = argv[2];
    
    Module setup(python_path , setup_file,  argc, argv);
    std::string outfile;
    setup.getValue("outfile",  outfile, "");
    
    // separate lists of names, ra's, and dec's
    std::vector<std::string> names;
    std::vector<double> ras, decs;
    std::vector<astro::SkyDir> directions;
    std::vector< SourceLikelihood*> likelihoods;
    std::vector< SourceLikelihood*>::iterator like_i;
    std::vector< SourceLikelihood*> bg_likelihoods;
    std::vector< SourceLikelihood*>::iterator like_j;
    std::vector<std::string> types;
    std::vector<double> npars;
    std::vector<double> dofits;
    std::vector<std::string> initvals;
    double bgROI=0.5;
    
    setup.getList("name", names);
    setup.getList("ra", ras);
    setup.getList("dec", decs);
    setup.getList("srctype", types);
    setup.getList("npar", npars);
    setup.getList("init", initvals);
    setup.getList("fit", dofits);
    setup.getValue("bgROI", bgROI);
    
    // flag, to designate first candidate as a central value
    int first_is_center(0);
    setup.getValue("first_is_center", first_is_center, 0);
    
    // flag, if present, to run sigma/gamma fitter
    int check_sigma(0);
    setup.getValue("check_sigma", check_sigma, check_sigma);
    setup.getValue("fitPSF", check_sigma, check_sigma);
    
    // use the  Data class to create the PhotonData object
    Data healpixdata(setup);
    
    // print out summary of the the data used?
    healpixdata.info(); 
    
    // define all parameters used by PointSourceLikelihood
    SourceLikelihood::setParameters(setup);
    
    std::ostream* out = &std::cout; 
     	ResultsFile* results = 0;
    	if( !outfile.empty() ) {
     	    results=new ResultsFile(outfile,healpixdata,names.size());
     	}
    
    (*out) << std::left << std::setw(20) <<"name" << "     TS   error    ra     dec\n";
    size_t n=0;
    for( ; n< names.size(); ++n){
      astro::SkyDir dir(ras[n], decs[n]);
      std::string name(names[n]);
      std::string type(types[n]);
      
      unsigned int npar((unsigned int)(npars[n]+0.5));
      double dofit=dofits[n];
      std::string inistr=initvals[n];
      std::vector<double> srcpar(npar,0);
      std::stringstream inistream(inistr);
      for(unsigned int i=0;i<npar;i++) inistream >> srcpar[i];
      
      //	    for(int i=0;i<npar;i++) std::cout<<srcpar[i]<<std::endl;
      
      // fit the point: create the fitting object
      SourceLikelihood* plike =new SourceLikelihood(healpixdata, name, dir, type, srcpar);
      bg_likelihoods.push_back(plike);
      if (dofit>0) likelihoods.push_back(plike);
    }
    
    for (like_i = likelihoods.begin(); like_i != likelihoods.end() ; like_i++){
      //delete all background sources from previous fit
      (*like_i)->clearBackgroundPointSource();
      //make all sources except the current one part of the background
      for (like_j = bg_likelihoods.begin(); like_j != bg_likelihoods.end() ; like_j++)
	if((*like_i)->name()!=(*like_j)->name()) {
	  (*like_j)->maximize();
	  double angle=std::acos( (*like_i)->dir()().dot((*like_j)->dir()()) );
	  //		    std::cout<< (*like_j)->name()<<" "<<(angle*180./M_PI)<<std::endl;
	  if(angle<bgROI) (*like_i)->addBackgroundPointSource(*like_j);
	}    
      //run fit
      SourceLikelihood& like =*(*like_i) ; 
      // initial fit to all levels at current point
      like.maximize(); 
      // now localize it, return error circle radius

      double sigma =like.localize();
      
      // add entry to table with name, total TS, localizatino sigma, fit direction
      (*out) << std::left << std::setw(20) << like.name() 
	     << std::setprecision(2) << std::setw(8) << std::fixed << std::right
	     << like.TS() 
	     << std::setprecision(4) 
	     << std::setw(10) << sigma
	     << std::setw(10) << like.dir().ra() 
	     << std::setw(10) << like.dir().dec() 
	     << std::endl;
      
      directions.push_back(like.dir());
      if( n>0 && first_is_center!=0){
	like.addBackgroundPointSource(likelihoods[0]);
      }
     
 	    if(results) results->fill(like);
      
      directions.push_back(like.dir());
    }

    if( check_sigma){
#ifdef OLD
      int minlevel(6), maxlevel(13);
      ParamOptimization so(healpixdata,directions,out,minlevel,maxlevel);
      
#if 0
      so.compute(ParamOptimization::SIGMA);
      so.compute(ParamOptimization::GAMMA);
#else
      so.compute();
#endif
      
#endif

    }
	
   if(results) results->writeAndClose();


  } catch(const std::exception& e){
    std::cerr << "Caught exception " << typeid(e).name() 
	      << " \"" << e.what() << "\"" << std::endl;
    help();
    rc=1;
  }
  return rc;
}
