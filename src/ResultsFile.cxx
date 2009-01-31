#include "pointlike/ResultsFile.h"
#include "skymaps/Band.h"
#include "pointlike/SpectralFitter.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace pointlike{

  ResultsFile::ResultsFile (const std::string& filename,const Data& datafile,int nsources)
     :srcTab(0)
     ,fileSvc(tip::IFileSvc::instance())
     ,fields_created(false) {
      
    std::list<skymaps::Band>::const_iterator SpecIt = datafile.map().begin();
    for(;SpecIt!= datafile.map().end(); SpecIt++){      
      emin.push_back((*SpecIt).emin());
      emax.push_back((*SpecIt).emax());
    }; 
    if (emin.size()<2) throw std::runtime_error("Cannot fit spectrum with only 1 energy bin.");
    levels=datafile.map().size();
    eratio=emax[0]/emin[0];      
    
    std::stringstream efmtstream;	efmtstream<<levels<<"E";  std::string efmt=efmtstream.str();
    std::stringstream jfmtstream;	jfmtstream<<levels<<"J";  std::string jfmt=jfmtstream.str();
      
    fileSvc.createFile(filename);
    fileSvc.appendTable(filename,"SOURCES");
    srcTab = fileSvc.editTable(filename,"SOURCES");
    
    if(srcTab==0) throw std::runtime_error("Could not create SOURCES table.");;
    
    tip::Header& srcHeader = srcTab->getHeader();
    srcTab->appendField("NAME", "30A");
    srcTab->appendField("TYPE", "10A");
    srcTab->appendField("RA", "1E");
    srcTab->appendField("DEC", "1E");
    srcTab->appendField("L", "1E");
    srcTab->appendField("B", "1E");
    srcTab->appendField("SIGMAX", "1E");
    srcTab->appendField("SIGMAY", "1E");
    srcTab->appendField("R0", "1E");
    srcTab->appendField("ERR_R0", "1E");
    srcTab->appendField("SUMTS", "1E");
    srcTab->appendField("SRCPAR", "5E");
    srcTab->appendField("ERR_SRCP", "5E");
    srcTab->appendField("EMIN", efmt);
    srcTab->appendField("EMAX", efmt);
    srcTab->appendField("ALPHA", efmt);
    srcTab->appendField("ERR_ALPHA", efmt);
    srcTab->appendField("NPHOTON", jfmt);
    srcTab->appendField("NSIGNAL", efmt);
    srcTab->appendField("ERR_NSIG", efmt);
    srcTab->appendField("TS", efmt);
    srcTab->appendField("EVENT_CLASS", jfmt);
    srcTab->appendField("BG", efmt);
    srcTab->appendField("MINOS_PARABOLIC","7E");
    srcTab->appendField("MINOS_PLUS","7E");
    srcTab->appendField("MINOS_MINUS","7E");
    srcTab->appendField("GAMMA", efmt);
    srcTab->appendField("SIGMA", efmt);
    srcTab->setNumRecords(nsources);
    srcHeader["EMIN"].set(emin[0]); 
    srcHeader["ERATIO"].set(eratio); 
    srcHeader["LEVELS"].set(levels); 
    srcHeader["NSOURCES"].set(nsources); 


    srcTabItor=srcTab->begin();
    
  }; 
  
  void ResultsFile::fill(SourceLikelihood& like){

    if(srcTab==0) throw std::runtime_error("Error in Fill: file already closed.");
      
    std::cout<<"setting values"<<std::endl;
    (*srcTabItor)["NAME"].set(like.name());
    (*srcTabItor)["TYPE"].set(like.type());
    (*srcTabItor)["RA"].set(like.dir().ra());
    (*srcTabItor)["DEC"].set(like.dir().dec()); 
    (*srcTabItor)["L"].set(like.dir().l());
    (*srcTabItor)["B"].set(like.dir().b());
    (*srcTabItor)["SIGMAX"].set(like.errorX());
    (*srcTabItor)["SIGMAY"].set(like.errorY());
    (*srcTabItor)["EMIN"].set(emin);
    (*srcTabItor)["EMAX"].set(emax);
    (*srcTabItor)["SRCPAR"].set(like.sourceParameters());
    (*srcTabItor)["ERR_SRCP"].set(like.sourceParErrors());
    (*srcTabItor)["R0"].set(like.sourceParameters()[0]);
    (*srcTabItor)["ERR_R0"].set(like.sourceParErrors()[0]);
 
    std::map< std::string, std::vector<double> > minosErrors=like.errorsMINOS();
    (*srcTabItor)["MINOS_PARABOLIC"].set(minosErrors["parabolic"]);
    (*srcTabItor)["MINOS_PLUS"].set(minosErrors["plus"]);
    (*srcTabItor)["MINOS_MINUS"].set(minosErrors["minus"]);

    std::vector<double> alpha(levels,0);
    std::vector<double> erra(levels,0);
    std::vector<double> flux(levels,0);
    std::vector<double> errflux(levels,0);
    std::vector<double> exposure(levels,0);
    std::vector<int> nphoton(levels,0);
    std::vector<double> nsig(levels,0);
    std::vector<double> errns(levels,0);
    std::vector<double> ts(levels,0);
    std::vector<double> evclass(levels,0);
    std::vector<double> bg(levels,0);
    std::vector<double> sigma(levels,0);
    std::vector<double> gamma(levels,0);

    int i=0;
    double sumTS=0.;
    for( SourceLikelihood::const_iterator it = like.begin(); it!=like.end() && i<levels; ++it,++i){
      ExtendedLikelihood& levellike = **it;
      while ((emin[i]+1e-5)<levellike.band().emin()) i++;
      
      std::pair<double,double> a= levellike.maximize();
      // Fix problem with exposure energy range
      //std::pair<double,double> f= levellike.flux();
      alpha[i]=a.first;
      erra[i]=a.second;
      nphoton[i] = levellike.photons();
      nsig[i] = levellike.photons()*a.first/levellike.psfIntegral();
      errns[i] = levellike.photons()*a.second/levellike.psfIntegral();
      ts[i] = levellike.TS();
      evclass[i] = levellike.band().event_class();
      bg[i] = levellike.average_b();
      gamma[i] = levellike.gamma();
      sigma[i] = levellike.sigma();
      sumTS+=ts[i];
    };

    (*srcTabItor)["SUMTS"].set(sumTS);
    (*srcTabItor)["ALPHA"].set(alpha);
    (*srcTabItor)["ERR_ALPHA"].set(erra);
    (*srcTabItor)["NPHOTON"].set(nphoton);
    (*srcTabItor)["NSIGNAL"].set(nsig);
    (*srcTabItor)["ERR_NSIG"].set(errns);
    (*srcTabItor)["TS"].set(ts);
    (*srcTabItor)["EVENT_CLASS"].set(evclass);
    (*srcTabItor)["BG"].set(bg);
    (*srcTabItor)["GAMMA"].set(gamma);
    (*srcTabItor)["SIGMA"].set(sigma);

    srcTabItor++;
  };

  void ResultsFile::fill(SourceLikelihood& like,
			 SpectralFitter& fitter){

    this->fill(like);
    srcTabItor--;
    
    if(!fields_created) srcTab->appendField("MODEL", "20A");
    (*srcTabItor)["MODEL"].set(fitter.getModel()->get_spec_type());
    
    std::vector<double> range(2);
    range[0]=fitter.getModel()->get_lower_bound();
    range[1]=fitter.getModel()->get_upper_bound();
    if(!fields_created) srcTab->appendField("FIT_RANGE","2E");
    (*srcTabItor)["FIT_RANGE"].set(range);

    std::stringstream nparstream;
    nparstream<<fitter.getModel()->get_npar()<<"E";
    std::string npar=nparstream.str();
    
    std::stringstream nbinstream;
    nbinstream<<fitter.getModel()->get_model_energies().size()<<"E";
    std::string nbin=nbinstream.str();
    
    if(!fields_created) srcTab->appendField("LOGLIKESUM","1E");
    (*srcTabItor)["LOGLIKESUM"].set(fitter.getModel()->get_logLikeSum());
  
    if(!fields_created) srcTab->appendField("PAR",npar);
    (*srcTabItor)["PAR"].set(fitter.getModel()->get_params());
  
    if(!fields_created) srcTab->appendField("PAR_ERR",npar);
    (*srcTabItor)["PAR_ERR"].set(fitter.getModel()->get_param_errors());

    if(!fitter.get_covar_entries().empty()){
      std::stringstream ncovarstream;
      ncovarstream<<fitter.get_covar_entries().size()<<"E";
      std::string ncovar=ncovarstream.str();
      if(!fields_created) srcTab->appendField("COVAR_ENTRIES",ncovar);
      (*srcTabItor)["COVAR_ENTRIES"].set(fitter.get_covar_entries());
    }

    if(fitter.getModel()->get_spec_type()=="POWER_LAW"){
      if(!fields_created) srcTab->appendField("DECORRELATION_ENERGY","1E");
      (*srcTabItor)["DECORRELATION_ENERGY"].set(fitter.getModel()->get_scale());
    }
    
    // Record model integrated flux 100 MeV < E < 300 GeV [ ph/cm^2/s ]
    if(!fields_created) srcTab->appendField("INTEGRAL_FLUX","1E");
    (*srcTabItor)["INTEGRAL_FLUX"].set(fitter.getModel()->get_flux(100,3.e5));
    
    // Record flux errors (high/low)
    if(!fields_created) srcTab->appendField("INTEGRAL_FLUX_ERR","2E");
    (*srcTabItor)["INTEGRAL_FLUX_ERR"].set(fitter.getModel()->get_flux_errors(100,3.e5));
    
    if(!fields_created) srcTab->appendField("MODEL_ENERGY",nbin);
    (*srcTabItor)["MODEL_ENERGY"].set(fitter.getModel()->get_model_energies());

    if(!fields_created) srcTab->appendField("MODEL_EXPOSURE",nbin);
    (*srcTabItor)["MODEL_EXPOSURE"].set(fitter.getModel()->get_exposures());

    if(!fitter.getFluxUpperLimits().empty()){
      if(!fields_created) srcTab->appendField("UPPER_LIMIT_RANGE","2E");
      (*srcTabItor)["UPPER_LIMIT_RANGE"].set(fitter.getUpperLimitRange());

      std::stringstream nclstream;
      nclstream<<fitter.getConfidenceLimits().size()<<"E";
      std::string ncl=nclstream.str();
      if(!fields_created) srcTab->appendField("CONFIDENCE_LIMIT",ncl);
      if(!fields_created) srcTab->appendField("FLUX_UPPER_LIMIT",ncl);
      (*srcTabItor)["CONFIDENCE_LIMIT"].set(fitter.getConfidenceLimits());
      (*srcTabItor)["FLUX_UPPER_LIMIT"].set(fitter.getFluxUpperLimits());

      std::stringstream nbandstream;
      nbandstream<<fitter.getBandUpperLimits().size()<<"E";
      std::string nband=nbandstream.str();
      if(!fields_created) srcTab->appendField("BAND_UPPER_LIMIT",nband);
      if(!fields_created) srcTab->appendField("ENERGY_UPPER_LIMIT",nband);
      if(!fields_created) srcTab->appendField("EXPOSURE_UPPER_LIMIT",nband);
      (*srcTabItor)["BAND_UPPER_LIMIT"].set(fitter.getBandUpperLimits());
      (*srcTabItor)["ENERGY_UPPER_LIMIT"].set(fitter.getEnergyUpperLimits());
      (*srcTabItor)["EXPOSURE_UPPER_LIMIT"].set(fitter.getExposureUpperLimits());
    }
    fields_created=true;
    srcTabItor++;

  };

  void ResultsFile::writeAndClose(){
    delete srcTab;
    srcTab=0;
  };

};
