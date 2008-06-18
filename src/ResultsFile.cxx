#include "pointlike/ResultsFile.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace pointlike{

ResultsFile::ResultsFile (const std::string& filename,const Data& datafile,int nsources):srcTab(0),fileSvc(tip::IFileSvc::instance()) {
      
//       emin = datafile.map().energyBins();
//       emax.resize(emin.size());
//       if (emin.size()<2) throw std::runtime_error("Cannot fit spectrum with only 1 energy bin.");
//       levels=datafile.map().levels();

      std::stringstream efmtstream;	efmtstream<<levels<<"E";  std::string efmt=efmtstream.str();
      std::stringstream jfmtstream;	jfmtstream<<levels<<"J";  std::string jfmt=jfmtstream.str();

//       eratio=emin[1]/emin[0]; 
//       std::copy(++emin.begin(),emin.end(),emax.begin());
//       emax[emin.size()-1]=eratio*emin[emin.size()-1];	

      fileSvc.createFile(filename);
      fileSvc.appendTable(filename,"SOURCES");
      srcTab = fileSvc.editTable(filename,"SOURCES");

      if(srcTab==0) throw std::runtime_error("Could not create SOURCES table.");;

      tip::Header& srcHeader = srcTab->getHeader();
      srcTab->appendField("NAME", "20A");
      srcTab->appendField("TYPE", "10A");
      srcTab->appendField("RA", "1E");
      srcTab->appendField("DEC", "1E");
      srcTab->appendField("L", "1E");
      srcTab->appendField("B", "1E");
      srcTab->appendField("SIGMAX", "1E");
      srcTab->appendField("SIGMAY", "1E");
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
      srcTab->appendField("BG", efmt);
      srcTab->setNumRecords(nsources);
//       srcHeader["EMIN"].set(emin[0]); 
      srcHeader["ERATIO"].set(eratio); 
      srcHeader["LEVELS"].set(levels); 
//       srcHeader["MINLEVEL"].set(datafile.map().minLevel()); 
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
//       (*srcTabItor)["EMIN"].set(emin);
//       (*srcTabItor)["EMAX"].set(emax);
      (*srcTabItor)["SRCPAR"].set(like.sourceParameters());
      (*srcTabItor)["ERR_SRCP"].set(like.sourceParErrors());

//       std::vector<double> alpha(levels);
//       std::vector<double> erra(levels);
//       std::vector<int> nphoton(levels);
//       std::vector<double> nsig(levels);
//       std::vector<double> errns(levels);
//       std::vector<double> ts(levels);
//       std::vector<double> bg(levels);

//       int i=0;
//       for( SourceLikelihood::const_iterator it = like.begin(); it!=like.end(), i<levels; ++it,++i){
//            ExtendedLikelihood& levellike = *it->second;

// 	   std::pair<double,double> a= levellike.maximize();
// 	   alpha[i]=a.first;
// 	   erra[i]=a.second;
// 	   nphoton[i] = levellike.photons();
// 	   nsig[i] = levellike.photons()*a.first/levellike.psfIntegral();
// 	   errns[i] = levellike.photons()*a.second/levellike.psfIntegral();
// 	   ts[i] = levellike.TS();
// 	   bg[i] = levellike.average_b();
//       };

//      (*srcTabItor)["ALPHA"].set(alpha);
//      (*srcTabItor)["ERR_ALPHA"].set(erra);
//      (*srcTabItor)["NPHOTON"].set(nphoton);
//      (*srcTabItor)["NSIGNAL"].set(nsig);
//      (*srcTabItor)["ERR_NSIG"].set(errns);
//      (*srcTabItor)["TS"].set(ts);
//      (*srcTabItor)["BG"].set(bg);
//       srcTabItor++;
 };

 void ResultsFile::writeAndClose(){
      delete srcTab;
      srcTab=0;
   };

};
