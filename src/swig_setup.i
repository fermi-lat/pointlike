%module(docstring="Interface to pointlike") pointlike
// $Header: /nfs/slac/g/glast/ground/cvs/users/echarles/healpix_changes/pointlike/src/swig_setup.i,v 1.4 2015/03/05 19:58:43 echarles Exp $
#define PySwigIterator pointlike_PySwigIterator
%{
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <utility>

// EAC, added ProjBase base class
#include "astro/ProjBase.h"
#include "astro/SkyProj.h"
#include "astro/Photon.h"
#include "astro/PointingHistory.h"
#include "astro/PointingInfo.h"
#include "astro/EarthCoordinate.h"
#include "astro/GPS.h"

#include "healpix/Healpix.h"
#include "healpix/HealPixel.h"

#include "embed_python/Module.h"

#include "pointlike/PointSourceLikelihood.h"
#include "pointlike/SimpleLikelihood.h"
#include "pointlike/SourceFinder.h"
#include "pointlike/Data.h"
#include "pointlike/Draw.h"
#include "pointlike/SourceList.h"
#include "pointlike/Alignment.h"
#include "pointlike/LeastSquaresFitter.h"
#include "pointlike/ConfidenceLevel.h"

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/EulerAngles.h"

#include "CLHEP/Vector/ThreeVector.h"

%}
%include stl.i
%exception {
   try {
      $action
   } catch (const std::exception & eObj) {
      if( strcmp(eObj.what(),"StopIteration")==0 ){
          PyErr_SetString(PyExc_StopIteration, const_cast<char *>(eObj.what()));
      } else if(strcmp(eObj.what(),"IndexError")==0 ){
          PyErr_SetString(PyExc_IndexError, const_cast<char *>(eObj.what()));
      } else {
          PyErr_SetString(PyExc_RuntimeError, const_cast<char*>(eObj.what()));
      }
      return NULL;
   }
}

%template(DoublePair) std::pair<double, double>;
%template(StringVector) std::vector<std::string>;
%template(DoubleVector) std::vector<double>;
%template(FloatVector) std::vector<float>;
%template(LongVector) std::vector<long>;
%template(IntVector) std::vector<int>;

%include CLHEP/Vector/ThreeVector.h
//%include $(CLHEPROOT)/include/CLHEP/Vector/Rotation.h
namespace CLHEP {
 class HepRotation {
public:
};
}
%extend CLHEP::Hep3Vector{
// for convenience: make it behave like array of 3 elements
   double __getitem__(size_t i) {
      switch (i){
      case 0: return self->x();
      case 1: return self->y();
      case 2: return self->z();
      case 3: throw std::range_error("StopIteration"); //must be exactly this string
      default: throw std::range_error("IndexError");
      }
   }
   size_t __len__() {      return 3;       }
}

%extend CLHEP::HepRotation{
   double __getitem__(size_t i){
   switch(i){
      case 0: return self->xx(); case 1: return self->xy(); case 2: return self->xz();
      case 3: return self->yx(); case 4: return self->yy(); case 5: return self->yz();
      case 6: return self->zx(); case 7: return self->zy(); case 8: return self->zz();
      case 9: throw std::range_error("StopIteration"); //must be exactly this string
      default: throw std::range_error("IndexError");
      }
   }
   size_t __len__() {return 9;}
   // need this, cannot parse the many includes
   CLHEP::HepRotation inverse(){
      return self->inverse();
      }
}

//%template(IntPair) std::pair<int,int>;
//%template(IntPairVector) std::vector<std::pair<int,int> >;
//%template(WeightedSkyDirVector) std::vector<skymaps::WeightedSkyDir>;

%extend pointlike::PointSourceLikelihood{
   pointlike::SimpleLikelihood * __getitem__(size_t i){ 
      if( i == (*self).size() ) throw std::range_error("StopIteration");
      if( i<0 || i > self->size() ) throw std::range_error("IndexError");
      return (*self)[i]; 
   }
   size_t __len__() {      return self->size();       }
}

// EAC, added ProjBase base class
%include astro/ProjBase.h
%include astro/SkyProj.h
%include astro/Photon.h
%include astro/PointingHistory.h
%include astro/PointingInfo.h
%include astro/EarthCoordinate.h
%include astro/GPS.h
%include astro/Quaternion.h
%include astro/JulianDate.h


%include healpix/Healpix.h

// fails now
//%include $(HEALPIXROOT)/healpix/HealPixel.h

%include embed_python/Module.h

%include pointlike/SimpleLikelihood.h

%include pointlike/PointSourceLikelihood.h
%include pointlike/Data.h
%include pointlike/SourceFinder.h
%include pointlike/Draw.h
//%include $(POINTLIKEROOT)/pointlike/SimpleTSmap.h
%include pointlike/SourceList.h

%include pointlike/LeastSquaresFitter.h


%extend pointlike::SourceList{
   pointlike::Source & __getitem__(size_t i){ 
      if( i == (*self).size() ) throw std::range_error("StopIteration");
      if( i<0 || i > self->size() ) throw std::range_error("IndexError");
      //return (*self)[i]; 
      pointlike::SourceList::iterator it= self->begin();
      for(int j(0); j!=i; ++j, ++it);
      return (*it); // note address of
   }
   size_t __len__() {      return self->size();       }
   
   void append( pointlike::Source& s){
       self->push_back(s); // note makes a copy
   }
}
//fails %include $(POINTLIKEROOT)/pointlike/ParamOptimization.h
%include pointlike/Alignment.h

%template(SourceVector)  std::vector<pointlike::Source>;
%template(StringDoubleMap) std::map<std::string,std::vector<double> >;
// these attempts fail
//%template(SkyDirIntPairVector) std::vector<std::pair<astro::SkyDir,int> >;

%include pointlike/ConfidenceLevel.h



