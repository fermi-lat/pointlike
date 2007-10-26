// $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/mainpage.h,v 1.2 2007/06/14 20:01:45 burnett Exp $
// Mainpage for doxygen

/*! \mainpage package pointlike

   \author  Toby Burnett, Marshall Roth, Bruce Lesnick


  This package has classes implemeting a simple, fast point-source likelihood, 
  class PointSourceLikelihood.

  The model it implements is a PSF (see pointlike::PsfFunction for the definition of the GLAST standard power-law function)
  plus a diffuse background. The background may be uniform, or will be taken from an instance of the class pointlike::DiffuseFunction.

  The data is encapsulated by the Data class, a healpix map of photon data, defined in map_tools::PhotonMap, 
  created from a fits file or list of fits files, see the class pointlike::Data.

  The following applications are defined. Each uses embed_python::Module to parse a python file containing parameters
  
  - pointfit - fit to a specified list of sources
  - pointfind - search a specified region of the sky, finding all sources

    <hr>
  \section notes release notes
  release.notes
  \section requirements requirements
  \include requirements

*/

