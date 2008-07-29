// $Header: /nfs/slac/g/glast/ground/cvs/pointlike/src/mainpage.h,v 1.5 2007/12/03 00:39:40 burnett Exp $
// Mainpage for doxygen

/*! \mainpage package pointlike

   \author  Toby Burnett, Marshall Roth, Bruce Lesnick


  This package has classes implemeting a simple, fast point-source likelihood, 
  class pointlike::PointSourceLikelihood.

  The model it implements is a PSF (see pointlike::PsfFunction for the definition of the GLAST standard power-law function)
  plus a background. The background may be uniform, or will be taken from an instance of the class pointlike::DiffuseFunction.

  The data is encapsulated by the Data class, a healpix map of photon data, defined in skymaps::BinnedPhotonData, 
  created from a fits file or list of fits files, see the class pointlike::Data.

  The following applications are defined. Each uses embed_python::Module to parse a python file containing parameters
  
  - pointfit - fit to a specified list of sources
  - pointfind - search a specified region of the sky, finding all sources

  Additionally, a lot of code is in the python folder.

    <hr>
\section intro Discussion
Philosophy:

According to the Cramer-Rao bound, a lower bound on the variance of any unbiased estimator is given by maximum likelihood. This implies that the most sensitive detection of point sources, if the background is understood, is unbinned maximum likelihood.
At the UW, we have been interested in point source detection for some time, and previously developed a variable-width wavelet technique. We switched to likelihood when we realized how to overcome technical problems: speed of the non-linear fitting, and scaling to large data sets. Our solutions are

- binning: Binning means that one will not achieve the Cramer-Rao bound, but if the bin sizes are chosen to be smaller than features in the associated variable, the effect will be small. We use nested HEALpix bins, correlated with the energy. We typically bin the entire sky at once, saving the result as a file. This avoids the problems with the poles and variable-size pixels. The SC2 obssim data for a year, 16M photons, is 22 MB in this format; it took a fiew minutes to create.
- speed: We employ Newton-Raphson methods for the two types of optimizations, magnitude and position, using specialized code with analytic derivatives. Convergence is very fast.
- Combining front and back: We originally used front only, but the back contributes 1/2 the statistics, with a PSF only 1.6 or so worse, which is not a problem at high energy. We gained this statistical power by including back photons with front in bins with the same PSF. [ This obviously distorts the spectrum; we will have to account for, or just separate the bins when we get to that.]
The result is fast enough, that we can afford a brute-force search for sources, basically examining every pixel in the sky.

Feature Summary:

- Likelihood: single point source + arbitrary background
counts in upto 9 energy bands with 1-parameter optimization
position - 2-parameter optimization of the total likelihood
significance TS
PSF parameters from collection of bright sources (or allgamma)
includes dispersion

- Detection
Find sources, down to TS=10 in all sky (3 hours to detect 3200 in 1 year @ 10% spurious)

- LAT Alignment
Measure the rotation angles to a few arcmin in a day. 

- Measure PSF parameters

Class hierarchy: (see the astro and skymaps packages)

- astro::SkyFunction - abstract base class defining a real function on the sphere (a function of a astro::SkyDir)
- skymaps::SkyImage - implement SkyFunction with FITS image; also create from a SkyFunction
- skymaps::SkySpectrum - abstract, allow specification of energy spectrum at any point
- skymaps::DiffuseFunction - interpolate a FITS cube. Used for the background for point source fits.
- skymaps::PhotonMap - pixelized photon data (wraps map_tools::PhotonMap for now)
- pointlike::PointSourceLikelihood - represent a fit
- skymaps::Convolution - convolution of a SkySpectrum object with another SkySpectrum, perhaps a PSF.
- skymaps::CompositeSkySpectrum - linear combination of SkySpectrum objects. Used to combine the galactic diffuse with nearby (< 1deg) strong sources
- skymaps::Exposure - Integrate an exposure cube over the acceptance to define the exposure at any point.


\section notes release notes
  release.notes
\section requirements requirements
\include requirements

*/

