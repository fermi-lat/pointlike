"""  A module to provide simple and standard access to pointlike fitting and spectral analysis.  The
     relevant parameters are fully described in the docstring of the constructor of the SpectralAnalysis
     class.
    
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointspec.py,v 1.14 2008/11/17 21:37:48 kerrm Exp $

    author: Matthew Kerr
"""
version='$Revision: 1.14 $'.split()[1]
import os
from numpy import array, arange
from pixeldata import PixelData

class SpectralAnalysis(object):
    """ 
    Interface to the spectral analysis code

    """
    
    def __init__(self,  event_files, history_files, diffuse_file, **kwargs):
        """
        call signature::

  sa = SpectralAnalysis(event_files, history_files, diffuse_file, **kwargs)


Create a new spectral analysis object.  

    event_files: a list of event files (FT1) to process (or a single file). 
    history_files: a list of spacecraft history files (FT2) to process (or a single file).
    diffusefile:  Full path to a galactic diffuse file, like  GP_gamma_conventional.fits  (note isotropic)                                 

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================

  =========   KEYWORDS CONTROLLING DATA BINNING AND LIVETIME CALCULATION
  roi_dir     [None] direction to use if exposure calculation will be limited; otherwise use all sky                                                                                  
  roi_radius  [ 25]  radius (deg) to use if calculate exposure or a ROI. (180 for full sky)                                                                                
  livetimefile [None] Exposure file: if specified and not found, create from FT2/FT1 info  
  zenithcut   [105]  Maximum spacecraft pointing angle with respect to zenith to allow
  align       [False] if True, perform boresigh alignment correction on pre-September flight data  
  datafile    [None]  HEALpix data file: if specified and not found, create from FT1 info                                                                               
  class_level [ 3]  select class level (set 0 for gtobssim                                                             
  binsperdecade [4] When generating Bands from the FT1 data.
  use_mc_psf  [False] Use PSF determined by MC analysis if true; otherwise as defined by data

  =========   KEYWORDS CONTROLLING EXPOSURE CALCULTION
  CALDB       [None] If not specified, will use environment variable
  irf         ['P6_v1_diff'] Used for effective area
  
  =========   KEYWORDS CONTROLLING BACKGROUND
  isotropic   [(1.5e-5,2.1)] tuple of flux>100 MeV, spectral index for isotropic diffuse to add to diffuse
  galactic_scale [1.] scale factor for galactic diffuse
  
  =========   KEYWORDS CONTROLLING SPECTRAL ANALYSIS
  event_class  [-1] Select event class (-1: all, 0: front, 1:back)                                                       
  emin        [100] Minimum energy                                                                                 
  emax        [None] Maximum energy: if not specified no limit
  extended_likelihood [False] Use extended likelihood
  maxROI      [25] maximum ROI for PointSourceLikelihood to use
  minROI      [None] minimum ROI

  =========   KEYWORDS OF MISCELLANY
  quiet       [False] Set True to suppress (some) output
  verbose     [FALSE] More output
  =========   =======================================================
  """
        self.quiet       = False
        self.verbose     = False
        self.emin        = 100
        self.emax        = None
        self.extended_likelihood=False
        self.event_class  = -1        
        self.maxROI      = 25 # for PointSourceLikelihood
        self.minROI     = 0
        self.__dict__.update(kwargs)

        self.pixeldata =  PixelData(event_files,history_files,**kwargs)
        from wrappers import ExposureWrapper,BackgroundWrapper
        self.exposure  =  ExposureWrapper(self.pixeldata,**kwargs)
        self.background = BackgroundWrapper(diffuse_file,self.exposure,**kwargs)        
        self.setup()

        
    def setup(self):
        from pointlike import SimpleLikelihood,PointSourceLikelihood
        PointSourceLikelihood.set_maxROI(self.maxROI)
        PointSourceLikelihood.set_minROI(self.minROI)
        SimpleLikelihood.enable_extended_likelihood(self.extended_likelihood) 
        PointSourceLikelihood.set_background(self.background())
        
        #select bands to fit
        #self.energies = fromiter((band.emin() for band in self.dmap if band.event_class()==0\
        #          and band.emin()>=self.emin
        #          and (self.emax is None or band.emax()<self.emax) ),float)
        #print "selected energies: " , floor(self.energies+0.5)

    class Fitter(object):
        """ manage spectral fitting"""

        def __init__(self, spectralanalysis, name, src_dir):
            """

            """
            from pointlike import PointSourceLikelihood
            from wrappers  import PointSourceLikelihoodWrapper

            # create PointSourceLikelihood object and do fits
            PointSourceLikelihood.set_energy_range(spectralanalysis.emin) #kluge
            self.src_dir = src_dir
            self.psl = PointSourceLikelihood(spectralanalysis.pixeldata.dmap, name, self.src_dir)
            self.psl.maximize()
            self.pslw = PointSourceLikelihoodWrapper(self.psl,spectralanalysis.exposure,**spectralanalysis.__dict__)
            print 'TS= %6.1f' % self.pslw.TS() #TS consistent with energy selection for spectral analysis
            self.models = []
            
        def printSpectrum(self):
            """Display results from spatial (and/or extended) likelihood analysis."""
            self.psl.printSpectrum()

        def add_source(self, other):
            """Add another source to the background of the one in question."""
            self.psl.addBackgroundPointSource(other.psl)
            self.pslw.update()

        def add_models(self,models,**kwargs):
            """models: a list of spectral models, e.g. ['PowerLaw,'BrokenPowerLaw',...]
               A list of available models is available via Fitter.available_models.
               Keyword arguments are passed to the model constructor."""
            for model in models:
                if model in self.available_models(printout=False):
                    if model not in [m.name for m in self.models]:
                        exec('from Models import %s'%model)
                        exec('self.models += [%s(**kwargs)]'%model)
                else:
                    print 'Warning! %s not a valid spectral model.'%model
                    self.available_models()

        def available_models(self,printout=True):
            """List available spectral models.
               printout - if False, return a list of the names"""
            
            from Models import DefaultModelValues
            if printout:
               print 'Available spectral models:'
               for name in DefaultModelValues.names: print name
            else: return DefaultModelValues.names
           

        def fit(self, models=['PowerLaw'],**kwargs):
            """ Fit all models that have been added, including via the optional keyword argument here.
                models: a list of spectral models, e.g. ['PowerLaw,'BrokenPowerLaw',...]
            """
            from specfitter import SpectralModelFitter
            self.add_models(models,**kwargs)
            for model in self.models:
               if not model.good_fit: SpectralModelFitter.poisson(self.pslw,model)

        def plot(self, fignum=1, date_tag=True, **kwargs):
            """Plot the flux or counts spectrum for the point source.      

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  fignum      [1] the label for the matplotlib figure
  date_tag    [True] place the date and time on the plot
  sed         [True] if True, plot the flux density; else, plot the counts spectrum
  e_weight    [2] factor of energy by which to multiply dN/dE
  cgs         [False] boolean; True to use ergs for abscissa energy; False uses ordinate units
  mev         [True] boolean; True to use MeV for ordinate units; False to use GeV
  =========   =======================================================
            """

            from plotters import EnergyAnalysis
            from pylab import figure,axes
            figure(fignum)
            p = EnergyAnalysis(self.pslw,models=self.models,**kwargs)
            p.show()
            if date_tag:
               import fermitime
               fermitime.date_tag()


    def fitter(self, name, source_dir):
        """
        return a SpecralAnalysis.Fitter object

        name        name to use for output
        source_dir  SkyDir object specifying where to start


        """
        return SpectralAnalysis.Fitter(self, name, source_dir)
