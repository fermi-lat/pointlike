"""  A module to provide simple and standard access to pointlike fitting and spectral analysis.  The
     relevant parameters are fully described in the docstring of the constructor of the SpectralAnalysis
     class.
    
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointspec.py,v 1.28 2009/07/28 13:01:20 burnett Exp $

    author: Matthew Kerr
"""
version='$Revision: 1.28 $'.split()[1]
import os
import sys

from numpy import array, arange
from pixeldata import PixelData

class AnalysisEnvironment(object):

   def init(self,glast_ts):

      from os.path import join

      if glast_ts:

         self.diffdir = r'f:/glast/data/galprop'
         self.catdir  = r'f:/glast/data/kerr'
         self.CALDB   = 'f:\\glast\\caldb\\v0r7p1\\CALDB\\data\\glast\\lat'

   def __init__(self,glast_ts=True,**kwargs):
      self.init(glast_ts)

class ConsistentBackground(object):

   def __init__(self,analysis_environment,background='ems_ring',**kwargs):

      from wrappers import Singleton
      from skymaps import DiffuseFunction,IsotropicSpectrum
      from psmanager import CatalogManager

      diffdir = analysis_environment.diffdir
      catdir  = analysis_environment.catdir

      self.background = background

      if background == 'nms_galprop':

         cat = os.path.join(catdir,r'gll_psc9month_v2r2.fit')

         singl = Singleton(DiffuseFunction,'gf',os.path.join(diffdir,'gas_mapcube_54_77Xvarh7S_P6_v3_diff_front.fits'))
         singl.add(        DiffuseFunction,'gb',os.path.join(diffdir,'gas_mapcube_54_77Xvarh7S_P6_v3_diff_back.fits'))
         singl.add(        DiffuseFunction,'if',os.path.join(diffdir,'ics_isotropic_mapcube_54_77Xvarh7S_P6_v3_diff_front.fits'))
         singl.add(        DiffuseFunction,'ib',os.path.join(diffdir,'ics_isotropic_mapcube_54_77Xvarh7S_P6_v3_diff_back.fits'))
         singl.add(        IsotropicSpectrum,'iso',os.path.join(diffdir,'Total_b30_EGBfree.txt'))

         gas_f  = singl('gf')
         gas_b  = singl('gb')
         ic_f   = singl('if')
         ic_b   = singl('ib')
         iso    = singl('iso')

         self.dmodels = [ [gas_f,gas_b] , iso , [ic_f,ic_b] ]

      if background == 'ems_ring':

         cat = os.path.join(catdir,r'gll_psc11month_v1b.fit')

         singl = Singleton(DiffuseFunction,'gf',os.path.join(diffdir,'gll_iem_v02_P6_v3_diff_front.fits'))
         singl.add(        DiffuseFunction,'gb',os.path.join(diffdir,'gll_iem_v02_P6_v3_diff_back.fits'))
         singl.add(        IsotropicSpectrum,'if',os.path.join(diffdir,'isotropic_iem_front_v02.txt'))
         singl.add(        IsotropicSpectrum,'ib',os.path.join(diffdir,'isotropic_iem_back_v02.txt'))

         gal_f  = singl('gf')
         gal_b  = singl('gb')
         iso_f  = singl('if')
         iso_b  = singl('ib')

         self.dmodels = [ [gal_f,gal_b], [iso_f,iso_b] ]
      
      singl.add(CatalogManager,'cat',cat)
      self.cm = singl('cat')

   def get_bgmodels(self, models = None, lat = None):

      from roi_modules import ROIBackgroundModel
      from Models import Constant,PowerLaw

      gal_index = True if (lat is not None and abs(lat) < 20) else False

      if self.background == 'nms_galprop':

         if models is not None:
            gas_m,iso_m,ic_m = models
         else:
            gas_m = PowerLaw(p=[1,1],free=[True,gal_index],index_offset=1)
            iso_m = Constant(free=[True])
            ic_m  = Constant(free=[False])            

         gas_model  = ROIBackgroundModel(self.dmodels[0],
                                          gas_m, 'Gas Galactic Diffuse')
         iso_model  = ROIBackgroundModel(self.dmodels[1],
                                          iso_m, 'Isotropic Diffuse')
         ic_model   = ROIBackgroundModel(self.dmodels[2],
                                          ic_m, 'IC Galactic Diffuse')
         return [gas_model,iso_model,ic_model]

      if self.background == 'ems_ring':

         if models is not None:
            gal_m,iso_m = models
         else:
            gal_m = PowerLaw(p=[1,1],free=[True,gal_index],index_offset=1)
            iso_m = Constant(free=[True])

         gal_model  = ROIBackgroundModel(self.dmodels[0],
                                          gal_m, 'Galactic Diffuse')
         iso_model  = ROIBackgroundModel(self.dmodels[1],
                                          iso_m, 'Isotropic Diffuse')
         return [gal_model,iso_model]

class SpectralAnalysis(object):
    """ 
    Interface to the spectral analysis code.

    """
    
    def __init__(self,  event_files, history_files, analysis_environment, **kwargs):
        """
        call signature::

  sa = SpectralAnalysis(event_files, history_files, analysis_environment, **kwargs)


Create a new spectral analysis object.  

    event_files: a list of event files (FT1) to process (or a single file). 
    history_files: a list of spacecraft history files (FT2) to process (or a single file).
    diffusefile:  Full path to a galactic diffuse file, like  GP_gamma_conventional.fits  (note isotropic)                                 

Optional keyword arguments:

  =========    =======================================================
  Keyword      Description
  =========    =======================================================

  =========    KEYWORDS CONTROLLING DATA BINNING AND LIVETIME CALCULATION
  +roi_dir     [None] direction to use if exposure calculation will be limited; otherwise use all sky                                                                                  
  +roi_radius  [ 25]  radius (deg) to use if calculate exposure or a ROI. (180 for full sky)                                                                                
  +livetimefile [None] Exposure file: if specified and not found, create from FT2/FT1 info  
  +zenithcut   [105]  Maximum spacecraft pointing angle with respect to zenith to allow
  +align       [False] if True, perform boresight alignment correction on pre-September flight data  
  +datafile    [None]  HEALpix data file: if specified and not found, create from FT1 info                                                                               
  +class_level [ 3]  select class level (set 0 for gtobssim)
  +mc_src_id   [-1] set to select on MC_SRC_ID column in FT1
  +use_mc_energy [False] set True to use MC_ENERGY instead of ENERGY
  +binsperdecade [4] When generating Bands from the FT1 data.
  +tstart      [0] Default no cut on time; otherwise, cut on MET > tstart
  +tstop       [0] Default no cut on time; otherwise, cut on MET < tstop

  =========    KEYWORDS CONTROLLING INSTRUMENT RESPONSE
  +exp_irf     ['P6_v3_diff'] Which IRF to use for effective area
  +psf_irf     ['P6_v3_diff'] Which IRF to use for the point spread function
  
  =========    KEYWORDS CONTROLLING BACKGROUND
  +background  ['ems_ring'] - a choice of global model specifying a diffuse background (Galactic, isotropic, etc.) and a catalog
  
  ### below is obsolete for ROI code -- only old pointliek spectral analysis
  +isotropic   [(1.5e-5,2.1)] tuple of flux>100 MeV, spectral index for isotropic diffuse to add to diffuse
  +iso_file    [None] if provided, a MapCube representing the isotropic diffuse; if not, use default power law
  +iso_scale   [1.] scale for the optional MapCube
  +galactic_scale [1.] scale factor for galactic diffuse
  
  =========    KEYWORDS CONTROLLING SPECTRAL ANALYSIS
  +event_class [-1] Select event class (-1: all, 0: front, 1:back)                                                       
  +emin        [100] Minimum energy                                                                                 
  +emax        [5e5] Maximum energy -- note, hardwired, so can only be _less_ than default!
  +extended_likelihood [False] Use extended likelihood
  +maxROI      [25] maximum ROI for PointSourceLikelihood to use
  +minROI      [0] minimum ROI
  +use_pointlike [False] if set to true, skips some steps for older analysis setup

  =========    KEYWORDS OF MISCELLANY
  +quiet       [False] Set True to suppress (some) output
  +verbose     [False] More output
  =========   =======================================================
  """
        self.quiet       = False
        self.verbose     = False
        self.exp_irf     = 'P6_v3_diff'
        self.psf_irf     = 'P6_v3_diff'
        self.background  = 'ems_ring'
        self.CALDB       = analysis_environment.CALDB
        self.roi_dir     = None
        self.ae          = analysis_environment
        self.use_pointlike = False
        if self.use_pointlike:
            # deprecated: if use this interface must set these
            self.event_class = -1        
            self.maxROI      = 25 # for PointSourceLikelihood
            self.minROI      = 0
            self.extended_likelihood=False
            self.emin        = 100
            self.emax        = 5e5



        self.__dict__.update(kwargs)


         #TODO -- sanity check that BinnedPhotonData agrees with analysis parameters
        self.pixeldata = PixelData(event_files,history_files,**self.__dict__)

        from wrappers import ExposureWrapper,BackgroundWrapper,Singleton
        self.exposure = ExposureWrapper(self.pixeldata,**self.__dict__)

        if self.use_pointlike:
            self.background = BackgroundWrapper(analysis_environment.galactic,self.exposure,**self.__dict__)
            self.setup()

        
    def setup(self):
        from pointlike import SimpleLikelihood,PointSourceLikelihood
        PointSourceLikelihood.set_maxROI(self.maxROI)
        PointSourceLikelihood.set_minROI(self.minROI)
        SimpleLikelihood.enable_extended_likelihood(self.extended_likelihood) 
        PointSourceLikelihood.set_background(self.background())

    #is passing a function name really what is wanted here??
    #def new_data(self,event_files,history_files):
    #    self.pixeldata = PixelData(event_files,history_files,self.__dict__.update)

    class Fitter(object):
        """ manage spectral fitting"""

        def __init__(self, spectralanalysis, name, src_dir, **kwargs):
            """
            spectral analysis: an instance of SpectralAnalysis
            name: an appellation for the source
            src_dir: an instance of SkyDir giving the source location

            """
            from pointlike import PointSourceLikelihood
            from wrappers  import PointSourceLikelihoodWrapper

            # create PointSourceLikelihood object and do fits
            PointSourceLikelihood.set_energy_range(spectralanalysis.emin) #kluge
            self.src_dir = src_dir
            self.psl =  PointSourceLikelihood(spectralanalysis.pixeldata.dmap, name, self.src_dir)
            self.pslw = PointSourceLikelihoodWrapper(self.psl,spectralanalysis.exposure,**kwargs.update(spectralanalysis.__dict__))
            print 'TS= %6.1f' % self.pslw.TS() #TS consistent with energy selection for spectral analysis
            self.models = []
            
        def printSpectrum(self):
            """Display results from spatial (and/or extended) likelihood analysis."""
            self.psl.printSpectrum()

        def localize(self):
            """Localize the point source and recalculate the likelihood."""
            err = self.psl.localize()
            self.pslw.update()
            return err

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
            if 'prefit' in kwargs:
               prefit = kwargs['prefit']
               kwargs.pop('prefit')
            else:
               prefit = True
            from specfitter import SpectralModelFitter
            if type(models) != type([]): models = [models]
            self.add_models(models,**kwargs)
            for model in self.models:
               if not model.good_fit: SpectralModelFitter.poisson(self.pslw,model,prefit=prefit)

        def plot(self, fignum=None, date_tag=True, filename=None,**kwargs):
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
            from pylab import figure,axes,savefig
            if type(fignum) is type(2): figure(fignum)
            p = EnergyAnalysis(self.pslw,models=self.models,**kwargs)
            p.show()
            if date_tag:
               import fermitime
               fermitime.date_tag()
            if filename: savefig(filename)


    def fitter(self, name, source_dir, **kwargs):
        """
        return a SpecralAnalysis.Fitter object

        name        name to use for output
        source_dir  SkyDir object specifying where to start


        """
        return SpectralAnalysis.Fitter(self, name, source_dir, **kwargs)
        
    def __call__(self, name, source_dir):
        """
        return a SpecralAnalysis.Fitter object, by simple call

        name        name to use for output
        source_dir  SkyDir object specifying where to start

        """
        return SpectralAnalysis.Fitter(self, name, source_dir)

    def roi(self, point_sources = None, backgrounds = None, previous_fit = None, **kwargs):
        """
        return an ROIAnalysis object with default settings.

        point_sources    a list of PointSource objects to merge with a Catalog list
                         (if None, the nearest catalog source will be fit)

        backgrounds      a list of ROIBackgroundModels with which to replace the default
                         isotropic and Galactic backgrounds (optional)
        """
        
        from roi_modules  import ROIPointSourceManager,ROIBackgroundManager
        from roi_analysis import ROIAnalysis
        
        lat,models,noneigh = None,None,False
        if 'lat'     in kwargs.keys(): lat = kwargs.pop('lat')
        if 'models'  in kwargs.keys(): models = kwargs.pop('models')
        if 'noneigh' in kwargs.keys(): noneigh = kwargs.pop('noneigh')

        cb = ConsistentBackground(self.ae,self.background)

        if backgrounds is None:
            backgrounds = cb.get_bgmodels(models=models,lat=lat)

        if noneigh:
            source_list = point_sources
        else:
            source_list = cb.cm.merge_lists(self.roi_dir,self.maxROI+5,point_sources)

        if previous_fit is not None:
            from roi_analysis import read_fit
            source_list,bg = read_fit(previous_fit)
            for i in xrange(len(bg)):
                backgrounds[i].smodel = bg[i]

        ps_manager = ROIPointSourceManager(source_list,self.roi_dir)
        bg_manager = ROIBackgroundManager(self, models = backgrounds)

        if point_sources is None: #make closest catalog source one to fit
            for i in xrange(len(ps_manager.models[0].free)):
               ps_manager.models[0].free[i] = True

        return ROIAnalysis(ps_manager,bg_manager,self,**kwargs)


#--------------------------------------------------------

def help(msg=None):
    if msg: print '\nError: %s' % msg
    print __doc__
    print SpectralAnalysis.__init__.__doc__
    sys.exit(0)

def main():
    
    from getopt import getopt, GetoptError
    from skymaps import SkyDir
    from glob import glob

    options = ''#'b:w:v'
    
    #get keyword arguments from docstring... maintenance will be easier
    long_options = [x[1:]+'=' for x in SpectralAnalysis.__init__.__doc__.split() if x[0] == '+']
    print long_options
    
    try:    
        (opts, args) = getopt(sys.argv[1:], options, long_options )
    except GetoptError, msg:
        help(msg)
    
    kwargs = dict()
    for opt in opts:
      exec 'kwargs[opt[0][2:]] = %s'%opt[1] in locals()

    ev_files,history_files,diffuse_file,sourcelist = args
    if ev_files[0]=='@':
        ev_files = [line.strip() for line in file(ev_files[1:]) if len(line)>0 and line[0]!='#']
    else: ev_files = glob(ev_files)
    if history_files[0]=='@':
        history_files = [line.strip() for line in file(history_files[1:]) if len(line)>0 and line[0]!='#']
    else: history_files = glob(history_files)

    sa = SpectralAnalysis(ev_files,history_files,diffuse_file,**kwargs)
    
    tokens = [x.strip().split() for x in file(sourcelist) if x.strip()[0] != '#']
    fitters = [0]*len(tokens)
    for n,token in enumerate(tokens):
        name = token[0]
        d = SkyDir(token[1],token[2])
        fitters[n] = sa.fitter(name,d)
        if len(token) > 3:
           models = token[3:]
           fitters[n].fit(models=models)

    

#--------------------------------------------------------
    
if __name__=='__main__':
    try:
        main()
    except:
        help()

