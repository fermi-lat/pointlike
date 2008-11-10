"""  spectral fit interface class SpectralAnalysis
    
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointspec.py,v 1.12 2008/11/05 00:41:13 kerrm Exp $

"""
version='$Revision: 1.12 $'.split()[1]
import os
from numpy import *

class SpectralAnalysis(object):
    """ 
    Interface to Matthew's spectral analysis code

    """
    
    def __init__(self,  event_files, history_files, diffuse_file, **kwargs):
        """
        call signature::

  sa = SpectralAnalysis(event_files, history_files, diffuse_file, **kwargs)


Create a new spectral analysis wrapper instance.  

    event_files: a list of event files (FT1) to process (or a single file). 
    history_files: a list of spacecraft history files (FT2) to process (or a single file).
    diffusefile:  Full path to a galactic diffuse file, like  GP_gamma_conventional.fits  (note isotropic)                                 

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  roi_dir     [None] direction to use if exposure calculation will be limited; otherwise use all sky                                                                                  
  roi_radius  [ 25]  radius (deg) to use if calculate exposure or a ROI. (180 for full sky)                                                                                
  livetimefile [None] Exposure file: if specified and not found, create from FT2/FT1 info                                                                                 
  
  datafile    [None]  HEALpix data file: if specified and not found, create from FT1 info                                                                               
  event_class  [-1] Select event class (-1: all, 0: front, 1:back)                                                                                  
  class_level [ 3]  select class level (set 0 for gtobssim                                                             
  emin        [100] Minimum energy                                                                                 
  emax        [None] Maximum energy: if not specified no limit                                                                                
  binsperdecade [4] When generating Bands from the FT1 data.
  use_mc_psf  [False] Use PSF determined by MC analysis of true; otherwise as defined by data

  isotropic   [(1.5e-5,2.1)] tuple of flux>100 MeV, spectral index for isotropic diffuse to add to diffuse
                                                                           
  extended_likelihood [False] Use extended likelihood                                                                         
  CALDB       [None] If not specified, will use environment variable
  irf         ['P6_v1_diff'] Used for effective area                                                                         
  quiet       [False] Set True to suppress (some) output                                                                               
  maxROI      [25] maximum ROI for PointSourceLikelihood to use
  =========   =======================================================
  """

        self.event_files = event_files if type(event_files)==type([]) else [event_files] 
        self.history_files=history_files if type(history_files)==type([]) else [history_files]
        self.roi_dir     = None
        self.roi_radius  = 25
        self.livetimefile= None
        self.datafile    = None
        self.align       = True
        self.binsperdecade=4
        self.use_mc_psf  = False
         
        self.diffusefile = diffuse_file
        self.isotropic   = (1.5e-5,2.1) 
         
        self.emin        = 100
        self.emax        = None
        self.extended_likelihood=False
        self.event_class  = -1 
        self.CALDB       = None 
        self.irf         ='P6_v1_diff'
        self.quiet       = False
        self.verbose     = False

        self.class_level = 3  # select class level
        self.maxROI      = 25 # for PointSourceLikelihood
        
        self.__dict__.update(kwargs) 

        # check explicit files
        for filename in  [self.event_files[0], self.history_files[0], self.diffusefile, self.CALDB] :
            if filename is not None and not os.path.exists(filename):
                raise Exception('SpectralAnalysis setup: file name or path "%s" not found'%filename)
        self.setup()

    class WrapExposure(object):
        """ Wrap up the C++ Exposure class to support an front/back distinction."""
        def __init__(self, exposure):
           self.exposure = exposure;

        def value(self, sdir, energy, event_class):
           return self.exposure[event_class].value(sdir, energy)

        
    def setup(self):
        from pointlike import  SimpleLikelihood
        #from Response import ModelResponse
        #from Fitters import EnergyBands

        # process livetime files
        self.lt = self.get_livetime()
 
        # process data files
        self.data = self.get_data()
        self.dmap = self.data.map()

        # setup exposure and background
        self.set_background()

        # for fitting
        SimpleLikelihood.enable_extended_likelihood(self.extended_likelihood) 

        #select bands to fit
        self.energies = fromiter((band.emin() for band in self.dmap if band.event_class()==0\
                  and band.emin()>=self.emin
                  and (self.emax is None or band.emax()<self.emax) ),float)
        print "selected energies: " , floor(self.energies+0.5)
        
    def get_livetime(self,  zenithcut=105,  pixelsize=1.0):
        from skymaps import LivetimeCube, Gti, SkyDir
        if self.livetimefile is None or not os.path.exists(self.livetimefile):
            if self.roi_dir is None:
                # no roi specified: use full sky
                self.roi_dir = SkyDir(0,0)
                self.roi_radius = 180
            lt = LivetimeCube(
                cone_angle=self.roi_radius,\
                dir=self.roi_dir,\
                pixelsize=pixelsize,\
                zcut=cos(radians(zenithcut)),\
                quiet=self.quiet )
            for n, hist in enumerate(self.history_files):
                lt.load(hist,  Gti(self.event_files[n]))
            if self.livetimefile is not None: lt.write(self.livetimefile)
        else: 
            lt = LivetimeCube(self.livetimefile)
            print 'loaded LivetimeCube %s ' % self.livetimefile
        return lt
 
    def get_data(self):
        from skymaps import PhotonBinner
        from pointlike import Data

        Data.set_class_level(self.class_level)

        if self.datafile is None or not os.path.exists(self.datafile):
            self.binner = PhotonBinner(self.binsperdecade) # per decade
            Data.setPhotonBinner(self.binner)
            if self.align: 
                # boresight alignment, only needed for August
                Data.set_alignment('default')
                for h in self.history_files:
                    Data.setHistoryFile(h)

            data = Data(self.event_files)
            if self.datafile is not None:
                print 'saving datafile %s for subsequent use' % self.datafile
                data.write(self.datafile)
        else:
            data = Data(self.datafile)
            print 'loaded datafile %s ' % self.datafile

        # modify the psf parameters in the band objects, which SimpleLikelihood will then use
        from psf import PSF
        import math
        self.psf = PSF(use_mc=self.use_mc_psf)
        if not self.quiet: print 'setting PSF parameters (use_mc=%d):\n  energy class  gamma sigma(deg)'%self.use_mc_psf
        for band in data.map():
             e = (band.emin()*band.emax())**0.5
             cl = band.event_class()
             gamma = self.psf.gamma(e,cl)
             sigma = self.psf.sigma(e,cl)
             if not self.quiet: print '%6.0f%5d%10.1f%10.2f' % (e, cl, gamma, sigma)
             band.setGamma(gamma)
             band.setSigma(math.radians(sigma))

        
        return data
        #self.data.info()
      

    def set_background(self):
        from skymaps import Exposure, EffectiveArea, Background,\
            DiffuseFunction, CompositeSkySpectrum, IsotropicPowerLaw, HealpixDiffuseFunc
        from pointlike import PointSourceLikelihood
        PointSourceLikelihood.set_maxROI(self.maxROI)
        if self.CALDB is not None: EffectiveArea.set_CALDB(self.CALDB)
        inst = ['front', 'back']
        self.ea  = [EffectiveArea(self.irf+'_'+x) for x in inst]
        if not self.quiet: print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        self.exposure = [Exposure(self.lt,ea) for ea in self.ea]

        import pyfits
        q = pyfits.open(self.diffusefile, memmap=1)
        if q[1].name == 'SKYMAP': # first extension name: is it healpix?
            self.galactic_diffuse = HealpixDiffuseFunc(self.diffusefile)
        else:
            self.galactic_diffuse = DiffuseFunction(self.diffusefile)
        self.isotropic_diffuse = IsotropicPowerLaw(self.isotropic[0],self.isotropic[1])
        self.diffuse = CompositeSkySpectrum(self.galactic_diffuse);
        self.diffuse.add(self.isotropic_diffuse)
        self.background = Background(self.diffuse, self.exposure[0], self.exposure[1]) # array interface does not work
        PointSourceLikelihood.set_background(self.background)
   

    class Fitter(object):
        """ manage spectral fitting"""

        def __init__(self, globaldata, name, src_dir):
            """

            """
            from pointlike import PointSourceLikelihood
            from wrappers  import PointSourceLikelihoodWrapper

            # create PointSourceLikelihood object and do fits
            PointSourceLikelihood.set_energy_range(globaldata.emin) #kluge
            self.src_dir = src_dir
            self.psl = PointSourceLikelihood(globaldata.dmap, name, self.src_dir)
            print 'TS= %6.1f' % self.psl.maximize()
            self.pslw = PointSourceLikelihoodWrapper(self.psl,globaldata.WrapExposure(globaldata.exposure),\
                        emin=globaldata.emin,quiet=globaldata.quiet)
            self.models = []
            
        def printSpectrum(self):
            self.psl.printSpectrum()

        def add_source(self, other):
            """ add another source to the background for this one """
            self.psl.addBackgroundPointSource(other.psl)
            self.pslw.update()

        def fit(self, model='PowerLaw',**kwargs):
            """ model: one of ['PowerLaw', 'BrokenPowerLaw', ...]
            """
            exec('from Models import %s'%model)
            exec('self.models += [self.pslw.poisson(%s(**kwargs))]'%model)

        def plot(self, fignum=1, date_tag=True, sed=True,e_weight=2,cgs=False):
            from wrappers import PointspecPlotter
            from pylab import figure,axes
            figure(fignum)
            mode = 'sed' if sed else 'counts'
            p = PointspecPlotter(self.pslw,e_weight=e_weight,cgs=cgs,models=self.models,mode=mode)
            axes(p.ax)
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
