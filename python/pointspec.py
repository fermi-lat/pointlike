"""  spectral fit interface class SpectralAnalysis
    
    $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/pointspec.py,v 1.3 2008/10/29 21:10:24 burnett Exp $

"""
version='$Revision$'
import os
from numpy import *

class SpectralAnalysis(object):
    """ 
    Interface to Matthew's spectral analysis code


    """
    
    
    def __init__(self,  event_files, history_files, diffuse_file, **kwargs):
        """
        call signature::

  sa = SpectralAnalysis(event_files, history_files,  **kwargs)


Create a new spectral analysis wrapper instance.  

    event_files: a list of event files (FT1) to process (or a single file)
    history_files: a list of spacecraft history files (FT2) to process (or a single file).
    diffusefile:  Full path to a diffuse file, like  GP_gamma_conventional.fits                                   

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
  
  method      ['MP'] Fit method                                                                                
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
        self.diffusefile = diffuse_file
        self.emin        = 100
        self.emax        = None
        self.binsperdecade=4
        self.method      = 'MP'
        self.extended_likelihood=False
        self.event_class  = -1 
        self.CALDB       = None #'f:\\glast\\packages\\ScienceTools-v9r7p1\\irfs\\caldb\\v0r7p1\\CALDB\\data\\glast\\lat'
        self.irf         ='P6_v1_diff'
        self.quiet       = False
        self.class_level = 3  # select class level
        self.maxROI      = 25 # for PointSourceLikelihood
        self.__dict__.update(kwargs) 

        # check explicit files
        for filename in  [self.event_files[0], self.history_files[0], self.diffusefile, self.CALDB] :
            if filename is not None and not os.path.exists(filename):
                raise Exception('SpectralAnalysis setup: file name or path "%s" not found'%filename)
        self.setup()

    class WrapExposure(object):
        """ create exposure required by GlobalData usage"""
        def __init__(self, exposure):
            self.exposure = exposure;

        def value(self, sdir, energy, event_class):
            return self.exposure[event_class].value(sdir, energy)

        def __call__(self, sdir, energies, event_class):
            if event_class==-1:
                return append(fromiter((self.value(sdir,e,0) for e in energies),float),\
                              fromiter((self.value(sdir,e,1) for e in energies),float))
            return fromiter((self.value(sdir,e,event_class) for e in energies),float)
        
    def setup(self):
        from pointlike import  SimpleLikelihood
        from Response import ModelResponse
        from Fitters import EnergyBands

        # process livetime files
        self.lt = self.get_livetime()
 
        # process data files
        self.data = self.get_data()
        self.dmap = self.data.map()

        # setup exposure and background
        self.set_background()

        # for fitting
        SimpleLikelihood.enable_extended_likelihood(self.extended_likelihood) 

        # required to implement GlobalData use by Source
        #select bands to fit
        self.energies = fromiter((band.emin() for band in self.dmap if band.event_class()==0\
                  and band.emin()>=self.emin
                  and (self.emax is None or band.emax()<self.emax) ),float)
        print "selected energies: " , floor(self.energies+0.5)

        self.sources=[] # needed for unknown reasons
        self.ebands=EnergyBands( self.energies )
        self.response = ModelResponse(self.ebands, SpectralAnalysis.WrapExposure(self.exposure))
        self.mask = lambda : [True]*(len(self.energies))*2 # for front and back
        self.bands= lambda : len(self.energies)

        
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
            data = Data(self.event_files)
            if self.datafile is not None:
                print 'saving datafile %s for subsequent use' % self.datafile
                data.write(self.datafile)
        else:
            data = Data(self.datafile)
            print 'loaded datafile %s ' % self.datafile
        return data
        #self.data.info()
        

    def set_background(self):
        from skymaps import Exposure, EffectiveArea, Background, DiffuseFunction
        from pointlike import PointSourceLikelihood
        PointSourceLikelihood.set_maxROI(self.maxROI)
        if self.CALDB is not None: EffectiveArea.set_CALDB(self.CALDB)
        inst = ['front', 'back']
        self.ea  = [EffectiveArea(self.irf+'_'+x) for x in inst]
        print ' -->effective areas at 1 GeV: ', ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
        self.exposure = [Exposure(self.lt,ea) for ea in self.ea]

        self.diffuse = DiffuseFunction(self.diffusefile)
        self.background = Background(self.diffuse, self.exposure[0], self.exposure[1]) # array interface does not work
        PointSourceLikelihood.set_background(self.background)
   

    class Fitter(object):
        """ generate spectrum"""

        def __init__(self, globaldata, name, src_dir):
            from pointlike import PointSourceLikelihood

            # create PointSourceLikelihood object and do fits
            PointSourceLikelihood.set_energy_range(globaldata.emin) #kluge
            self.src_dir = src_dir
            self.psl = PointSourceLikelihood(globaldata.dmap, name, self.src_dir)
            print 'TS= %6.1f' % self.psl.maximize()
            self.globaldata = globaldata
            
        def printSpectrum(self):
            self.psl.printSpectrum()

        def add_source(self, other):
            """ add another source to the background for this one """
            self.psl.addBackgroundPointSource(other.psl)

        def fit(self, **kwargs):
            """ model: one of ['PowerLaw', 'BrokenPowerLaw', ...]
            """
            from SourceLib import Source
            self.fitter = Source(self.psl, self.globaldata)
            self.fitter.fit(**kwargs)

        def plot(self, fignum=30, date_tag=True,
                flag='sed', cgs=False, sedweight=2  # relevant if sed only.
                ):
            """ spectral plot
            
            """
            import SourcePlot
            import fermitime
            #sed(self.fitter, models = self.fitter.MP.models)
            models = self.fitter.MP.models
            SourcePlot.spectrum(self.fitter,models=models,flag=flag,fignum=fignum,cgs=cgs,sedweight=sedweight)
            if date_tag: fermitime.date_tag()

    def fitter(self, name, source_dir):
        return SpectralAnalysis.Fitter(self, name, source_dir)
