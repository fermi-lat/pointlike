"""  spectral fit interface class SpectralAnalysis
    
    $Header$

"""
import os
from numpy import *

class SpectralAnalysis(object):
    """ 
    Interface to Matthew's spectral analysis code
    """
    
    
    def __init__(self,  event_files, history_files,  **kwargs):

        self.event_files, self.history_files = event_files, history_files
        self.roi_dir     =None
        self.roi_radius  = 25
        self.livetimefile=None
        self.datafile    =None 
        self.diffusefilename=r'f:\glast\data\galprop\GP_gamma_conventional.fits'
        self.emin        = 100
        self.emax        = None
        self.binsperdecade=4
        self.method      = 'MP'
        self.extended_likelihood=False
        self.event_class  =-1 
        self.CALDB = 'f:\\glast\\packages\\ScienceTools-v9r7p1\\irfs\\caldb\\v0r7p1\\CALDB\\data\\glast\\lat'
        self.irf         ='P6_v1_diff'
        self.quiet       = False
        self.class_level = 3  # select class level
        self.maxROI      =25 # for PointSourceLikelihood
        self.__dict__.update(kwargs) 

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
        from skymaps import LivetimeCube, Gti
        if self.livetimefile is None or not os.path.exists(self.livetimefile):
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

        self.diffuse = DiffuseFunction(self.diffusefilename)
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
            from spectrum import Source
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
