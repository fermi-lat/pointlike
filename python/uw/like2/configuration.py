"""
Manage the analysis configuration

$Header$

"""
import os, sys, types
import numpy as np
from . import dataset
import skymaps
from uw.utilities import keyword_options
from uw.like import pypsf

        
class ExposureManager(object):
    """A small class to handle the trivial combination of effective area and livetime.
    
    Also handles an ad-hoc exposure correction
    """

    def __init__(self, dataset, **datadict): 
        """
        Parameters
        ----------
        dataset :  DataSet object
            for CALDB, aeff, some parameters
            
        datadict['exposure-correction'] : list of strings defining functions of energy
            the correction factors to apply to front, back 
        """

        def make_exposure():
            if dataset.exposure_cube is not None:
                ## use pregenerated gtexpcube2 cube; turn off interpolation
                return [skymaps.DiffuseFunction(f,1000.,False) for f in dataset.exposure_cube]
                
            skymaps.EffectiveArea.set_CALDB(dataset.CALDBManager.CALDB)
            skymaps.Exposure.set_cutoff(np.cos(np.radians(dataset.thetacut)))
            inst = ['front', 'back']
            aeff_files = dataset.CALDBManager.get_aeff()
            ok = [os.path.exists(file) for file in aeff_files]
            if not all(ok):
                raise DataSetError('one of CALDB aeff files not found: %s' %aeff_files)
            self.ea  = [skymaps.EffectiveArea('', file) for file in aeff_files]
            if dataset.verbose: print ' -->effective areas at 1 GeV: ', \
                    ['%s: %6.1f'% (inst[i],self.ea[i](1000)) for i in range(len(inst))]
            
            if dataset.use_weighted_livetime and hasattr(dataset, 'weighted_lt'):
                return [skymaps.Exposure(dataset.lt,dataset.weighted_lt,ea) for ea in self.ea]
            else:
                return  [skymaps.Exposure(dataset.lt,ea) for ea in self.ea]
                
        self.exposure = make_exposure()

        correction = datadict.pop('exposure_correction', None)
        if correction is not None:
            self.correction = map(eval, correction)
            energies = [100, 1000, 10000]
            if not self.quiet:  'Exposure correction: for energies %s ' % energies
            for i,f in enumerate(self.correction):
                if not self.quiet:  ('\tfront:','\tback: ','\tdfront:', '\tdback')[i], map( f , energies)
        else:
            self.correction = lambda x: 1.0, lambda x: 1.0
            
    def value(self, sdir, energy, event_type):
        return self.exposure[event_type].value(sdir, energy)*self.correction[event_type](energy)
        
    def __call__(self, event_type):
        """Return a SkySppectrum-compatible object of the exposure for the given event tupe (e.g., front or back)
        """
        class Exposure(object):
            def __init__(self, eman, event_type):
                self.eman = eman
                self.et =event_type
                self.energy=1000.
            def __repr__(self):
                return 'Exposure SkySpectrum for event_type %d, default energy set to %.0f MeV' % (self.et, self.energy)
            def __call__(self, sdir, e=None):
                if e is None: e=self.energy
                return self.eman.value(sdir, e, self.et)
            def setEnergy(self, e):
                self.energy=e
        return Exposure(self, event_type)
        
class PSFmanager(object):
    """ manage the PSF
    This is an interface to the original pointlike code, to hide some of that complexity,
    providing a clean interface. Especially, it is the only place that the event_type is handled
    """
    def __init__(self, dataset):
        """dataset : DataSet object, used to extract CALDBmanager
        """
        cdbman = dataset.CALDBManager 
        self.cdb_psf = pypsf.CALDBPsf(cdbman)
        self.irfname = cdbman.irf
        
    def __repr__(self):
        return '%s.%s: IRF "%s"' % (self.__module__, self.__class__.__name__, self.irfname)
        
    
    def __call__(self, event_type):
        """Return a PSF functor of distance for the given event type (e.g., front or back)
            the energy is specifed by a setEnergy method
        """
        cdb_psf = self.cdb_psf
        
        class PSF(object):
            def __init__(self,  event_type):
                self.event_type =event_type
                self.energy=1000.
            def __repr__(self):
                return 'PSF for event_type %d, current energy  %.0f MeV' % (self.event_type, self.energy)
            def setEnergy(self, e):
                self.energy=e
 
            def __call__(self, delta):
                class Band(object):
                    def __init__(self, energy, event_type, **kwargs):
                        self.e, self.ct = energy, event_type
                        self.__dict__.update(kwargs)
                        
                bpsf = pypsf.BandCALDBPsf(cdb_psf, Band(self.energy, self.event_type),
                        override_en=False,adjust_mean=False)
                return bpsf(delta)
               
        return PSF(event_type)

        
class ExposureCorrection(object):
    """ logarithmic interpolation function
    """
    def __init__(self, a,b, ea=100, eb=300):
        self.c = (b-a)/np.log(eb/ea)
        self.d =  a -self.c*np.log(ea)
        self.a, self.b = a,b
        self.ea,self.eb = ea,eb
    def __call__(self, e):
        if e>self.eb: return self.b
        if e<self.ea: return self.a
        return self.c*np.log(e) + self.d
    def plot(self, ax=None, **kwargs):
        import pylab as plt
        if ax is None: 
            ax = plt.gca()
        dom = np.logspace(1.5, 3.5, 501) 
        ax.plot(dom, map(self, dom), lw=2, color='r', **kwargs)
        plt.setp(ax, xscale='log', xlim=(dom[0], dom[-1]))
 

class Configuration(object):
    """
    Manage an analysis configuration, including the IRF and exposure and data
    """
    defaults =(
        ('irf', None,  'Set to override value in config.txt : may find in custom_irf_dir'),
        ('extended', None, 'Set to override value in config.txt'),
        ('nocreate', True, 'Do not allow creation of a binned photon file'),
        ('quiet', False, 'set to suppress most output'),
        ('postpone', False, 'set true to postpone loading data'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, modeldir='.', **kwargs):
        """ 
        parameters
        ----------
        modeldir: folder containing configuration definition, file config.txt
        """
        keyword_options.process(self, kwargs)
        if not self.quiet:  'Defining configuration in folder: %s', os.path.join(os.getcwd(),modeldir)
        
        # extract parameters config
        
        config = eval(open(os.path.expandvars(modeldir+'/config.txt')).read())
        for key in 'extended irf'.split():
            if self.__dict__[key] is None or self.__dict__[key]=='None': 
                if not self.quiet:
                    print '%s: %s set from config: "%s"' %(key, kwargs.get(key,None), config.get(key,None))
                self.__dict__[key]=config.get(key, None)

        # set up IRF
        
        irf = config['irf']
        if 'CUSTOM_IRF_DIR' not in os.environ and os.path.exists(os.path.expandvars('$FERMI/custom_irfs')):
            os.environ['CUSTOM_IRF_DIR'] = os.path.expandvars('$FERMI/custom_irfs')
        custom_irf_dir = os.environ['CUSTOM_IRF_DIR']
        
        # set up dataset
        
        datadict = config.get('datadict',None)
        assert isinstance(datadict, dict), 'Did not find datadict dictionary in config.txt'
        self.dataset = dataset.DataSet(datadict['dataname'], interval=datadict.get('interval',None),
                nocreate = self.nocreate,
                irf = irf,
                quiet = self.quiet,
                postpone=self.postpone,
                )
        if not self.quiet:  print self.dataset
        
        # use dataset to extract psf and exposure, set up respective managers
        
        exposure_correction=datadict.pop('exposure_correction', None)        
        self.exposureman = ExposureManager(self.dataset, exposure_correction=exposure_correction)
        
        self.psfman = PSFmanager(self.dataset)
        
  

