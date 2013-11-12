"""
Manage the analysis configuration

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/configuration.py,v 1.6 2013/11/10 20:05:08 burnett Exp $

"""
import os, sys, types
import numpy as np
from . import ( dataset, exposure, psf)
import skymaps
from uw.utilities import keyword_options
        

class Configuration(object):
    """
    Manage an analysis configuration: the IRF, exposure and data
    """
    defaults =(
        ('event_class_names', ('front', 'back'), 'names to describe event classes'),
        ('bins_per_decade', 4, 'number of bins per decade'),
        ('irf', None,  'Set to override value in config.txt : may find in custom_irf_dir'),
        ('extended', None, 'Set to override value in config.txt'),
        ('nocreate', True, 'Do not allow creation of a binned photon file'),
        ('quiet', False, 'set to suppress most output'),
        ('postpone', False, 'set true to postpone loading data'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, configdir='.', **kwargs):
        """ 
        parameters
        ----------
        configdir: folder containing configuration definition, file config.txt
        It is a Python dictionary, and must have the keys:
            irf : a text string name for the IRF
            diffuse: a dict defining the diffuse components
            datadict: a dict with at least a key 'dataname'
        """
        keyword_options.process(self, kwargs)
        # make absolute path
        if configdir.startswith('$'):
            self.configdir = os.path.expandvars(configdir)
        elif not configdir.startswith('/'):
            self.configdir = os.path.join(os.getcwd(), configdir)
        else:
            self.configdir = configdir
        
        if not os.path.exists(self.configdir):
            raise Exception( 'Configuration folder %s not found' % self.configdir)
        if not self.quiet: print 'Defining configuration in folder: %s' % self.configdir
        
        # extract parameters config
        config = eval(open(os.path.expandvars(self.configdir+'/config.txt')).read())
        for key in 'extended irf'.split():
            if self.__dict__.get(key, None) is None: 
                if not self.quiet:
                    print '%s : override config.txt, "%s" => "%s"' %(key, kwargs.get(key,None), config.get(key,None))
                self.__dict__[key]=config.get(key, None)

        # set up IRF
        
        irf = config['irf']
        if 'CUSTOM_IRF_DIR' not in os.environ and os.path.exists(os.path.expandvars('$FERMI/custom_irfs')):
            os.environ['CUSTOM_IRF_DIR'] = os.path.expandvars('$FERMI/custom_irfs')
        custom_irf_dir = os.environ['CUSTOM_IRF_DIR']
        
        # define diffuse dictionary for constucting spatial models
        
        self.diffuse = config['diffuse']
        
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
        self.exposureman = exposure.ExposureManager(self.dataset, exposure_correction=exposure_correction)
        
        self.psfman = psf.PSFmanager(self.dataset)
      
    def __repr__(self):
        return '%s.%s: %s' %(self.__module__, self.__class__.__name__, self.configdir)
        
    def get_bands_for_ROI(self, index):
        return self.get_bands(skymaps.Band(12).dir(index), minROI=5)
    
    def get_bands(self, roi_dir, **kwargs):
    
        """ return a list of dataset.ROIband objexts for an ROI
        """
        dset = self.dataset
        dset.load()
        band_kwargs = dict(emin=dset.emin, emax=dset.emax, minROI=dset.minROI, maxROI=dset.maxROI)
        band_kwargs.update(kwargs)
        radius = band_kwargs['minROI'] # fixed radius now
        bandlist = []
        for band in dset.dmap:
            emin,emax, event_type =  band.emin(), band.emax(), band.event_class()&5
            if (emin + 1) < band_kwargs['emin'] or (emax - 1) >band_kwargs['emax']: continue
            #print int(emin), event_class
            energy= np.sqrt(emin*emax)
            bandlist.append( ROIBand(band, self.psfman(event_type,energy), self.exposureman(event_type,energy), 
                roi_dir, radius))
        return np.asarray(bandlist)
 
class ROIBand(object):
    """Data access class
        indludes the appropriate PSF and exposure functions for this energy band
    """

    def __init__(self, cband, psf, exposure,  skydir, radius):
        """ cband : SWIG-wrapped skymaps::Band C++ object
        """
    
        self.psf, self.exposure, self.skydir = psf, exposure, skydir
        self.radius       = radius
        self.event_type   = cband.event_class()
        self.emin         = cband.emin()
        self.emax         = cband.emax()
        self.energy       = 0
        self.set_energy(np.sqrt(self.emin * self.emax))
        self.integrator   = self.exposure.integrator(self.skydir, self.emin, self.emax)
        #self.exposureint  = self.exposure(self.sd, self.energy)*(self.emax-self.emin) ### FIXME
        ## data stuff
        self.wsdl = skymaps.WeightedSkyDirList(cband, skydir, self.radius_in_rad,False)
        self.pix_counts   = np.asarray([x.weight() for x in self.wsdl]) if len(self.wsdl)>0 else []
        self.npix         = self.wsdl.total_pix()
        self.has_pixels   = self.wsdl.counts()>0
        self.pixel_area   = cband.pixelArea()
        self.solid_angle  = self.npix*self.pixel_area #ragged edge
    @property
    def radius_in_rad(self):
        return np.radians(self.radius) 
    @property
    def sd(self): # deprecated
        return self.skydir
        
    def __repr__(self):
        return '%s.%s: %s, %d photons' \
            % (self.__module__,self.__class__.__name__, self.title, np.sum(self.pix_counts))
        
    def set_energy(self, energy):
        """set the energy for the psf and exposure objects
        """
        if energy==self.energy: return
        self.energy=energy
        self.psf.setEnergy(energy)
        self.exposure.setEnergy(energy)
    def solid_angle(self):
        return np.pi*self.radius_in_rad**2 
    @property
    def title(self):
        return '%.0f-%.0f MeV event_type %d' % (self.emin, self.emax, self.event_type)


 
class Bandlite(object):
    """ behaves like ROIBand, but does not require data
    Default initialization is for back, first energy bin above 100 MeV
    """
    def __init__(self, roi_dir, config,  roi_index=None  , 
            radius=5, event_type=1, emin=10**2, emax=10**2.25):
        self.skydir=roi_dir if roi_dir is not None else  skymaps.Band(12).dir(roi_index)
        self.radius =5
        self.event_type = event_type
        self.emin, self.emax = emin, emax
        self.energy = energy =  np.sqrt(emin*emax)
        self.psf=config.psfman(event_type, energy)
        self.exposure = config.exposureman(event_type, energy)
        self.has_pixels=False
        self.integrator = self.exposure.integrator(self.skydir, self.emin, self.emax)

    def __repr__(self):
        return '%s.%s: %s' \
            % (self.__module__,self.__class__.__name__, self.title)
    def set_energy(self, energy):
        self.psf.setEnergy(energy)
        self.exposure.setEnergy(energy)
        self.energy=energy
    def load_data(self, cband):
        self.wsdl = skymaps.WeightedSkyDirList(cband, skydir, self.radius_in_rad, False)
        self.has_pixels = len(self.wsdl)>0
    @property
    def radius_in_rad(self): return np.radians(self.radius)
    @property #alias, for compatibilty, but deprecated
    def sd(self): return self.skydir
    @property
    def solid_angle(self):
        return np.pi*self.radius_in_rad**2 
    @property
    def title(self):
        return '%.0f-%.0f MeV event_type %d' % (self.emin, self.emax, self.event_type)

class BandSet(list):
    """ manage the list of bands
    """
    def __init__(self, config, roi_index, max=None):
        """fill the list at the nside=12 indes"""
        self.config = config
        energybins = np.logspace(2,5.5,15) # default 100 MeV to 3.16 GeV
        self.roi_index = roi_index
        for emin, emax  in zip(energybins[:-1], energybins[1:]):
            for et in range(2):
                self.append(Bandlite(None, config, roi_index, event_type=et, emin=emin,emax=emax))
    
    def __repr__(self):
        return '%s.%s : %d bands %d-%d MeV for ROI %d' % (
            self.__module__,self.__class__.__name__,
            len(self), self[0].emin, self[-1].emax, self.roi_index)
    
    def load_data(self):
        """load the current dataset, have the respective bands fill pixels from it"""
        dset = self.config.dataset
        dset.load()
        for cband in dset.dmap:
            emin, emax, event_type =  cband.emin(), cband.emax(), cband.event_class()
            print emin, emax, event_type


        