"""
manage band classes

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/bands.py,v 1.8 2014/03/12 15:52:09 burnett Exp $
"""
import os
import numpy as np
import skymaps

energybins = np.logspace(2,5.5,15) # default 100 MeV to 3.16 GeV, 4/decade
  
class EnergyBand(object):
    """ Combine three concepts:
    * definition of the binning in solid angle and energy
    * PSF and exposure for the given configuration, appropriate to this band
    * The pixel data extracted from the binned photon data
    """
    def __init__(self, config,  roi_dir, 
            radius=5, event_type=1, emin=10**2, emax=10**2.25):
        """
        
        """
        # define bin boundaries

        self.skydir=roi_dir 
        self.radius =radius
        self.event_type = event_type
        self.emin, self.emax = emin, emax
        self.energy = energy =  np.sqrt(emin*emax)
        
        # save appropriate psf and exposure
        
        self.psf=config.psfman(event_type, energy)
        self.exposure = config.exposureman(event_type, energy)
        
        # used by client to integrate a function of energy over exposure
        self.integrator = self.exposure.integrator(self.skydir, self.emin, self.emax) 
  
        # these changed if data loaded -- see load_data
        self.pixel_area=0
        self.wsdl = None
        self.pix_counts=[]

    def __repr__(self):
        return '%s.%s: %s' % (self.__module__,self.__class__.__name__, self.title)
    def set_energy(self, energy):
        self.psf.setEnergy(energy)
        self.exposure.setEnergy(energy)
        self.energy=energy
    def load_data(self, cband):
        # cband is a C++ skymaps.Band object
        self.wsdl = skymaps.WeightedSkyDirList(cband, self.skydir, self.radius_in_rad, False)
        self.pix_counts = np.asarray([x.weight() for x in self.wsdl]) if len(self.wsdl)>0 else []
        self.pixel_area = cband.pixelArea()
    @property
    def radius_in_rad(self): return np.radians(self.radius)
    @property #alias, for compatibilty, but deprecated
    def sd(self): return self.skydir
    @property
    def solid_angle(self):
        return np.pi*self.radius_in_rad**2 
    @property
    def has_pixels(self): return self.wsdl is not None and self.pixels>0
    @property 
    def pixels(self):
        return len(self.wsdl) if self.wsdl is not None else 0
    @property
    def events(self):
        return sum(w.weight() for w in self.wsdl) if self.wsdl is not None else 0
    @property
    def title(self):
        ret = '%.0f-%.0f MeV event_type %d' % (self.emin, self.emax, self.event_type)
        if self.wsdl is not None:
            ret += ', %d pixels, %d events' % (self.pixels, self.events)
        return ret
        
class BandSet(list):
    """ manage the list of EnergyBand objects
    """
    def __init__(self, config, roi_index, load=False, radius=5):
        """create a list of energy bands
        
        config : configuration.Configuration object
        roi_index : int | None
            specify direction, either as nside=12 index, or from config.roi_spec
        load : bool
            if False, do not load data into the pixels
        """
        self.config = config
        if roi_index is None or roi_index<0:
            self.roi_dir = config.roi_spec.pos
            self.radius = config.roi_spec.radius
            roi_index = None
        else:
            ### HEALPix case
            self.roi_index = roi_index
            self.roi_dir = skymaps.Band(12).dir(roi_index) # could be defined otherwise
            self.radius=radius
        for emin, emax  in zip(energybins[:-1], energybins[1:]):
            for et in range(2):
                self.append(EnergyBand(config, self.roi_dir,  event_type=et, radius=self.radius, emin=emin,emax=emax))
        self.has_data = False
        
        if load:
            self.load_data()
    
    def __repr__(self):
        ret = '%s.%s : %d bands %d-%d MeV for ROI %d' % (
            self.__module__,self.__class__.__name__,
            len(self), self[0].emin, self[-1].emax, self.roi_index)
        if self.has_data:
            ret += ', %d pixels, %d events' % (self.pixels, self.events)
        return ret
            
    def find(self, emin, event_type):
        for band in self:
            if round(emin)==round(band.emin) and event_type==band.event_type :
                return band
        raise Exception('Band with emin, event_type = %d,%d not found in %s' % (round(emin), event_type, self))
    
    def load_data(self):
        """load the current dataset, have the respective bands fill pixels from it
        """
        dset = self.config.dataset
        dset.load()
        found = 0
        for i,cband in enumerate(dset.dmap):
            emin, emax, event_type =  cband.emin(), cband.emax(), cband.event_class()
            try:
                band = self.find(emin, event_type)
                found +=1
            except:
                continue
            band.load_data(cband)
        if found!=len(self):
            raise Exception('%s: Did not load all bands' % self.__repr__())
        self.has_data = True

    @property
    def pixels(self): return sum( band.pixels for band in self) if self.has_data else 0
    @property
    def events(self): return sum( band.events for band in self) if self.has_data else 0
        
