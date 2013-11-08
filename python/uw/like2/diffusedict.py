"""
Manage the diffuse sources

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/diffusedict.py,v 1.8 2013/10/29 14:09:17 burnett Exp $

author:  Toby Burnett
"""
import os, types, pyfits, collections, zipfile, pickle
import numpy as np
import pandas as pd 
from .pub import healpix_map #make this local in future
from uw.utilities import keyword_options, convolution
from . import convolution

import skymaps #from Science Tools: for SkyDir, DiffuseFunction, IsotropicSpectrum
from uw.like import  Models
def PowerLaw(*pars, **kw):   return Models.PowerLaw(p=pars, **kw)

class DiffuseException(Exception):pass


class DiffuseBase(object):
    """Base class for global diffuse sources
    expect subclasses to implement SkySpectrum interface
    """
    def setupfile(self, filename):
        """ filename: string 
        """
        self.filename = os.path.expandvars(filename)
        if not os.path.exists(self.filename):
            self.filename = os.path.expandvars(os.path.join('$FERMI','diffuse',self.filename))
        assert os.path.exists(self.filename), 'DiffuseFunction file "%s" not found' % self.filename
        self.loaded =  False

    def load(self): 
        assert not hasattr(self, 'filename'), 'Logic error; class %s must implement load' %self.name
        pass
    
    def __repr__(self):
        t= self.__class__.__name__ + ': '+self.filename
        if hasattr(self, 'kw') and self.kw is not None:
            t+= '\n'+ '\n'.join(['\t    %-10s: %s' % item for item in self.kw.items() if item[0]!='filename'])
        return t
    
    def grid_generator(self, band):# skydir, psf, exposure):
        """Return a GridGenerator object appropriste for this model, and the psf, exposre
         """
        return convolution.GridGenerator(band, self) #psf, exposure, skydir, self)
        
    @property
    def name(self):
        return self.__class__.__name__
        
    def show(self, title=None, scale='log', **kwargs):
        """make an AIT image for testing
        """
        from uw.utilities import image
        ait = image.AIT(self, **kwargs)
        ait.imshow(title=self.name if title is None else title, scale=scale)
        return ait.axes.figure
       

class Isotropic(DiffuseBase, skymaps.IsotropicSpectrum):
    """Implement the isotropic diffuse, by wrapping skymaps.IsotopicSpectrum
    """
    def __init__(self, filename):
        self.filename=filename
        super(Isotropic, self).__init__(filename)
        self.loaded=True
        
    def load(self):
        pass


class MapCube(DiffuseBase, skymaps.DiffuseFunction):
    """ wrapper for eventual invocation of skymaps.DiffuseFunction, which interprets a map cube
    load() must be called before use
    """
    def __init__(self, filename):
        """ filename: string or dict
        """
        self.setupfile( filename)
            
    def load(self, interpolate=False):
        if  self.loaded: return
        self.loaded=True
        if not interpolate: 
            print 'loading diffuse file %s: warning, not interpolating' %self.filename
        super(MapCube,self).__init__(self.filename, 1000., interpolate)

class Healpix(DiffuseBase):
    """Diffuse map using HEALPix representation.
    Presumes that columns have eneregies (found in hdu#3) which exactly
    correspond to bands
    """
    def __init__(self, filename):
        self.setupfile(filename)
        self.energy=0
        self.eindex=-1
        self.fits = None
        
    def load(self):
        if self.fits is not None: return
        self.fits = pyfits.open(self.filename)
        table = self.fits[1]
        self.columns = table.data
        self.col_names = [t.name for t in table.get_coldefs()]
        self.energies = self.fits[2].data.field('MeV')

    def setEnergy(self, energy):
        self.load()
        if energy==self.energy: return
        self.eindex = -1
        for i,e in enumerate(self.energies):
            if abs(e-energy)<0.01*e:
                self.eindex= i
                self.energy=energy
                break
        assert self.eindex>=0, 'Energy %.0f not found' % energy
        cname = self.col_names[self.eindex]
        self.hpm = healpix_map.HParray(cname, self.columns.field(cname))
        
    def __call__(self, skydir, energy=None):
        if energy is not None and energy !=self.energy: self.setEnergy(energy)
        return self.hpm(skydir)
        
    @property
    def name(self):
        return '%s %s' % (self.__class__.__name__, 
            self.col_names[self.eindex] if self.eindex>-1 else 'not loaded')
    
class IsotropicSpectralFunction(DiffuseBase):
    """ wrapper for using a standard spectral function with isotropic
    """
    def __init__(self, expression):
        """
        expression: 
        """
        try:
            self.expression  = expression.split('_')[-1]
            self.spectral_function = eval(self.expression)
        except Exception, msg:
            print 'Failure to evaluate IsotropicSpectralFunction %s : %s' % (self.expression, msg)
            raise
        self.energy=1000
    def __repr__(self):
        return '%s: %s' % (self.__class__.__name__, self.expression )
    def __call__(self, skydir, energy=None):
        return self.spectral_function(self.energy if energy is None else energy)
    def setEnergy(self, energy): self.energy=energy


class CachedMapCube(DiffuseBase):
    """ for compatibility with previous models"""
    def __init__(self, zipfilename):
        self.setupfile(zipfilename)
    
    def load(self):
        z = zipfile.ZipFile(self.filename)
        t = z.namelist()
        if len(t)==1729: # if ziped with foldername
            t = t[1:]
        self.files = sorted(t) 
        assert len(self.files)==1728, 'wrong number of files: expected 1728, found %d' % len(files)
        self.opener = z.open

    def grid_generator(self, band): #skydir, psf, exposure):
        """Return a GridGenerator object
         """
                
        return convolution.CachedGridGenerator(band, band.sd, self) #psf, exposure, skydir,self, quiet=False, **self.kw)

    # these not implemented, and not needed for this class
    def __call__(self, skydir, energy=None):
        return np.nan
    def setEnergy(self, energy):
        pass

class DiffuseList(list):
    """A list of  subtype, or front/back list of DiffuseBase objects If only one, applied to all
    """
    def __init__(self, inlist):
        super(DiffuseList,self).__init__(inlist)
        self.gg = None
    def __repr__(self):
        return  '%s.%s\n\t' % (self.__module__, self.__class__.__name__) + '\n\t'.join([t.__repr__() for t in self])
        
    def __getitem__(self, index):
        if len(self)==1: index=0
        return super(DiffuseList,self).__getitem__(index) 
    def load(self):
        for x in self:
            x.load()
            
    def grid_generator(self, band): #skydir, event_type, configuration):
        """return a grid generator for the given skydir and event type index. Gets psf and exposure from config.
        """
        return self[band.event_type].grid_generator(band) #skydir, configuration.psfman(event_type), configuration.exposureman(event_type))
                
    def grid_for_band(self, band):
        """return convolved grid for an Band object, which has the event_type, position, psf, esposure, and energy
        """
        self.load()
        grid = self[band.event_type].grid_generator(band) #band.sd, band.psf, band.exposure)
        return grid(band.energy)
    
def diffuse_factory(value):
    """
    Create a DiffuseList object from a text specification
    """
    isdict = issubclass(value.__class__, dict)
    if  not hasattr(value, '__iter__') or isdict:
        value  = (value,)
    
    if isdict:
        files = DiffuseList([val['filename'] for val in value])
        type = value[0].pop('type', None)
        kws = value
    else:
        files = DiffuseList(value)
        kws = type = None

    # checking only first element, and only upto a comma, if not with '('
    f = files[0].split(',')[0] if files[0].find('(')<0 else files[0]
    ext = os.path.splitext(f)[-1]
    if type is not None:
        try:
            dfun = eval(type)
        except Exception, msg:
            raise Exception('Diffuse type specification "%s" failed: %s'%(type, msg))
    else:
        # inferr class to use from file type
        try:
            dfun = {'.txt': Isotropic, 
                '.fit': MapCube, '.fits': MapCube,
                '.zip': CachedMapCube,
                ')': IsotropicSpectralFunction, 
                }[ext if ext[-1]!=')' else ')']
        except Exception, msg:
            raise Exception('File type, "%s", for diffuse not recognized, from "%s":%s (message :%s)'\
                % (ext, files, ext, msg))
    
    if dfun==IsotropicSpectralFunction:
        diffuse_source = map(dfun,files)
    else:
        full_files = map( lambda f: os.path.expandvars(os.path.join('$FERMI','diffuse',f)), files)
        check = map(lambda f: os.path.exists(f) or f[-1]==')', full_files) 
        assert all(check), 'not all diffuse files %s found' % full_files
        diffuse_source= map(dfun, full_files) 
    # if a dict, add keywords to the objects
    if kws is not None:
        for x,kw in zip(diffuse_source, kws):
            x.kw = kw
    return DiffuseList(diffuse_source)


class DiffuseDict(collections.OrderedDict):
    """ create a dictionary of global diffuse objects
        key:   a string defined by the filename following an underscore, or key in input to init
        value: (both,) or (front,back) or (...) of  diffuse objects determined by the extension:
            txt: IsotropicSpectrum
            fit or fits: DiffuseFunction
            ): IsotropicSpectralFunction
            none: expect that the key is an expression to be evaluated to create a spectral model function
    
    Also maintain a list of spectral models for each ROI: the attribute models
    """
    def __init__(self, modeldir='.'):
        """
        modeldir : string 
            folder containing a file config.txt
        """
        try:
            self.spec=eval(open(os.path.join(modeldir, 'config.txt')).read()
                #.replace('dict(','collections.OrderedDict(')
                )['diffuse']
        except Exception, msg:
            print 'Failed to open model at %s: %s' % (modeldir, msg)
            raise
        super(DiffuseDict, self).__init__()
        for key, value in self.spec.items():
            #print key,value
            self[key] = diffuse_factory(value) 
        self.models=[dict() for i in range(1728)] # list of spectral models for each ROI

    def __repr__(self):
        r = self.__class__.__name__ + ':'
        for key, values in self.items():
            r += '\n    '+key
            for value in values:
                r += '\n\t' + str(value)
        return r
        
    def add_model(self, index, name, model):
        """ add a model to the dict with ROI index index
        """
        assert name in self.keys()
        model.background=True # for printing
        self.models[index][name]=model
        
    def get_sources(self, index, GlobalSource):
        """ return a list of GlobalSource objects for the given index
        """
        global_sources = []
        for name in self.keys():
            model = self.models[index].get(name, None)
            if model is None: continue
            s = GlobalSource(name=name, model=model, skydir=None, index=index) 
            s.dmodel = self[name]
            global_sources.append(s)
        return global_sources

def test(model_dir='.', energy=1333, event_type=0):
    """ Test creation of diffuse list from config file 
    """
    from skymaps import SkyDir
    from uw.like2 import configuration
    gc, npole = SkyDir(0,0, SkyDir.GALACTIC), SkyDir(0,90, SkyDir.GALACTIC)
    
    # setup the configuration, and load the diffuse dictioinary
    config = configuration.Configuration(model_dir, quiet=True, postpone=True)
    dd = DiffuseDict(model_dir)
    class Bandlite(object):
        def __init__(self, roi_dir, event_type=event_type, energy =energy):
            self.event_type=event_type
            self.energy=energy
            self.sd=roi_dir
            self.psf=config.psfman(event_type,energy)
            self.exposure = config.exposureman(event_type,energy)
            self.radius = 5
        def set_energy(self, energy):
            self.psf.setEnergy(energy)
            self.exposure.setEnergy(energy)
            
    band = Bandlite(roi_dir=gc)
   
    for dname, dlist in dd.items():
        dm = dlist[event_type]
        print '\n%-10s: %s' % (dname, dm)
        dm.load() # check load
        print '\tvalue at GC, NP: %.2e, %.2e'%  tuple([dm(dir, energy) for dir in gc, npole])
        # check convolution at GC
        gg = dlist.grid_generator(band)
        grid = gg(energy)
        cvals = grid.cvals
        cvmean = cvals.mean()
        print '\tconvolved grid: mean=%.1f counts, std/mean=%.3f' % ( cvmean, cvals.std()/cvmean if cvmean>0 else np.nan)

    
if __name__=='__main__':
    test()