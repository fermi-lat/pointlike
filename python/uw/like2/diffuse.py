"""
Manage the diffuse sources

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/diffuse.py,v 1.37 2014/01/31 14:56:10 burnett Exp $

author:  Toby Burnett
"""
import os, types, pyfits, collections, zipfile, pickle
import numpy as np
import skymaps #from Science Tools: for SkyDir, DiffuseFunction, IsotropicSpectrum

from .pub import healpix_map #make this local in future
from . import (response, sources    )

def PowerLaw(*pars, **kw):   return sources.PowerLaw(*pars, **kw)

class DiffuseException(Exception):pass


class DiffuseBase(object):
    """Base class for global diffuse sources
    expect subclasses to implement SkySpectrum interface
    
    This uses a modified singleton pattern: 
    any subclass will have only one instance associated with the constructor argument, usually a filename.
    see http://stackoverflow.com/questions/42558/python-and-the-singleton-pattern
    """
    _instance_dict = dict()
    def __new__(cls, filename, **kwargs):
        if filename not in cls._instance_dict:
            cls._instance_dict[filename] = super(DiffuseBase, cls).__new__(
                                cls, filename, **kwargs)
            #print 'diffuse: loaded ', filename
        return cls._instance_dict[filename]
    
    def setupfile(self, filename):
        """ filename: string 
        """
        self.filename=filename
        self.fullfilename = os.path.expandvars(filename)
        if not os.path.exists(self.fullfilename):
            self.fullfilename = os.path.expandvars(os.path.join('$FERMI','diffuse',self.filename))
        assert os.path.exists(self.fullfilename), 'DiffuseFunction file "%s" not found' % self.fullfilename
        self.loaded =  False

    def load(self): 
        assert not hasattr(self, 'filename'), 'Logic error; class %s must implement load' %self.name
        pass
    
    def __repr__(self):
        t= self.__class__.__name__ + ': '+self.filename
        if hasattr(self, 'kw') and self.kw is not None:
            t+= '\n'+ '\n'.join(['\t    %-10s: %s' % item for item in self.kw.items() if item[0]!='filename'])
        return t
        
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
       
    def plot_spectra(self, ax=None, glat=0, glon=(0, 2,-2,  30, -30, 90, -90), title=None, label=None):
        """ show spectral for give glat, various glon values """
        from matplotlib import pylab as plt
        if not self.loaded: self.load()
        ee = np.logspace(1.5,6,101)
        if ax is None:
            fig, ax = plt.subplots(1,1, figsize=(5,5))
        else:
            fig = ax.figure
        for x in glon:
            sd = skymaps.SkyDir(glat,x, skymaps.SkyDir.GALACTIC)
            ax.loglog(ee, map(lambda e: e**2*self(sd,e), ee), 
                        label='b=%.0f'%x if label is None else label)
            if hasattr(self, 'energies'):
                et = self.energies
                ax.loglog(et, map(lambda e: e**2*self(sd,e), et), '+k')
        
        ax.grid(); 
        if len(glon)>1:
            ax.legend(loc='lower left', prop=dict(size=10))
        plt.setp(ax, xlabel=r'$\mathrm{Energy\ [MeV]}$', ylabel=r'$\mathrm{E^2\ dN/dE\ [Mev\ cm^{-2}\ s^{-1}\ sr^{-1}]}$')
        ax.set_title('Galactic diffuse at l=%.0f'%glat if title is None else title, size=12)
        return fig
        
    def plot_map(self, energy=1000, title=None, cbtext='', **kwargs):
        """show a map at the given energy"""
        self.setEnergy(energy)
        if not self.loaded: self.load()
        fig = self.show(cbtext='log10(flux)', **kwargs)
        fig.axes[0].set_title('%.0f MeV'%self.energy if title is None else title, size=12)
        return fig

class Isotropic(DiffuseBase): 
    """Implement the isotropic diffuse
    """
    def __init__(self, filename):
        self.setupfile(filename)
        try:
            t = np.loadtxt(self.fullfilename)
        except Exception,msg:
            raise Exception('Fail to load file %s: %s' % (self.fullfilename, msg))
        self.energies=t[:,0]
        self.spectrum = t[:,1]
        self.loaded=True
        self.logeratio = np.log(self.energies[1]/self.energies[0])
        self.setEnergy(1000.)

        
    def load(self):
        pass

    def setEnergy(self, energy): 
        # set up logarithmic interpolation
        self.energy=energy
        r = np.log(energy/self.energies[0])/self.logeratio
        self.energy_index = i = max(0, min(int(r), len(self.energies)-2))
        self.energy_interpolation = r-i

    def __call__(self, skydir, energy=None):
        if energy is not None and energy!=self.energy: 
            self.setEnergy(energy)
        a,i = self.energy_interpolation, self.energy_index
        return np.exp(    np.log(self.spectrum[i])   * (1-a) 
                        + np.log(self.spectrum[i+1]) * a     ) 



class MapCube(DiffuseBase, skymaps.DiffuseFunction):
    """ wrapper for eventual invocation of skymaps.DiffuseFunction, which interprets a map cube
    load() must be called before use
    """
    def __init__(self, filename):
        """ filename: string or dict
        """
        if not self.__dict__.get('loaded', False): #allows for invokation of singleton
            self.setupfile( filename)
            
    def load(self, interpolate=False):
        if  self.loaded: return
        self.loaded=True
        if not interpolate: 
            pass
            #print 'loading diffuse file %s: warning, not interpolating' %self.filename
        super(MapCube,self).__init__(self.filename, 1000., interpolate)

class HealpixCube(DiffuseBase):
    """ Jean-Marc's vector format, or the column version
    """
    def __init__(self,filename):
        """ filename : string
                Name of a FITS file, perhaps a gz. Either vector or column format
        """
        if not self.__dict__.get('loaded', False):
            self.setupfile( filename)
            
    def load(self):
        try:
            hdus = pyfits.open(self.fullfilename)
            self.energies = hdus[2].data.field(0)
            self.vector_mode = len(hdus[1].columns)==1
            if self.vector_mode:
                # one vector column, expect 2d array with shape (12*nside**2, len(energies))
                self.spectra = hdus[1].data.field(0)
                self.nside = int(np.sqrt(self.spectra.shape[0]/12.))
                assert self.spectra.shape[1]==len(self.energies), 'shape inconsistent with number of energies'
            else:
                # one column per energy: expect len(energies) columns
                hdu1 = hdus[1]
                assert len(hdu1.columns)==len(self.energies) , 'wrong number of columns'
                self.data = hdu1.data
                self.nside = int(np.sqrt(self.data.field(0).flatten().shape[0]/12.))
                
            self.loaded=True
            self.indexfun = skymaps.Band(self.nside).index
        except Exception, msg:
            print 'bad file or unexpected FITS format, file %s: %s' % (self.fullfilename, msg)
            raise
        self.logeratio = np.log(self.energies[1]/self.energies[0])
        self.setEnergy(1000.)
        
    
    def __call__(self, skydir, energy=None):
        if energy is not None and energy!=self.energy: 
            self.setEnergy(energy)
        skyindex = self.indexfun(skydir)
        a = self.energy_interpolation
        return np.exp( np.log(self.eplane1[skyindex]) * (1-a) 
                     + np.log(self.eplane2[skyindex]) * a      )


    def setEnergy(self, energy): 
        # set up logarithmic interpolation
        self.energy=energy
        r = np.log(energy/self.energies[0])/self.logeratio
        self.energy_index = i = max(0, min(int(r), len(self.energies)-2))
        self.energy_interpolation = r-i
        if self.vector_mode:
            self.eplane1 = self.spectra[:,i]
            self.eplane2 = self.spectra[:,i+1]
        else:
            self.eplane1 = np.ravel(self.data.field(i))
            self.eplane2 = np.ravel(self.data.field(i+1))
            
    
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
        self.col_names = [t.name for t in table.columns]
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
        self.loaded=True
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
        assert len(self.files)==1728, 'wrong number of files: expected 1728, found %d' % len(self.files)
        self.opener = z.open

    def __repr__(self):
        return '%s.%s: %s' % (self.__module__, self.__class__.__name__, self.filename)
    #def grid_generator(self, band): #skydir, psf, exposure):
        #"""Return a GridGenerator object
        # """
        #        
        #return convolution.CachedGridGenerator(band, band.sd, self) #psf, exposure, skydir,self, quiet=False, **self.kw)

    # these not implemented, and not needed for this class
    def __call__(self, skydir, energy=None):
        return np.nan
    def setEnergy(self, energy):
        pass

class DiffuseList(list):
    """A list of  event type list of DiffuseBase objects. If only one, applied to all
    """
    def __init__(self, inlist, event_type_names):
        super(DiffuseList,self).__init__(inlist)
        self.gg = None
        self.event_type_names = event_type_names
    def __repr__(self):
        return  '%s.%s\n\t' % (self.__module__, self.__class__.__name__) + '\n\t'.join([t.__repr__() for t in self])
        
    def __getitem__(self, index):
        if len(self)==1: index=0
        return super(DiffuseList,self).__getitem__(index) 
    def load(self):
        for x in self:
            x.load()
    @property
    def type(self):
        """return the diffuse class name implementing the global source"""
        return self[0].__class__.__name__
        
    def plot_spectra(self, glat=0, glon=(0, 2,-2,  30, -30, 90, -90), title=None):
        import pylab as plt
        fig,ax = plt.subplots(1,1, figsize=(5,5))
        for i,ds in enumerate(self):
            ds.plot_spectra(ax=ax, glat=glat, glon=glon, label=self.event_type_names[i], title=title)
        ax.grid()
        ax.legend(prop=dict(size=10))
        return fig

                
def diffuse_factory(value, event_type_names=('front', 'back')):
    """
    Create a DiffuseList object from a text specification
    value : [string | list | dict ]
        if string: a single filename to apply to all event types; if contains a '*', expand it
            to a list with the list of event type names, e.g., front, back
        if list: a set of filesnames corresponding to the list of event types
        if dict: must have keyword "filename", may have "type" to specify the type, 
            which must be the name of a class in this module inheriting from DiffuseBase
    
    If the keyword "type" is not specified, filenames are examined for extensions:
        txt : Isotropic
        zip : CachedMapCube
        fit for fits : MapCube
    A special case is ')' : IsotropicSpectralFunction
    """
    isdict = issubclass(value.__class__, dict)
    if isinstance(value, str) and  '*' in value:
        value = [value.replace('*',et) for et in event_type_names]
        
    elif  not hasattr(value, '__iter__') or isdict:
        value  = (value,)
    
    if isdict:
        try:
            files = DiffuseList([val['filename'] for val in value], event_type_names)
        except KeyError:
            raise DiffuseException('expected "filename" key in dict')
        type = value[0].pop('type', None)
        kws = value
    else:
        files = DiffuseList(value,event_type_names)
        kws = type = None

    # checking only first element, and only upto a comma, if not with '('
    f = files[0].split(',')[0] if files[0].find('(')<0 else files[0]
    ext = os.path.splitext(f)[-1]
    if type is not None:
        try:
            dfun = eval(type)
        except Exception, msg:
            raise DiffuseException('Diffuse type specification "%s" failed: %s'%(type, msg))
    else:
        # inferr class to use from file type
        try:
            dfun = {'.txt': Isotropic, 
                '.fit': MapCube, '.fits': MapCube,
                '.zip': CachedMapCube,
                ')': IsotropicSpectralFunction, 
                }[ext if ext[-1]!=')' else ')']
        except Exception, msg:
            raise DiffuseException('File type, "%s", for diffuse not recognized, from "%s":%s (message :%s)'\
                % (ext, files, ext, msg))
    
    if dfun==IsotropicSpectralFunction:
        diffuse_source = map(dfun,files)
    else:
        full_files = map( lambda f: os.path.expandvars(os.path.join('$FERMI','diffuse',f)), files)
        check = map(lambda f: os.path.exists(f) or f[-1]==')', full_files) 
        if not all(check):
            raise DiffuseException('not all diffuse files %s found' % full_files)
        diffuse_source= map(dfun, full_files) 
    # if a dict, add keywords to the objects
    if kws is not None:
        for x,kw in zip(diffuse_source, kws):
            x.kw = kw
    return DiffuseList(diffuse_source, event_type_names=event_type_names)


