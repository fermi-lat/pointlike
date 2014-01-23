"""
Source classes
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sources.py,v 1.41 2013/12/27 16:09:58 burnett Exp $

"""
import os, copy
import numpy as np
from skymaps import SkyDir
from uw.like import  Models
from . import response

# convenience adapters 
def LogParabola(*pars, **kw):return Models.LogParabola(p=pars, **kw)
def PowerLaw(*pars, **kw):   return Models.PowerLaw(p=pars, **kw)
def ExpCutoff(*pars, **kw):  return Models.ExpCutoff(p=pars, **kw)
def PLSuperExpCutoff(*pars, **kw): return Models.PLSuperExpCutoff(p=pars, **kw)
def Constant(*pars, **kw):   return Models.Constant(p=pars, **kw)
def FBconstant(f,b, **kw): return Models.FrontBackConstant(f,b, **kw)
    
def ismodel(model):
    """ check that model is an instance of Models.Model"""
    return isinstance(model, Models.Model)

def set_default_bounds( model, force=False):
    """
    Handy utility to set bounds for a model from like.Models
    force=True to override previously set bounds.
    """
    if not force and hasattr(model, 'bounds'):
        # model has bounds. Were they set? check to see if all are None
        notset =  np.all(np.array([np.all(b ==[None,None]) for b in model.bounds]))
        if not notset: return
    bounds=[]
    def to_internal(fun, values):
        return [fun(value) if value is not None else None for value in values]
    for pname, mp in zip(model.param_names, model.mappers):
        plim = (None,None)
        try:
            plim = dict(
                Index=(-0.5, 5), 
                Norm=(10**-16, 10**-7),
                Scale=(0.001, 4.0),
                beta=(-0.1, 5.), 
                Cutoff=(100., 1e5),
                )[pname.split('_')[0]]
        except: pass
        bounds.append( to_internal(mp.tointernal, plim) )
    model.bounds = np.array(bounds) # convert to array so can mask with free

class Source(object):
    """ base class for all pointlike/like2 sources
    Subclasses are:
        PointSource
        ExtendedSource
        GlobalSource
        
    All instances have the folloiwng properties:
    * model, a Models.Model object
    * skydir : [skymaps.Skydir  | None]
        
    Subclasses must implement a function response(band), which, given a BandLite parameter, 
        returns a Response object appropriate for the source. This provides the angular dependence 
        of the response specific the band energy and event type.
    
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.changed = False # flag for bandlike
        assert self.name is not None, 'bad source name'
        self.name = str(self.name) # force to be a string
        if self.skydir is None:
            # global source: keep original model
            self.free = self.model.free.copy()  # save copy of initial free array to restore
            return
        elif hasattr(self.skydir, '__iter__'): #allow a tuple of (ra,dec)
            self.skydir = SkyDir(*self.skydir)
        if 'model' not in kwargs or self.model is None:
            self.model=LogParabola(1e-14, 2.2, 1e-3, 1e3)
            self.model.free[2:]=False
        self.free = self.model.free.copy()
        if self.model.name=='PowerLaw':
            par,sig = self.model.statistical()
            self.model = LogParabola(*(list(par)+[0, self.model.e0]))
            self.model.free[2:]=False
 
        elif self.model.name=='ExpCutoff':
            try:
                print 'converting %s to PLSuperExpCutoff' %self.name
                self.model = self.model.create_super_cutoff()
            except FloatingPointError:
                print 'Failed'
                
        elif self.model.name=='PowerLawFlux':
            f, gamma = self.model.get_all_parameters() #10**self.model.p
            emin = self.model.emin
            try:
                self.model=LogParabola(f*(gamma-1)/emin, gamma, 0, emin)
            except Exception, msg:
                print 'Failed to create LogParabola for source %s, pars= %s'% (self.name, (f,gamma,emin))
                raise
            self.model.free[2:]=False
        elif self.model.name=='LogParabola':
            if hasattr(self, 'free') and len(self.free)>3: self.free[3]=False
            self.model.free[-1]=False # make sure Ebreak is always frozen (bug in handling extended sources?)
            #assert not self.model.free[3], 'LogParabola model for source %s hae Ebreak free' % self.name
            if False: ##########TURN OFF FOR NOW not self.model.free[2] and self.model[2]!=0.:
                oldbeta=self.model[2]
                self.model[2]=0
                self.model.internal_cov_matrix[2,:]=0
                self.model.internal_cov_matrix[:,2]=0
                print 'Warning: set fixed beta=0 for source %s, beta was %.2f' %(self.name, oldbeta)
        if self.model.name not in ['LogParabola','PLSuperExpCutoff','ExpCutoff', 'Constant']:
            raise Exception('model %s not supported' % self.model.name)
        if not hasattr(self.model, 'npar'):
            raise Exception('model %s for source %s was not converted to new format'\
                    % (self.model.name, self.name))
        # finally, add bounds to the models object, ignoring similar capability in Models.
        set_default_bounds( self.model )
           
            
    def get_spectral_model(self):
        return self.model
    def set_spectral_model(self, newmodel):
        t =self.model
        self.model = newmodel
        return t
    spectral_model = property(get_spectral_model, set_spectral_model)

    def freeze(self, parname, value=None):
        self.model.freeze(parname)
        if value is not None: self.model.setp(parname, value)
        self.changed=True
        assert sum(self.model.free)>0, 'cannot freeze all parameters this way'

    def thaw(self, parname):
        self.model.freeze(parname, freeze=False)
        self.changed = True

    def __str__(self):
        return self.name + ' '+ self.skydir.__str__() +' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')
    def __repr__(self):
        return '%s.%s: %s' % (self.__module__,self.__class__.__name__ , self.name)
        
    @property
    def isextended(self):
        return hasattr(self, 'dmodel') and not self.isglobal

    @property
    def isglobal(self):
        return self.skydir is None

class PointSource(Source):
    def __init__(self, **kwargs):
        kwargs.update(spatial_model=None) # allow test for extent (no extent!)
        super(PointSource, self).__init__(**kwargs)
    def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)
    def copy(self):
        """ return a new PointSource object, with a copy of the model, others"""
        ret = PointSource(**self.__dict__)
        ret.model = self.model.copy()
        return ret
    def response(self, band, **kwargs):
        return response.PointResponse(self, band, **kwargs)


class ExtendedSource(Source):

    def __str__(self):
        return self.name + ' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')  
  
    def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)
        
    def copy(self):
        """ return a new ExtendSource object, with a copy of the model object"""
        ret = ExtendedSource(**self.__dict__)
        ret.model = self.model.copy()
        if ret.model.name=='LogParabola':
            ret.model.free[-1]=False # make sure Ebreak is frozen
        return ret
         
    def response(self, band, **kwargs):
        """ return a Respose object, which, given a band, can create a convolved image
        and calculate expected counts
        """
        return response.ExtendedResponse(self, band, **kwargs)

        
class GlobalSource(Source):
    def __init__(self, **kwargs):
        super(GlobalSource, self).__init__(**kwargs)
        self.dmodel= kwargs.get('dmodel', None)
        assert self.skydir is None # used as a flag

    def copy(self):
        """ return a new PointSource object, with a copy of the model, others"""
        ret = GlobalSource(**self.__dict__)
        ret.model = self.model.copy()
        return ret

    def response(self, band, **kwargs):
        """ return a Response class for the band"""
        assert self.dmodel, 'Need DiffuseBase object to determine response'
        try:
            resp_class =  dict(
                Isotropic =response.IsotropicResponse,
                MapCube   =response.DiffuseResponse,
                CachedMapCube=response.CachedDiffuseResponse,
                Healpix   =response.DiffuseResponse,
                HealpixCube = response.DiffuseResponse,
                IsotropicSpectralFunction = response.IsotropicResponse,
                )[self.dmodel.type]
        except Exception, msg:
            raise Exception('Could not find a response class for source %s:"%s"' %(self,msg))
        return resp_class(self,band, **kwargs) 
    

