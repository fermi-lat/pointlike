"""
Source descriptions for SkyModel
$Header$

"""
import os, pickle, glob, types, copy
import numpy as np
from skymaps import SkyDir,  IsotropicSpectrum, DiffuseFunction
from uw.like import  Models
from uw.like import pointspec_helpers


class Source(object):
    """ base class for various sources
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        if self.skydir is None:
            # global source: keep original model
            self.free = self.model.free.copy()  # save copy of initial free array to restore
            return
        
        if 'model' not in kwargs:
            self.model=Models.LogParabola(p=[1e-14, 2.2, 1e-3])
            self.model.free[2:]=False
        if self.model.name=='PowerLaw':
            par,sig = self.model.statistical()
            self.model = Models.LogParabola(p=list(par)+[1e-3, self.model.e0])
            self.model.free[2:]=False
        elif self.model.name=='PLSuperExpCutoff':
            par,sig=self.model.statistical()
            self.model = Models.ExpCutoff(p=par[:-1])
        elif self.model.name=='LogParabola':
            self.model.free[-1]=False
        elif self.model.name=='PowerLawFlux':
            f, gamma = 10**self.model.p
            emin = self.model.emin
            self.model=Models.LogParabola(p=[f*(gamma-1)/emin, gamma, 1e-3, emin])
            self.model.free[2:]=False
        if self.model.name not in ['LogParabola','ExpCutoff','Constant']:
            raise Exception('model %s not supported' % self.model.name)
            
        self.free = self.model.free.copy()  # save copy of initial free array to restore
    def freeze(self, freeze):
        self.model.free[:] = False if freeze else self.free
    def __str__(self):
        return self.name + ' '+ self.skydir.__str__() +' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')
 
class PointSource(Source):
    def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)
        
class GlobalSource(Source):
    def __init__(self, **kwargs):
        super(GlobalSource, self).__init__(**kwargs)
        assert self.skydir is None # used as a flag
    
class ExtendedSource(Source):
    def __str__(self):
        return self.name + ' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')    
    def near(self, otherdir, distance=10):
        return self.skydir.difference(otherdir) < np.radians(distance)

class Singleton(object):
    _instance={}
    def __init__(self,constructor):
        self.constructor = constructor
        self.key=str(constructor)
    def set_instance(self,  *pars):
        Singleton._instance[self.key]= self.constructor(*pars)
    def instance(self):
        try:
            return Singleton._instance[self.key]
        except KeyError:
            print 'SkyModel: Global source %s not initialized' % self.key
            raise
            
class Diffuse(Singleton):
    """ manage a skymaps.DiffuseFunction object, create only when new fits file found to open
    """
    _dfile = None
    locked = False
    key = None
    def __init__(self, dfile, lock=False):
        super(Diffuse,self).__init__(DiffuseFunction)
        if Diffuse.locked or dfile==Diffuse._dfile: return
        Diffuse._dfile=dfile
        self.set_instance(dfile)
        print 'loaded new diffuse map, %s with lock=%s' %(dfile, lock)
        Diffuse.locked = lock
        
class Isotropic(Singleton):
    """ manage a skymaps.IsotropicSpectrum object, create only when new fits file found to open
    """
    _dfile = None
    locked = False
    key = None
    def __init__(self, dfile, lock=False):
        super(Isotropic,self).__init__(IsotropicSpectrum)
        if Isotropic.locked or dfile==Isotropic._dfile: return
        Isotropic._dfile=dfile
        self.set_instance(dfile)
        print 'loaded new isotropic spectrum, %s, with lock=%s' %(dfile, lock)
        Isotropic.locked = lock

class ExtendedCatalog( pointspec_helpers.ExtendedSourceCatalog):
    """ subclass to add this lookup function """

    def __init__(self, *pars, **kwargs):
        """ initialize by also filling an array with all source spectral models"""
        super(ExtendedCatalog,self).__init__(*pars, **kwargs)
        self.sources = self.get_sources(SkyDir(),180)
        assert len(self.sources)==len(self.names), 'inconsistent list lengths'
        
    def lookup(self, name):
        """ return an ExtendedSource object, None if not found """
        for source in self.sources:
            if name==source.name:
                #source.free = source.model.free.copy()
                #if source.model.name=='LogParabola': source.free[-1]=False # E_break not free
                #return source
                # make a new object copied from original
                extsource= ExtendedSource(name=name, skydir=source.skydir,
                    model = source.model, 
                    spatial_model = source.spatial_model,
                    smodel= source.smodel,      # these reference copies needed
                    dmodel= source.spatial_model
                    )
                if extsource.model.name=='LogParabola': extsource.free[-1]=False # E_break not free
                return extsource    
        return None #raise Exception( ' extended source %s not found' % name)
  
def validate( ps):
    """ validate a Source: if not OK, give standard parameters, disable all but small flux level
    """
    model = ps.model
    if model.name=='LogParabola':
        norm, alpha, beta, eb = 10**model.p
        check = norm>1e-16 and alpha>0.5 and alpha<5 and beta<3
        if check: return
        print 'SkyModel warning for %-20s: out of range found %s' %(ps.name, model.p)
        model.p[:] = [-15, 0.4, -3, 3]
        ps.free[1:] = False
        model.cov_matrix[:] = 0 
    elif model.name=='ExpCutoff':
        pass
    else:
        print 'Skymodel warning: model name %s for source %s not recognized'%(model.name, ps.name)
    if np.any(np.diag(ps.model.cov_matrix)<0):
        print 'SkyModel warning for %-20s: invalid cov matrix ' %ps.name
        ps.model.cov_matrix[:] = 0 
  
