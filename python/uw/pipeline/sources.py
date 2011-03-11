"""
Source descriptions for SkyModel
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/sources.py,v 1.6 2011/03/07 00:07:44 burnett Exp $

"""
import os, pickle, glob, types, copy
import numpy as np
from skymaps import SkyDir, Band, IsotropicSpectrum, DiffuseFunction
from uw.like import  Models
from uw.like import pointspec_helpers

# convenience adapters 
def LogParabola(*pars):return Models.LogParabola(p=pars)
def PowerLaw(*pars):   return Models.PowerLaw(p=pars)
def ExpCutoff(*pars):  return Models.ExpCutoff(p=pars)
    
class Source(object):
    """ base class for various sources
    """
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        assert self.name is not None, 'bad source name'
        if self.skydir is None:
            # global source: keep original model
            self.free = self.model.free.copy()  # save copy of initial free array to restore
            return
        if 'model' not in kwargs or self.model is None:
            self.model=LogParabola(1e-14, 2.2, 1e-3, 1e3)
            self.model.free[2:]=False
        if type(self.model)==types.StringType:
            self.model = eval(self.model)
        if self.model.name=='PowerLaw':
            par,sig = self.model.statistical()
            self.model = LogParabola(list(par)+[1e-3, self.model.e0])
            self.model.free[2:]=False
        elif self.model.name=='PLSuperExpCutoff':
            par,sig=self.model.statistical()
            self.model = ExpCutoff(par[:-1])
        elif self.model.name=='LogParabola':
            if 'index_offset' not in self.model.__dict__:
                self.model.index_offset = 0
            self.model.free[-1]=False
            if self.model.cov_matrix[3,3]<0:  
                self.model.cov_matrix[3,3]= 100.
                print 'fix covariance matrix for source %s' % self.name
        elif self.model.name=='PowerLawFlux':
            f, gamma = self.model.get_all_parameters() #10**self.model.p
            emin = self.model.emin
            self.model=LogParabola(f*(gamma-1)/emin, gamma, 1e-3, emin)
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
    def __init__(self, **kwargs):
        kwargs.update(spatial_model=None) # allow test for extent (no extent!)
        super(PointSource, self).__init__(**kwargs)
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
        self.alias = kwargs.pop('alias', dict())
        super(ExtendedCatalog,self).__init__(*pars, **kwargs)
        self.sources = self.get_sources(SkyDir(),180)
        assert len(self.sources)==len(self.names), 'inconsistent list lengths'

    def realname(self, cname):
        """ cname was truncated"""
        if cname in self.names: return cname
        for name in self.names:
            assert name is not None, 'bad name'
            t = name.replace(' ','')
            if t==cname: return name
        assert 'compressed name %s not found in list of names, %s' %(cname,self.names)

    def lookup(self, name):
        """ return an ExtendedSource object, None if not found """
        aname = self.alias.get(name,name) #alias will be the new name
        #if aname != name: print 'Renaming extended source %s to %s' % (name, aname)
        for source in self.sources:
            if source.name in (name, aname, aname.replace(' ',''), name.replace(' ','')):
                #source.free = source.model.free.copy()
                #if source.model.name=='LogParabola': source.free[-1]=False # E_break not free
                #return source
                # make a new object copied from original
                if source.model.name=='BrokenPowerLaw': #convert this
                    model = Models.LogParabola()
                else: model = source.model
                extsource= ExtendedSource(name=self.realname(aname), skydir=source.skydir,
                    model = model, 
                    spatial_model = source.spatial_model,
                    smodel= source.smodel,      # these reference copies needed
                    dmodel= source.spatial_model
                    )
                if extsource.model.name=='LogParabola': extsource.free[-1]=False # E_break not free
                return extsource    
        return None #raise Exception( ' extended source %s not found' % name)
  
def validate( ps, nside, filter):
    """ validate a Source: if not OK, reset to standard parameters, disable all but small flux level
    """
    if filter is not None:
        ret = filter(ps)
        if not ret: 
            print 'SkyModel: removed source %s' % ps.name
            return ret
    model = ps.model
    if '_p' not in model.__dict__:
        model.__dict__['_p'] = model.__dict__.pop('p')  # if loaded from old representation
    hpindex = lambda x: Band(nside).index(x)
    if model.name=='LogParabola':
        norm, alpha, beta, eb = model.get_all_parameters() #10**model.p
        if beta<0.01: # linear
            check = norm>1e-18 and norm< 1e-4 and alpha>0.25 and alpha<5 
            if check: return True
            print 'SkyModel warning for %-20s(%d): out of range, norm,alpha=%.2e %.2f' %(ps.name, hpindex(ps.skydir),norm,alpha)
            model.set_all_parameters( [-15, 0.4, -3, 3], internal=True)
            ps.free[1:] = False
            model.cov_matrix[:] = 0 
        else: #log parabola
            check = norm>1e-18 and alpha>1e-4 and alpha<10 and beta<10
            if check: return True
            print 'SkyModel warning for %-20s(%d): out of range, norm,alpha=%.2e %.2f' %(ps.name, hpindex(ps.skydir),norm,alpha)
            model.set_all_parameters(  [-15, 0.4, -3, 3], internal=True)
            ps.free[2:] = False
            model.cov_matrix[:] = 0 
        
    elif model.name=='ExpCutoff':
        norm, gamma, ec = model.get_all_parameters() #10**model.p
        if np.any(np.diag(model.cov_matrix)<0): model.cov_matrix[:]=0 
        check = norm>1e-18 and gamma>1e-5 and gamma<5 and ec>100
        if check: return True
        print 'SkyModel warning for %-20s(%d): out of range, ressetting from %s' %(ps.name, hpindex(ps.skydir),model.p)
        model.set_all_parameters( [-11, 0, 3], internal=True)
        model.cov_matrix[:] = 0 
    else:
        print 'Skymodel warning: model name %s for source %s not recognized'%(model.name, ps.name)
    if np.any(np.diag(ps.model.cov_matrix)<0):
        print 'SkyModel warning for %-20s: invalid cov matrix ' %ps.name
        ps.model.cov_matrix[:] = 0 
    return True
  
