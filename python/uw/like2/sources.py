"""
Source descriptions for SkyModel
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sources.py,v 1.26 2013/10/11 16:10:39 burnett Exp $

"""
import os, pickle, glob, types, copy
import numpy as np
from skymaps import SkyDir, Band, IsotropicSpectrum
import skymaps 
from uw.like import  Models
from uw.like import pointspec_helpers
from uw.utilities import parmap 

# convenience adapters 
def LogParabola(*pars, **kw):return Models.LogParabola(p=pars, **kw)
def PowerLaw(*pars, **kw):   return Models.PowerLaw(p=pars, **kw)
def ExpCutoff(*pars, **kw):  return Models.ExpCutoff(p=pars, **kw)
def PLSuperExpCutoff(*pars, **kw): return Models.PLSuperExpCutoff(p=pars, **kw)
    
def convert_model(oldmodel):
    """ convert the original version to the new one with parameter mappers
    """
    if hasattr(oldmodel, 'mappers'): return oldmodel #new, or already converted
    pars = 10**oldmodel._p
    # absorb index_offset into index (but avoid exactly zero for now)
    if oldmodel.name=='PowerLaw':
        pars[1] -= oldmodel.index_offset
        if pars[1]==0: pars[1]=1e-3
    elif oldmodel.name=='FrontBackConstant':
        # different constructor for this guy
        # assume still use log
        return Models.FrontBackConstant(pars[0],pars[1])
    model = eval('Models.%s()' % oldmodel.name)
    model.set_all_parameters(pars)
    model.free = oldmodel.free
    # try to set complete state: e0 tricky
    if hasattr(oldmodel, 'e0') and 'e0' in model.default_extra_params : 
        model['e0'] = oldmodel.e0 
    # jacobian for conversion from all log10 to default
    j =  (np.log(10)*pars) / model.dexternaldinternal()
    model.internal_cov_matrix = j * oldmodel.cov_matrix * j.T
    #model.internal_cov_matrix = oldmodel.cov_matrix
    return model
    

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
        if self.model.name=='PowerLaw':
            par,sig = self.model.statistical()
            self.model = LogParabola(*(list(par)+[1e-3, self.model.e0]))
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
                self.model=LogParabola(f*(gamma-1)/emin, gamma, 1e-3, emin)
            except Exception, msg:
                print 'Failed to create LogParabola for source %s, pars= %s'% (self.name, (f,gamma,emin))
                raise
            self.model.free[2:]=False
        elif self.model.name=='LogParabola':
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
            
            
        self.free = self.model.free.copy()  # save copy of initial free array to restore
    def freeze(self, freeze):
        self.model.free[:] = False if freeze else self.free
    def __str__(self):
        return self.name + ' '+ self.skydir.__str__() +' '+ self.model.name \
                +  (' (free)' if np.any(self.model.free) else ' (fixed)')
    def __repr__(self):
        return '%s.%s: %s' % (self.__module__,self.__class__.__name__ , self.name)

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
        
class GlobalSource(Source):
    def __init__(self, **kwargs):
        super(GlobalSource, self).__init__(**kwargs)
        assert self.skydir is None # used as a flag
    
class GlobalSourceList(list):
    """ a list, indexed by ROI number, of GLobalSource lists
        each element is a list if the GlobalSource objects, includeing Models, in the ROI
    """
    def __repr__(self):
        return '%s.%s: %d elements' % (self.__module__, self.__class__.__name__, len(self))

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
 

class ExtendedCatalog( pointspec_helpers.ExtendedSourceCatalog):
    """ subclass to add this lookup function """

    def __init__(self, *pars, **kwargs):
        """ initialize by also filling an array with all source spectral models"""
        self.alias = kwargs.pop('alias', dict())
        super(ExtendedCatalog,self).__init__(*pars, **kwargs)
        # create list of sources using superclass, for lookup by name
        self.sources = [self.get_sources(self.dirs[i], 0.1)[0] for i in range(len(self.names))]
        for source in self.sources:
            model = source.model
            if model.mappers[0].__class__.__name__== 'LimitMapper':
                #print 'converting mappers for model for source %s, model %s' % (source.name, model.name)
                source.model = eval('Models.%s(p=%s)' % (model.name, list(model.get_all_parameters())))


    def realname(self, cname):
        """ cname was truncated"""
        if cname in self.names: return cname
        for name in self.names:
            assert name is not None, 'bad name'
            t = name.replace(' ','')
            if t==cname: return name
        assert 'compressed name %s not found in list of names, %s' %(cname,self.names)

    def lookup(self, name):
        """ return an ExtendedSource object by name, None if not found """
        aname = self.alias.get(name,name) #alias will be the new name
        try:
            i = list(self.names).index(name)
        except ValueError:
            return None
        source = self.sources[i]
        # make a new object copied from original
        if source.model.name=='BrokenPowerLaw': #convert this
            model = Models.LogParabola()
        else: model = source.model
        if model.name=='LogParabola': model.free[-1]=False # E_break ne free
        ### seems to be necessary for some models created from 
        if model.mappers[0].__class__.__name__== 'LimitMapper':
            print 'wrong mappers: converting model for source %s, model %s' % (name, model.name)
            model = eval('Models.%s(p=%s)' % (model.name, list(model.get_all_parameters())))
        extsource= ExtendedSource(name=self.realname(aname), skydir=source.skydir,
            model = model, 
            spatial_model = source.spatial_model,
            smodel= model,      # these reference copies needed
            dmodel= source.spatial_model
            )
        if extsource.model.name=='LogParabola': extsource.free[-1]=False # E_break not free
        return extsource  

  
def validate( ps, nside, filter):
    """ validate a Source: if not OK, reset to standard parameters, disable all but small flux level
    """
    if filter is not None:
        ret = filter(ps)
        if not ret: 
            print 'SkyModel: removed source %s' % ps.name
            return ret
    model = ps.model
    assert hasattr(model, 'npar'), 'no npar: %s ' % model.__dict__
    if '_p' not in model.__dict__:
        model.__dict__['_p'] = model.__dict__.pop('p')  # if loaded from old representation
    hpindex = lambda x: Band(nside).index(x)
    if model.name=='LogParabola':
        norm, alpha, beta, eb = model.get_all_parameters() #10**model.p
        #if norm<1e-18: model[0]=1e-18 #quietly prevent too small
        if beta<0.01: # linear
            check =  norm< 1e-4 and alpha>-0.5 and alpha<=5 
            if check: return True
            print 'SkyModel warning for %-20s(%d): out of range, norm,alpha=%.2e %.2f' %(ps.name, hpindex(ps.skydir),norm,alpha)
            #assert False, 'debug'
            #model[:]= [1e-15, 2.4, 1e-3, 1000]
            #ps.free[1:] = False
            #model.cov_matrix[:] = 0 
        else: #log parabola
            check =  alpha>=-0.5 and alpha<100 and beta<100
            if check: return True
            print 'SkyModel warning for %-20s(%d): out of range, norm,alpha,beta=%.2e %.2f %.2f %.2f'\
                    %(ps.name, hpindex(ps.skydir),norm,alpha,beta, eb)
            #ps.free[1:] = False
            #model.cov_matrix[:] = 0 
        
    elif model.name=='PLSuperExpCutoff':
        norm, gamma, ec, b = model.get_all_parameters() #10**model.p
        #if np.any(np.diag(model.cov_matrix)<0): model.cov_matrix[:]=0 
        if norm<1e-18: model[0]=1e-18 #quietly prevent too small
        check =  gamma>=-0.5 and gamma<5 and ec>100
        if check: return True
        print 'SkyModel warning for %-20s(%d): out of range, norm, gamma,ec %s' %(ps.name, hpindex(ps.skydir),model.get_all_parameters())
        #model[:] = [1e-15, 2.2, 2000.]
        #model.cov_matrix[:] = 0 
    else:
        print 'Skymodel warning: model name %s for source %s not recognized'%(model.name, ps.name)
    if np.any(np.diag(ps.model.internal_cov_matrix)<0):
        print 'SkyModel warning for %-20s: invalid cov matrix ' %ps.name
        #ps.model.cov_matrix[:] = 0 
    return True
  
