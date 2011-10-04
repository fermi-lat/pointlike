"""
Source descriptions for SkyModel
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/sources.py,v 1.1 2011/09/28 16:56:26 burnett Exp $

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
            self.model = ExpCutoff(*par[:-1])
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
    def copy(self):
        """ return a new PointSource object, with a copy of the model, others"""
        ret = PointSource(**self.__dict__)
        ret.model = self.model.copy()
        return ret
        
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
    def copy(self):
        """ return a new ExtendSource object, with a copy of the model object"""
        ret = ExtendedSource(**self.__dict__)
        ret.model = self.model.copy()
        return ret
    
class DiffuseDict(dict):
    """ create a dictionary of global diffuse objects
        key:   a string defined by the filename preceding an underscore
        value: (both) or (front,back)  diffuse objects determined by the extension:
            txt: IsotropicSpectrum
            fit or fits: DiffuseFunction
    
    """
    def __init__(self, diffuse):
        """ diffuse: a list, where each entry is a file name or a tuple of one or two file names, for front and back
        """
        assert len(diffuse)<4, 'expect 2 or 3 diffuse names, or front/back tuples'
        # convert each single entry to a tuple: assume iterables are tuples of strings
        tuplelist = map( lambda x: (x,) if not hasattr(x,'__iter__') else x, diffuse)
        keys = map( lambda x: x[0].split('_')[0], tuplelist) # key or name from first one
        for key, files in zip(keys, tuplelist):
            full_files = map( lambda f: os.path.expandvars(os.path.join('$FERMI','diffuse',f)), files)
            check = map(lambda f: os.path.exists(f), full_files) 
            assert all(check), 'not all diffuse files %s found' % full_files
            ext = os.path.splitext(full_files[0])[-1]
            try:
                dfun = {'.txt':IsotropicSpectrum, '.fit':DiffuseFunction, '.fits':DiffuseFunction}[ext]
            except:
                raise Exception('File type, %s, for diffuse not recognized'% ext)
            self[key]= map(dfun, full_files) 


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
        #if norm<1e-18: model[0]=1e-18 #quietly prevent too small
        if beta<0.01: # linear
            check =  norm< 1e-4 and alpha>0.25 and alpha<5 
            if check: return True
            print 'SkyModel warning for %-20s(%d): out of range, norm,alpha=%.2e %.2f' %(ps.name, hpindex(ps.skydir),norm,alpha)
            model.set_all_parameters( [-15, 0.4, -3, 3], internal=True)
            ps.free[1:] = False
            model.cov_matrix[:] = 0 
        else: #log parabola
            check =  alpha>1e-4 and alpha<100 and beta<100
            if check: return True
            print 'SkyModel warning for %-20s(%d): out of range, norm,alpha=%.2e %.2f %.2f'\
                    %(ps.name, hpindex(ps.skydir),norm,alpha,beta)
            model[:]= [1e-15, 5.0, 10.0, 10000]
            ps.free[2:] = False
            model.cov_matrix[:] = 0 
        
    elif model.name=='ExpCutoff':
        norm, gamma, ec = model.get_all_parameters() #10**model.p
        if np.any(np.diag(model.cov_matrix)<0): model.cov_matrix[:]=0 
        if norm<1e-18: model[0]=1e-18 #quietly prevent too small
        check =  gamma>1e-10 and gamma<5 and ec>100
        if check: return True
        print 'SkyModel warning for %-20s(%d): out of range, ressetting from %s' %(ps.name, hpindex(ps.skydir),model.get_all_parameters())
        model[:] = [1e-15, 1.0, 500.]
        model.cov_matrix[:] = 0 
    else:
        print 'Skymodel warning: model name %s for source %s not recognized'%(model.name, ps.name)
    if np.any(np.diag(ps.model.cov_matrix)<0):
        print 'SkyModel warning for %-20s: invalid cov matrix ' %ps.name
        ps.model.cov_matrix[:] = 0 
    return True
  
