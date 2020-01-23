"""
Manage the diffuse sources

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/Attic/diffusedict.py,v 1.4 2013/10/14 15:11:42 burnett Exp $

author:  Toby Burnett
"""
import os, types, json, collections
#import  pickle, glob, copy, zipfile
#import numpy as np
#import pandas as pd
#from .pub import healpix_map make this local in future
#from uw.utilities import keyword_options, convolution

import skymaps #from Science Tools: for DiffuseFunction, IsotropicSpectrum
from uw.like import  Models
def PowerLaw(*pars, **kw):   return Models.PowerLaw(p=pars, **kw)

class DiffuseBase(object):
    """Base class for diffuse sources
    expect subclasses to implement SkySpectrum interface
    """
    def load(self): pass
    def __repr__(self):
        t= self.__class__.__name__ + ': '+self.filename
        if hasattr(self, 'kw') and self.kw is not None:
            t+= '\n'+ '\n'.join(['\t    %-10s: %s' % item for item in self.kw.items() if item[0]!='filename'])
        return t
    @property
    def name(self):
        return self.__class__.__name__
    

class Isotropic(DiffuseBase, skymaps.IsotropicSpectrum):
    """Implement the isotropic diffuse, by wrapping skymaps.IsotopicSpectrum
    """
    def __init__(self, filename):
        self.filename=filename
        super(Isotropic, self).__init__(filename)
    

class MapCube(DiffuseBase, skymaps.DiffuseFunction):
    """ wrapper for eventual invocation of skymaps.DiffuseFunction, which interprets a map cube
    load() must be called before use
    """
    def __init__(self, filename):
        """ filename: string or dict
                if dict, expect 'filename' key
        """
        self.filename = filename
        if not os.path.exists(self.filename):
            self.filename = os.path.expandvars(os.path.join('$FERMI','diffuse',self.filename))
        assert os.path.exists(self.filename), 'DiffuseFunction file "%s" not found' % self.filename
        self.loaded =  False
            
    def load(self, interpolate=False):
        if  self.loaded: return
        self.loaded=True
        if not interpolate: 
            print ('loading diffuse file %s: warning, not interpolating' %self.filename)
        super(MapCube,self).__init__(self.filename, 1000., interpolate)

  
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
        except Exception as msg:
            print ('Failure to evaluate IsotropicSpectralFunction %s : %s' % (self.expression, msg))
            raise
        self.energy=1000
    def __repr__(self):
        return '%s: %s' % (self.__class__.__name__, self.expression )
    def __call__(self, skydir, energy=None):
        return self.spectral_function(self.energy if energy is None else energy)
    def setEnergy(self, energy): self.energy=energy


class CachedMapCube(DiffuseBase):
    def __init__(self, zipfilename):
        self.filename = zipfilename
    def __call__(self, skydir, energy=None):
        return np.nan

class DiffuseList(list):
    """contain the subtype, or front/back list of DiffuseBase objects If only one, applied to all
    """
    def __init__(self, inlist):
        super(DiffuseList,self).__init__(inlist)
    def __repr__(self):
        return  '%s.%s' % (self.module, self.__class__.__name__)
        
    def __getitem__(self, index):
        if len(self)==1: index=0
        return super(DiffuseList,self).__getitem__(index) 
    def load(self):
        for x in self:
            x.load()
        
    
def diffuse_factory(value):
    """
    Create a DiffuseList object 
    """
    isdict = issubclass(value.__class__, dict)
    if  not hasattr(value, '__iter__') or isdict:
        value  = (value,)
    
    if isdict:
        files = DiffuseList([val['filename'] for val in value])
        kws = value
    else:
        files = DiffuseList(value)
        kws = None

    # checking only first element, and only upto a comma, if not with '('
    f = files[0].split(',')[0] if files[0].find('(')<0 else files[0]
    ext = os.path.splitext(f)[-1]
    try:
        dfun = {'.txt': Isotropic, 
            '.fit': MapCube, '.fits': MapCube,
            '.zip': CachedMapCube,
            ')': IsotropicSpectralFunction, 
            }[ext if ext[-1]!=')' else ')']
    except Exception as msg:
        raise Exception('File type, "%s", for diffuse not recognized, from "%s":%s'% (ext, files, ext))
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
        value: (both,) or (front,back)  diffuse objects determined by the extension:
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
        except Exception as msg:
            print ('Failed to open model at %s: %s' % (modeldir, msg))
            raise
        super(DiffuseDict, self).__init__()
        for key, value in self.spec.items():
            print (key,value)
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
            s.smodel = s.model
            global_sources.append(s)
        return global_sources

def main():
    pass
    
if __name__=='__main__':
    main()