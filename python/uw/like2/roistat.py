"""
Manage likelihood calculations for an ROI

mostly class ROIstat, which computes the likelihood and its derivative from the lists of
sources (see .sourcelist) and bands (see .bandlike)

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roistat.py,v 1.12 2011/10/01 13:35:06 burnett Exp $
Author: T.Burnett <tburnett@uw.edu>
"""
import sys
import numpy as np
from . import bandlike
from . import sourcelist

         
class ROIstat(object):
    """ manage statistical analysis of an ROI
    
    Initialize from an existing ROIAnalysis object, for now. (Penalty in that
    have to repeat the convolutions -- but see its skip_setup)
    
    Contains two lists:
       * all sources, in self.sources
         order is, for now, the same as for ROIAnlysis
       * a list of BandLike objects, one per band, attribute all_bands
          each manages the computation of the likelihood for its band
    Notes:
        * This is not intended to be a replacement for ROIAnalysis: complete user-level
          functionality will be the function of a client 
        * the constructor takes an object containing a list of ROIBand objects and models, 
        * computation of the likelihood and its derivative, the basic task for 
          this class, is easily limited to a subset of the original set of bands,
          see select_bands
        * localization and extended source fits are handled by clients of this class
          by specifying a single source to be variable with the freelist parameter for
          initialize and setting reset=True in the update method 
    """
    
    def __init__(self, roi, bandsel=lambda b: True, quiet=False):
        """
        roi : ROIsetup object
            use to setup sources, find bands
        bandsel : function object returns bool
            returns True for a selected band
        """
        self.name = roi.name
        self.roi_dir = roi.roi_dir
        self.sources = sourcelist.SourceList(roi) 
        self.all_bands = bandlike.factory(filter(bandsel, roi.bands), self.sources)
        self.selected_bands = self.all_bands # the possible subset to analyze
        self.calls=0
        self.call_limit=1000
        self.quiet=quiet
        self.prior= ModelPrior(self.sources, self.sources.free)
    
    def __str__(self):
        return 'ROIstat for ROI %s: %d bands %d sources (%d free)' \
                % (self.name, len(self.all_bands), len(self.sources), sum(self.sources.free))
    @property
    def parameter_names(self):
        return self.sources.parameter_names
    def get_parameters(self):
        return self.sources.get_parameters()
    def set_parameters(self,parameters):
        # todo: return if no change
        if not np.allclose(parameters,self.get_parameters(),rtol=0,atol=1e-6):
            self.sources.set_parameters(parameters)
            self.update()
    def get_external(self):
        return 10**self.sources.get_parameters()
    def set_external(self, par):
        self.sources.set_parameters(np.log10(par))
    parameters = property(get_external, set_external, doc='array of free parameters')
    @property
    def uncertainties(self): return self.sources.uncertainties
    @property
    def energies(self):
        """ array of energies for selected bands (may include front and back for a given energy)"""
        return sorted(list(set([band.energy for band in self.selected_bands])))
        
    def initialize(self, freelist=None):
        """ reinitialize a set of sources, setting up angular distributions
        freelist : None or array of Bool
            if None, use list of sources with at least one free spectral parameter
        """
        if freelist is None: freelist = self.sources.free
        for i,band in enumerate(self.all_bands):
            #if not self.quiet: 
            #    status_string = '...initializing band %2d/%2d'%(i+1,len(self.all_bands))
            #    print status_string;sys.stdout.flush()
            band.initialize(freelist)
        self.prior.initialize(freelist)

    def select_bands(self, bandsel=lambda b: True, indices=None):
        """ select a subset of the bands for analysis
        bandsel : function of a ROIBand that returns bool, like lambda b: b.e<200
        To restore, call with no arg
        Note that one could also replace selected_bands with a subset of all_bands
        
        """
        if indices is not None:
            self.selected_bands = self.all_bands[indices]
        else:
            self.selected_bands = np.array([bs for bs in self.all_bands if bandsel(bs.band)])
        if not self.quiet:print 'selected subset of %d bands for likelihood analysis' % len(self.selected_bands)
        
    def update(self, reset=False):
        """ perform update on all selected bands, and variable sources
        if reset is True, assume that must also reinitialize angular dependence of sources
        """
        map(lambda s: s.update(reset), self.selected_bands)
        self.prior.update()
        
    def log_like(self):
        """ return sum of log likelihood for all bands
        """
        return sum([blike.log_like() for blike in self.selected_bands]) + self.prior.log_like()
        
    def __call__(self, par):
        """ (negative) log likelihood as a function of the free parameters par 
        appropriate for minimizing
        """
        self.set_parameters(par)
        self.update()
        self.calls +=1
        if self.calls>self.call_limit:
            raise RuntimeError('function call limit, %d exceeded' %self.calls)
        return -self.log_like()

    def gradient(self, parameters=None):
        """ gradient of -log(like), or the call interface, with respect to parameters
            (note that the individual gradients assume -log(like))
        """
        if parameters is not None: 
            self.set_parameters(parameters)
        self.update()
        t = np.array([blike.gradient() for blike in self.selected_bands]).sum(axis=0)\
                + self.prior.gradient()
        par = self.get_parameters()
        assert len(t)==len(par), 'inconsistent number of free parameters'
        # this is required by the convention in all of the Models classes to use log10 for external
        jacobian= 10**par/np.log10(np.e)
        return t*jacobian
        
    def chisq(self):
        return sum([blike.chisq() for blike in self.selected_bands])
        
    def dump(self, **kwargs):
        map(lambda bs: bs.dump(**kwargs), self.selected_bands)
        
class NoPrior(object):
    def __call__(self,x): return 0
    def grad(self,x): return 0
    def delta(self,x): return 0
    def apply_limit(self,x): return x
    def __str__(self): return 'NoPrior'
    
class LimitPrior(NoPrior):
    def __init__(self,a,b,c=1e2):
        self.a,self.b,self.c=a,b,c
    def __str__(self):
        return 'LimitPrior(%f,%f)' % (self.a,self.b)
    def apply_limit(self,x):
        if   x<self.a: return self.a
        elif x>self.b: return self.b
        return x
    def __call__(self,x):
        if x>self.b:   return self.c*(x-self.b)**2
        elif x<self.a: return self.c*(x-self.a)**2
        return 0
    def grad(self,x):
        if   x>self.b: return 2*self.c*(x-self.b)
        elif x<self.a: return 2*self.c*(x-self.a)
        return 0
    def delta(self,x):
        if   x>self.b: return x-self.b
        elif x<self.a: return x-self.a
        return 0
    
        
class ModelPrior(object):
    """ manage addition of priors to the likelihood"""

    def __init__(self, sources, free=None):   
        self.sources = sources
        self.models = np.array([s.spectral_model for s in sources])
        self.initialize(free)
        self.update()
        self.enabled=True

    def initialize(self, free):
        self.free = free if free is not None else self.sources.free
        self.free_models = self.models[self.free]
        self.npar = np.sum([np.sum(m.free) for m in self.free_models])
        self.priors= []
        for model in self.free_models:
            for pname in np.array(model.param_names)[model.free]:
                if pname=='Index':
                   prior = LimitPrior(np.log10(1e-3), np.log10(3.0), 10.)
                elif pname=='Norm':
                    if model.name =='PowerLaw': #for ring
                        prior = LimitPrior(np.log10(0.5), np.log10(1.5), 10.)
                    else: prior = LimitPrior(-16, -8, 1e4)
                elif pname=='Scale': # used by iso
                   prior = LimitPrior(np.log10(1e-4), np.log10(1.5), 10.)
                else:
                    prior = NoPrior()
                self.priors.append(prior)
  
    def update(self, reset=False):
        self.pars = self.sources.get_parameters()

    def limit_pars(self, update=False):
        self.update()
        for i, prior in enumerate(self.priors):
            self.pars[i] = prior.apply_limit(self.pars[i])
        if update:
            self.sources.set_parameters(self.pars)
        
    def log_like(self):
        if not self.enabled: return 0
        t = np.sum(prior(x) for prior,x in zip(self.priors, self.pars))
        if False: # t>10.0: 
            n = np.arange(len(self.pars))[np.abs(self.gradient())>0]
            print 'likelihood penalty: %.1f %s' % (t,n)
        return -t
    
    def check(self):
        names = np.array(self.sources.parameter_names)
        deltas = np.array([prior.delta(x) for prior,x in zip(self.priors, self.pars)])
        sel= deltas!=0
        return np.vstack([names[sel], self.pars[sel], deltas[sel]])
        
    def gradient(self):
        if not self.enabled: return np.zeros(len(self.pars))
        return np.array([prior.grad(x) for prior,x in zip(self.priors, self.pars)])