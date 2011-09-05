"""
Manage likelihood calculations for an ROI

mostly class ROIstat, which computes the likelihood and its derivative from the lists of
sources (see .sourcelist) and bands (see .bandlike)

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roistat.py,v 1.7 2011/09/02 16:30:43 burnett Exp $
Author: T.Burnett <tburnett@uw.edu>
"""

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
        * the constructor takes, for now, an existing ROIAnalysis object, 
          from which it extracts the sources and bands.
        * computation of the likelihood and its derivative, the basic task for 
          this class, is easily limited to a subset of the original set of bands,
          see select_bands
        * the fit() method is equivalent to the same method in ROIanalysis,
          except that it allows specification of a subset of parameters to use
          for optimization
        * localization and extended source fits are handled by clients of this class
          by specifying a single source to be variable with the freelist parameter for
          initialize and setting reset=True in the update method 
        
    Functionality planned, but yet implemented: 
        * fits for SED plots: that is, fits combining front and back for each band
        * modifying list of sources in place
        ...
        Note that these should be implemented by client classes, and
        not added here.
    """
    
    def __init__(self, roi, bandsel=lambda b: True, quiet=False):
        """
        roi : ROIanalysis object
            use to setup sources, find bands
        bandsel : function object returns bool
            returns True for a selected band
        """
        self.name = roi.name
        self.sources = sourcelist.SourceList(roi) 
        self.all_bands = bandlike.factory(filter(bandsel, roi.bands), self.sources)
        self.selected_bands = self.all_bands # the possible subset to analyze
        self.calls=0
        self.quiet=quiet
    
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
    def energies(self):
        """ array of energies for selected bands (may include front and back for a given energy)"""
        return np.array(set([band.energy for band in self.selected_bands]))
        
    def initialize(self, freelist=None):
        """ reinitialize a set of sources, setting up angular distributions
        freelist : None or array of Bool
            if None, use list of sources with at least one free spectral parameter
        """
        if freelist is None: freelist = self.sources.free
        map(lambda s: s.initialize(freelist), self.all_bands )

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
        
    def log_like(self):
        """ return sum of log likelihood for all bands
        """
        return sum([blike.log_like() for blike in self.selected_bands])
        
    def __call__(self, par):
        """ (negative) log likelihood as a function of the free parameters par 
        appropriate for minimizing
        """
        self.set_parameters(par)
        self.update()
        self.calls +=1
        return -self.log_like()

    def gradient(self, parameters=None):
        """ gradient of -log(like), or the call interface, with respect to parameters
            (note that the individual gradients assume -log(like))
        """
        if parameters is not None: 
            self.set_parameters(parameters)
        self.update()
        t = np.array([blike.gradient() for blike in self.selected_bands]).sum(axis=0)
        # this is required by the convention in all of the Models classes to use log10 for external
        jacobian= 10**self.get_parameters()/np.log10(np.e)
        return t*jacobian
        
    def chisq(self):
        return sum([blike.chisq() for blike in self.selected_bands])
        
    def dump(self, **kwargs):
        map(lambda bs: bs.dump(**kwargs), self.selected_bands)
        
##### these are for testing, some may turn into methods
    
def par_scan(s, i, q=0.1,  dom=None): 
    """ scam parameter i, showing relative likelihood and gradient
    """
    if dom is None:
        dom = np.linspace(-q,q, 11)
    par = s.get_parameters().copy()
    d = np.zeros(len(par))
    d[i]=1.0
    for x in dom:  
        print '%6.3f %8.4f %8.4f' % (x,s(par+x*d)-s(par),s.gradient(par+x*d)[i])
    s.set_parameters(par)
    
def compare(s, i=0):
    def delta(name, a, b):
        print '%-15s %10.1f%10.1f%10.2f' % (name, a, b, a-b)
    bs = s.all_bands[i]
    #bs.dump()
    b = bs.bandmodels[i].band
    print 'quantity           pointlike  roistat  difference'
    delta('log likelihood', -b.logLikelihood(), bs.log_like())
    delta('pixel like',\
        (b.pix_counts * np.log(b.bg_all_pix_counts + b.ps_all_pix_counts)).sum(),\
         sum( bs.data * np.log(bs.model) ))
    bscounts = np.array([bm.counts for bm in bs.bandmodels])
    numbg = len(b.bg_counts)
    delta('diffuse sources', b.bg_all_counts, bscounts[:numbg].sum())
    delta('  galactic',      b.bg_counts[0], bscounts[0])
    delta('  isotropic',     b.bg_counts[1], bscounts[1])
    if numbg>2:
        delta('  extended',  b.bg_counts[2], bscounts[2])
    delta('point sources',   b.ps_all_counts, bscounts[numbg:].sum())
    delta( 'model sum',      b.bg_all_counts + b.ps_all_counts, bs.counts)
    delta( 'data sum',       sum(b.pix_counts), sum(bs.data))
    return b, bs
    
def pardump(self, eps=1e-1):
    print 'variable parameters'
    par = self.get_parameters()
    grad = self.gradient(par)
    f0 = self(par)
    for i, name in enumerate(self.parameter_names):
        dpar = par.copy()
        dpar[i] += eps/grad[i]
        print '%-20s%10.2e%10.2e %10.2e'% (name, par[i], grad[i], (self(dpar)-f0)*grad[i]/eps)
        self.set_parameters(par)

class ROIfit(object):
    """ adapt an ROI to the fitter interface """
    def __init__(self, roi):  self.roi = roi
    def __call__(self, par):  return self.roi.logLikelihood(par)
    def get_parameters(self): return self.roi.get_parameters()
    def set_parameters(self,par): self.roi.set_parameters(par)
    def gradient(self, par):   return self.roi.gradient(par)
    def fit(self, **kwargs):
        mm = fitter.Minimizer(self)
        return mm(**kwargs)
        
