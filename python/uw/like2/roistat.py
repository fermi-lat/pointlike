"""
Manage likelihood calculations for an ROI
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roistat.py,v 1.2 2011/08/18 16:46:41 burnett Exp $
Author: T.Burnett <tburnett@uw.edu>
"""

import numpy as np
from uw.utilities import fitter
from . import bandmodel
from . import sourcelist

         
class ROIstat(object):
    """ manage statistical analysis of an ROI
    Initialize from an existing ROIAnalysis object, for now. (Penalty in that
    have to repeat the convolutions)
    
    Contains two lists:
       * all sources, in self.sources
         order is, for now, the same as for ROIAnlysis
       * a list of BandModelStat objects, one per band, attribute all_bands
          each manages the computation of the likelihood for its band
    The constructor takes an existing ROIAnalysis object, 
    from which it extracts the sources and bands.
    the fit() method is equiavlent to the same method in ROIanalysis
    
    Not implemented yet: 
        * fits for SED plots: that is, fits combining front and back for each band
        * localization - probably just use or adapt code in ROiAnalysis 
        * computation of TS -- same
        * modifying list of sources in place
        * defining a subset of variables?
        ...
      
    """
    
    def __init__(self, roi, bandsel=lambda b: True):
        """
        roi : ROIanalysis object
            use to setup sources, find bands
        bandsel : function object returns bool
            returns True for a selected band
        """
        self.name = roi.name
        self.sources = sourcelist.SourceList(roi) 
        self.all_bands = bandmodel.factory(filter(bandsel, roi.bands), self.sources)
        self.selected_bands = self.all_bands # the possible subset to analyze
        self.calls=0
    
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

    def initialize(self):
        """ used only if free changes, to reinitialize all free models """
        map(lambda s: s.initialize(self.sources.free), self.all_bands )

    def select_bands(self, bandsel=lambda b: True):
        """ select a subset of the bands for analysis
        bandsel : function of a ROIBand that returns bool, like lambda b: b.e<200
        To restore, call with no arg
        Note that one could just replace selected_bands with a subset of all_bands
        """
        self.selected_bands = np.array([bs for bs in self.all_bands if bandsel(bs.band)])
        print 'selected subset of %d bands for likelihood analysis' % len(self.selected_bands)
        
    def update(self):
        map(lambda s: s.update(), self.selected_bands)
        
    def log_like(self):
        """ return sum of log likelihood for all bands
            (could easily set for a subset of the bands --TODO)
        """
        return sum([bstat.log_like() for bstat in self.selected_bands])
        
    def __call__(self, par):
        """ (negative) log likelihood as a function of the free (internal rep) parameters par """
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
        t = np.array([bstat.gradient() for bstat in self.selected_bands]).sum(axis=0)
        # this is required by the convention in all of the Models classes to use log10 for external
        jacobian= 10**self.get_parameters()/np.log10(np.e)
        return t*jacobian
        
    def chisq(self):
        return sum([bstat.chisq() for bstat in self.selected_bands])
        
    def fit(self, **kwargs):
        """ Perform fit, return fitter object to examine errors, or refit
        """
        fit_kw = dict(use_gradient=True, estimate_errors=True)
        fit_kw.update(kwargs)
        initial_value, self.calls = self.log_like(), 0
        mm = fitter.Minimizer(self)
        mm(**fit_kw)
        print '%d calls, likelihood improvement: %.1f' % (self.calls, self.log_like() - initial_value)
        if fit_kw['estimate_errors']:
            self.sources.set_covariance_matrix(mm.cov_matrix)
        return mm

    def dump(self, **kwargs):
        map(lambda bs: bs.dump(**kwargs), self.selected_bands)
        
##### these are for testing, some may turn into methods
def test(roi):
    r= ROIstat(roi, lambda b: b.e<200)
    print r
    return r
    
def source_scan(s, i=0, dom=np.linspace(0.5,1.5, 11)):
    """ scan likelihood as a function of the relative normalization of model component"""
    m = s.sources[i].model
    a = s.log_like()
    print s.sources[i].name, m
    print 'nominal:',a
    print 'factor  logl'
    norm = m[0]
    for r in dom:
        m[0]=norm*r
        s.update()
        print '%5.2f %9.2f ' %(r,  s.log_like()-a)
    m[0] = norm
    s.update()

def par_scan(s, i, q=0.1,  dom=None): #np.linspace(-0.05, 0.05, 11)):
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
    delta('  galactic',     b.bg_counts[0], bscounts[0])
    delta('  isotropic',     b.bg_counts[1], bscounts[1])
    if numbg>2:
        delta('  extended',     b.bg_counts[2], bscounts[2])
    delta('point sources', b.ps_all_counts, bscounts[numbg:].sum())
    delta( 'model sum',  b.bg_all_counts + b.ps_all_counts, bs.counts)
    delta( 'data sum',  sum(b.pix_counts), sum(bs.data))
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
    """ adapt an ROI to my fitter interface """
    def __init__(self, roi):  self.roi = roi
    def __call__(self, par):  return self.roi.logLikelihood(par)
    def get_parameters(self): return self.roi.get_parameters()
    def set_parameters(self,par): self.roi.set_parameters(par)
    def gradient(self, par):   return self.roi.gradient(par)
    def fit(self, **kwargs):
        mm = fitter.Minimizer(self)
        return mm(**kwargs)
        
