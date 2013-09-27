"""
Manage likelihood calculations for an ROI

mostly class ROIstat, which computes the likelihood and its derivative from the lists of
sources (see .sourcelist) and bands (see .bandlike)

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/roistat.py,v 1.23 2012/09/29 16:03:54 burnett Exp $
Author: T.Burnett <tburnett@uw.edu>
"""
import sys
import numpy as np
from . import bandlike, sourcelist, prior

         
class ROIstat(object):
    """ manage statistical analysis of an ROI
    
   
    Contains two lists:
       * all sources, in self.sources
         order is, for now, the same as for ROIAnlysis
       * a list of BandLike objects, one per band, attribute all_bands
          each manages the computation of the likelihood for its band
    Notes:
        * This is not intended to be a replacement for ROIAnalysis: complete user-level
          functionality is managed by a subclass, see main.py
        * the constructor takes an object containing a list of ROIBand objects and models, 
        * computation of the likelihood and its derivative, the basic task for 
          this class, is easily limited to a subset of the original set of bands,
          see select_bands
        * localization and extended source fits are handled by clients of this class
          by specifying a single source to be variable with the freelist parameter for
          initialize and setting reset=True in the update method 
    """
    
    def __init__(self, roi, bandsel=lambda b: True, **kwargs):
        """
        roi : ROIsetup object
            use to setup sources, find bands
        bandsel : function object returns bool
            returns True for a selected band
        """
        self.roi_setup = roi # save to allow creating a clone
        self.name = roi.name
        self.roi_dir = roi.roi_dir
        self.exposure = roi.exposure
        #try:
        self.sources = sourcelist.SourceList(roi) 
        #except Exception,msg:
        #    t = 'Failed to create ROI %s: %s' % (self.name, msg)
        #    raise Exception(t) 
            
        quiet = kwargs.pop('quiet', False)
        self.all_bands = bandlike.factory(filter(bandsel, roi.bands), self.sources , roi.exposure, quiet=quiet)
        self.selected_bands = self.all_bands # the possible subset to analyze
        self.saved_bands = None
        self.calls=0
        self.call_limit=1000
        self.quiet=quiet
        #self.prior= prior.ModelPrior(self.sources, self.sources.free)
        #self.prior.enabled=kwargs.pop('enable_prior', False)
        self.saved_pars = self.get_parameters()
    
    def __str__(self):
        return 'ROIstat for ROI %s: %d bands %d sources (%d free)' \
                % (self.name, len(self.all_bands), len(self.sources), sum(self.sources.free))
    @property
    def parameter_names(self):
        return self.sources.parameter_names
    def get_parameters(self):
        return self.sources.get_parameters()
    def set_parameters(self,parameters):
        self.sources.set_parameters(parameters)
        self.update()
    def get_external(self):
        return self.sources.get_parameters()
    def set_external(self, par):
        assert False, 'do I use this?'
        self.sources.set_parameters(np.log10(par))
    parameters = property(get_external, set_external, doc='array of free parameters')
    @property
    def model_parameters(self):
        return self.sources.model_parameters
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
        for band in self.all_bands:
            band.initialize(freelist)
        #self.prior.initialize(freelist)
        self.update()
        self.saved_pars = self.get_parameters()
        
    def reset(self):
        """ restore original set of parameters"""
        self.set_parameters(self.saved_pars)
        self.update()

    def select_bands(self, bandsel= None, emin=None): 
        """ select a subset of the bands for analysis
        parameteters
        ------------
            bandsel : function of a ROIBand that returns bool, like lambda b: b.e<200
        To restore, call with no arg
        Note that one could also replace selected_bands with a subset of all_bands
        """
        if bandsel is not None:
            self.selected_bands = np.array([bs for bs in self.all_bands if bandsel(bs.band)])
            if not self.quiet:print 'selected subset of %d bands for likelihood analysis' % len(self.selected_bands)
        elif emin is not None:
            self.selected_bands = np.array([bs for bs in self.all_bands if bs.band.emin>=emin])
            if not self.quiet:
                print 'selected subset of %d/%d bands for likelihood analysis'\
                    % (len(self.selected_bands), len(self.all_bands))
        else:
            self.selected_bands =  self.all_bands
            print '*** restoring %d bands ***' % len(self.all_bands)
        
    
    def select_source(self, sourcename):
        """ 
        Select a specific source, or all
        paramaters
            sourcename : string or None
                if None, all variable sources are selected
                otherwise will compute likelihood only for the given source
        Note: this does not change the set of variable sources defined by the SourceList: but should.
        It should also require that the source has variable parameters
        """
        if sourcename is None:
            self.initialize()
            return None
        self.source_mask = np.array([source.name==sourcename for source in self.sources])
        if sum(self.source_mask)!=1:
            raise Exception('select_source: source %s not found'%sourcename)
        self.initialize(self.source_mask)
        return np.array(self.sources)[self.source_mask][0]

    def update(self, reset=False):
        """ perform update on all selected bands, and variable sources
        if reset is True, assume that must also reinitialize angular dependence of sources
        """
        if reset:
            assert sum(self.source_mask)==1, 'expect only one source selected'
        map(lambda s: s.update(reset), self.selected_bands)
        #self.prior.update()
        
    def log_like(self):
        """ return sum of log likelihood for all bands
        """
        self.update()
        return sum([blike.log_like() for blike in self.selected_bands]) #+ self.prior.log_like()
        
    def __call__(self, par):
        """ (negative) log likelihood as a function of the free parameters par 
        appropriate for minimizing
        """
        self.set_parameters(par)
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
        return np.array([blike.gradient() for blike in self.selected_bands]).sum(axis=0)
        
        
    def chisq(self):
        return sum([blike.chisq() for blike in self.selected_bands])
        
    def dump(self, **kwargs):
        map(lambda bs: bs.dump(**kwargs), self.selected_bands)

    def select_band(self, energy, etype):
        """ Select and return a BandLike object based on energy and conversion type
        """
        # note: could be made a lot more efficient if needed
        sel = filter(lambda b:etype==b.band.ct and energy>=b.band.emin and energy<b.band.emax ,
            self.all_bands)
        assert len(sel)==1, 'Band with energy %s, type %s not found' % (energy, etype)
        return sel[0]
        
        