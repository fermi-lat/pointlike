"""
Manage spectral and angular models for an energy band to calculate the likelihood, gradient
   Currently delegates some computation to classes in modules like.roi_diffuse, like.roi_extended
   
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/bandlike.py,v 1.35 2013/11/18 00:03:24 burnett Exp $
Author: T.Burnett <tburnett@uw.edu> (based on pioneering work by M. Kerr)
"""

import sys, types
import numpy as np
from scipy import misc
from  uw.utilities import keyword_options

    
class BandLike(object):
    """ manage the likelihood calculation for a band 
    """
    defaults = (
        ('quiet', True, 'set False for info'),
        )
    @keyword_options.decorate(defaults)
    def __init__(self, band, sources, free=None, **kwargs):
        """
           band    : ROIband object
           sources : list of sources.Source objects
           free    : [array of bool | None]
                to select models with variable parameters
                If None, select all
        """
        keyword_options.process(self, kwargs)
        # make a list of the Response objects
        self.bandsources = np.array(map(lambda s: s.response(band, quiet=self.quiet), sources))

        self.band = band 
        self.exposure_factor = band.exposure.correction
        self.data = band.pix_counts  if band.has_pixels else []# data from the band
        self.pixels=len(self.data)
        
        self.initialize(free)
        self.update()
        # special code to unweight if galactic diffuse too large
        self.unweight = self.make_unweight()
         
    def make_unweight(self):
        """ return an unweighting factor <=1.0 to use to multiply the log likelihood
        assume that first BandSource object has the galactic diffuse
        
        systematic : float
            a fraction representing the relative systematic uncertainty in the galactic diffuse
        """
        try:
            if hasattr(self['ring'].source, 'systematic'): 
                systematic = self['ring'].source.systematic
            else:  return 1.0
        except: return 1.0
        if systematic==0: return 1.0 
        n = 1/systematic**2
        # m is the number of counts from the galactic diffuse in the footprint of a point source
        m = self['ring'].counts / (self.band.psf(0)[0]*self.band.solid_angle)
        u = min(1., n/m)
        return u

    def __str__(self):
        b = self.band
        return '%s.%s: %d models (%d free) applied to band %.0f-%.0f, %s with %d pixels, %d photons'\
                % (self.__module__,self.__class__.__name__,len(self.bandsources), sum(self.free), b.emin, b.emax, 
                 ('front back'.split()[b.event_type]), self.pixels, sum(self.data), )
                 
    def __repr__(self): return self.__str__()
    
    def __getitem__(self, i): 
        """ return a BandSource object refererence, either by index, or by source name"""
        if type(i)==types.StringType:
            t = list(self.bandsources)
            for bs in t:
                if i==bs.source.name:
                    return bs
            raise Exception('Source "%s" not found in band sources' %i)
        return self.bandsources[i]
        
    def initialize(self, free):
        """ should only call if free array changes.
            Saves the combined prediction from the models with fixed parameters
        """
        self.free = free if free is not None else np.ones(len(self.bandsources), bool)
        self.free_sources = self.bandsources[self.free]
        self.counts = self.fixed_counts = sum([b.counts for b in self.bandsources[ ~ self.free]])
        if not self.band.has_pixels: return # no data, should quit?
        self.fixed_pixels = np.zeros(self.pixels)
        for m in self.bandsources[ ~ free]:
            self.fixed_pixels += m.pix_counts
        self.model_pixels = self.fixed_pixels.copy()
        for m in self.free_sources:
            self.model_pixels += m.pix_counts
        
    def update(self, reset=False, **kwargs):
        """ assume that parameters have changed. Update only contributions 
        from models with free parameters. *must* be called before evaluating likelihood.
        reset: bool, optional
            if True, need to reinitialize variable source(s), for change of position or shape
        """
        if self.band.has_pixels: self.model_pixels[:]=self.fixed_pixels
        self.counts = self.fixed_counts
        force = kwargs.get('force', False) # set to defeat
        for bandsource in self.free_sources:
            if reset: 
                bandsource.initialize()
                bandsource.source.changed=False
            elif bandsource.source.changed or force:
                bandsource.update()
            if self.band.has_pixels: self.model_pixels += bandsource.pix_counts
            self.counts+= bandsource.counts
        self.weights = self.data / self.model_pixels

    def log_like(self):
        """ return the Poisson extended log likelihood """
        try:
            pix = np.sum( self.data * np.log(self.model_pixels) )  if self.pixels>0 else 0
            w = pix - self.counts * self.exposure_factor
            return self.unweight * w
        except FloatingPointError, e:
            print '%s: Floating point error %s evaluating likelihood for band %s' % (e,self,)
            raise
    
    def gradient(self):
        """ gradient of the likelihood with resepect to the free parameters
        """
        if not self.band.has_pixels: return None
        return self.unweight * np.concatenate(
                [m.grad(self.weights, self.exposure_factor) for m in self.free_sources]
            )
       
    def model_counts(self, sourcemask=None):
        """ return the model predicted counts for all or a subset of the sources
        sourcemask : array of bool
            select a set of sources. In order of sources from factory
        """
        if sourcemask is not None:
            assert len(sourcemask)==len(self), 'bad input to model_counts'
        t = np.array([s.counts for s in self])
        return sum(t) if sourcemask is None else sum(t[sourcemask])

    def add_source(self, source):
        """ add a new source 
            source: sources.Source object
        """
        # ugly but compact
        t = list(self.bandsources)
        t.append(source.response(self.band))
        self.bandsources = np.array(t)
        
    def del_source(self, source):
        """ remove the source """
        t = list(self.bandsources)
        for bs in t:
            if source.name==bs.source.name:
                t.remove(bs)
                self.bandsources = np.array(t)
                return
        raise Exception('source "%s" not found to delete' % source.name)
    
    def fill_grid(self, sdirs):
        """ fill a grid with values, which are counts/sr, so must be multiplied by the pixel size
        """
        t = np.zeros(len(sdirs))
        for m in self.bandsources:
            t+= m.fill_grid(sdirs)
        return t
        
    def dataframe(self, **kw):
        """ return a pandas.DataFrame for diagnostics """
        import pandas as pd
        df = pd.DataFrame(dict([(s.source.name, 
                dict(counts=s.counts.round(), 
                    overlap=s.overlap, 
                    free=self.free[i],
                    extended=s.source.skydir is not None and hasattr(s.source, 'dmodel'),
                    diffuse=s.source.skydir is None,
                    distance=np.degrees(s.band.skydir.difference(s.source.skydir)) if s.source.skydir is not None else 0 ),
                    )
                for i,s in enumerate(self)])
            ).T
        return df

    def counts_in_pixel(self, source_index, skydir):
        """ return a tuple of predicted signal and background counts in the pixel corresponding to skydir
        Note that if the pixel has no data, it will not have been computed for the model; instead
        this will return the average background
 
        source_index : int
            the index of the source
        skydir : SkyDir object
        """
        from skymaps import WeightedSkyDir
        band = self.band
        if not hasattr(self, 'pixel_indeces'): 
            self.pixel_indeces = list([band.b.index(x) for x in band.wsdl])
        source = self[source_index]
        hp_index = band.b.index(skydir)
        try:
            pixel_index = self.pixel_indeces.index(hp_index)
            signal = source.counts * source.pixel_values[pixel_index]
            back  = self.model_pixels[pixel_index] - signal
        except:
            # missing pixel; calculate the expected source counts
            # and estimate total by mean of model in ROI
            signal = source.flux_value(skydir) * source.counts
            back = self.model_pixels.mean()
        return signal, back 
       
         
class BandLikeList(list):
    
    def __init__(self, roi_bands, roi_sources):
        """ create a list, one per band, of BandLike objects for a given ROI 
        Provide some useful methods 
        
        parameters
        ----------        
        roi_bands : list of ROIBand objects
        roi_sources :  sourcelist.SourceList object
            a list of Source objects 
        """

        self.sources = roi_sources
        self.bands = roi_bands
        for band in roi_bands:
            self.append( BandLike(band, self.sources, self.sources.free) )
    @property
    def free_list(self):
        return self.sources.free
    
    def __repr__(self):
        return '%s.%s: \n\t%s\n\t%s\n\tparameters: %d/%d free' % (self.__module__,
            self.__class__.__name__,
            self.sources, self.bands, sum(self.free_list), len(self.free_list))

    def initialize(self):
        assert False, 'needed?'
        
    def log_like(self):
        return sum( b.log_like() for b in self)
        
    def update(self, **kwargs):
        for b in self: 
            b.update(**kwargs)
        self.sources.parameters.clear_changed()
        
    def gradient(self):
        return np.array([blike.gradient() for blike in self]).sum(axis=0) 
        
    def hessian(self, delta=1e-6):

        """ return a hessian matrix based on the current parameters
        For sigmas and correlation coefficients:
                cov =  self.hessian().inv()
                sigs = np.sqrt(cov.diagonal())
                corr = cov / np.outer(sigs,sigs)

        """
        parameters = self.sources.parameters
        gzero = self.gradient()
        fzero = self.log_like()
        parz = parameters.get_all()
        def dg(i, delta=delta):
            parameters[i] = parz[i]+delta
            self.update()
            fl = self.log_like()
            assert abs(fl- fzero)>1e-8, '%d %.3e' % (i, fl - fzero)
            gprime = self.gradient()
            ret= (gprime-gzero)/(2*delta)
            parameters[i] = parz[i]
            self.update()
            gcheck = self.gradient()
            assert np.all(np.abs(gcheck-gzero)<1e-5), gcheck-gzero
            return ret
        
        cov = np.matrix(map(dg, range(len(parz))))
        return cov
        
    def likelihood_plots(self, index = None):
        """ one, or all likelihood plots
        """
        import matplotlib.pyplot as plt
        if index is None:
            n = len(self.parameters)
            fig, axx = plt.subplots(n/5,5, figsize=(15,12), sharex=True, sharey=True)
            for i, ax in enumerate(axx.flatten()):
                self.likelihood_functor(i).plot(ax = ax, nolabels=True)
        else:
            fig = self.likelihood_functor(index).plot()
        return fig
        
    
    def likelihood_functor(self, indexlist, force=False ):
        """return a functor of one variable, for testing at the moment"""
        blike = self
        class LikelihoodFunctor(object):
            def __init__(self, indexlist, force):
                self.parameters = blike.sources.parameters
                self.aup = force
                self.k = indexlist
                
            def setpar(self,par):
                self.parameters[self.k] = par
                blike.update(force=self.aup)
                       
            def __call__(self, par):
                self.setpar( par )
                ret = blike.log_like()
                return ret
                
            def gradient(self, par, full=False):
                self.setpar(par)
                return blike.gradient()[self.k] if not full else blike.gradient()
            
            def get_quad(self, par, dx=1e-5):
                """estimate position of peak and sigma using scipy derivative, assuming quadratric
                """
                d1,d2 = [misc.derivative(self, par, n=n, dx=dx) for n in (1,2)]
                sig = np.sqrt(-2./d2)
                pmax = par - d1/d2
                return pmax, sig
            
            def derivative(self, par, n=1, dx=0.0001):
                """ numerical derivative; make sure state is same after"""
                pz = self.parameters[self.k]
                before = self(pz)
                d = misc.derivative(self, par, dx=dx, n=n)
                self.parameters[self.k] = pz
                assert self(pz)==before, 'failed to restore?'
                return d
                
            def plot(self, ax=None, nolabels=False , y2lim=(-10,10)):
                """make a plot showing the log likelihood and its derivative as a function of
                expected sigma, evaluated from the second derivative at the current point
                """
                import matplotlib.pyplot as plt
                func = self
                index=self.k
                pz =self.parameters[self.k]
                x0, sig = func.get_quad(pz)
                ref = func(x0)
                if ax is None:
                    fig, ax = plt.subplots( figsize=(3,3))
                else: fig = ax.figure
                plt.subplots_adjust(wspace=0.3)
                xsig = np.linspace(-3, 3, 27)
                x =  x0 + xsig * sig 
                ax.plot(xsig, map(func,x)-ref, '-')
                ax.plot(xsig, -((x-x0)/sig)**2, '--')
                ax.plot((pz-x0)/sig, func(pz)-ref, 'db')
                plt.setp(ax, ylim=(-9,0.5))
                if not nolabels: ax.set_ylabel('log likelihood')
                ax.set_title('#%d: %s' %(index,blike.sources.parameter_names[index]), size=10)
                ax.text( -2,-8, 'par %.3f\nsig %.3f' % (pz,sig), size=10)
                ax.axvline(0, color='k', ls = ':')
                
                ax2 = ax.twinx()
                gradvals = sig*np.array(map(func.gradient, x))
                ax2.plot(xsig, gradvals, '-r')
                ax2.axhline(0, color='r', ls=':')
                ax2.set_ylim( y2lim)
                if not nolabels: ax2.set_ylabel('derivative (sig units)')
                else: ax2.set_yticklabels([])
                
                self.setpar(pz) # restore when done
                return fig
                
                
        return LikelihoodFunctor(indexlist, force)