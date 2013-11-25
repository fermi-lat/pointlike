"""
Manage spectral and angular models for an energy band to calculate the likelihood, gradient
   
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/bandlike.py,v 1.42 2013/11/24 17:45:02 burnett Exp $
Author: T.Burnett <tburnett@uw.edu> (based on pioneering work by M. Kerr)
"""

import sys, types
import numpy as np
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
        Methods to calculate the likelihood, gradient, and hessian
        
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
            
        self.set_selected(self)# set selected for a subset?
        self.all_energies = self.energies[:]
        
    #  the band selection mechanism, used by log_like, update and gradient   
    def set_selected(self, values):
        """setter for the property selected, which must be a subset of self"""
        if values in self: # single guy
            self._selected = [values]
            return
        assert set(values).issubset(self), 'Improper selection'
        self._selected = values
    def get_selected(self):
        return self._selected
    selected = property(get_selected, set_selected)
        
    def select(self, index=None, event_type=None):
        """ Select an energy band or bands
        parameters:
        ----------
        index: None or integer
            an index into the list of energies; if None, select all bands
            and use the current spectral model, otherwise a powerlaw to 
            represent an model-independent flux over the band.
        event_type : None or integer
            if None, select both front and back, otherwise 0/1 for front/back
                
        """
        if index==None: #select all (initially selected) bands, use input model
            selected_bands = self[:]
        else:
            energy = self.all_energies[index]
            type_select = lambda x : True if event_type is None else x==event_type
            selected_bands = filter(lambda b: abs(b.band.energy-energy)<1 and type_select(b.band.event_type), self)
            if len(selected_bands)==0:
                raise Exception( 'did not find any bands for energy %.1f: %s are available' %( energy, self.energies))
        self.selected = selected_bands

    @property
    def free_sources(self):
        """ list of sources with currently free parameters
        """
        return self.sources.free
    
    @property
    def energies(self):
        return  np.sort(list(set([ sm.band.energy for sm in self._selected])))
    @property
    def emin(self):
        return np.array([b.band.emin for b in self._selected]).min()
    @property
    def emax(self):
        return np.array([b.band.emax for b in self._selected]).max()
        
    def __repr__(self):
        sel = '%d bands'%len(self.bands) if len(self.selected)==len(self) else '%d / %d selected bands'\
            %(len(self.selected),len(self))
        return '%s.%s: \n\t%s\n\t%s\n\tParameters: %d in %d/%d free sources' % (
            self.__module__, self.__class__.__name__,
            self.sources, sel, len(self.sources.parameters), 
            sum(self.free_sources),  len(self.free_sources))

    def initialize(self, sourcename):
        """ initialize the specifed source in  all selected bands
            and update the band"""
        for b in self._selected:
            b[sourcename].initialize()
            b.update()
        
    # the following methods sum over the current set of bands
    def log_like(self, summed=True):
        """log likelihood for current set of bands
        summed : bool, optional
        if false, return the array of likelihods for each band
        """
        r = np.array([b.log_like() for b in self._selected])
        return  sum(r) if summed else r
        
    def update(self, **kwargs):
        for b in self._selected: 
            b.update(**kwargs)
        self.sources.parameters.clear_changed()
        
    def gradient(self):
        return np.array([blike.gradient() for blike in self._selected]).sum(axis=0) 
        
    def hessian(self, mask=None, delta=1e-6):
        """ return a hessian matrix based on the current parameter set
        This makes a numerical derivative of the analytic gradient, so not exactly
        symmetric, but the the result must be (nearly) symmetric.
        
        mask : [None, array of bool]
            If present, must have dimension of the parameters, will generate a sub matrix
        
        For sigmas and correlation coefficients, invert to covariance
                cov =  self.hessian().I
                sigs = np.sqrt(cov.diagonal())
                corr = hess / np.outer(sigs,sigs)
        """
        # get the source parameter management object
        parameters = self.sources.parameters
        parz = parameters.get_parameters()
        if mask is None: mask = np.ones(len(parz),bool)
        else:
            mask = np.asarray(mask)
            assert len(mask)==len(parz)
        # initial values for the likelihood and gradient
        fzero = self.log_like()
        glast = gzero = self.gradient()[mask]
        t = []
        for i in np.arange(len(parz))[mask]:
            # increment current variable and get new gradient
            parameters[i] = parz[i]+delta
            self.update()
            gnow = self.gradient()[mask]
            # numerical derivative of gradient with respect to this parameter
            t.append( (gnow-glast)/(2*delta))
            glast = gnow
        hess = np.matrix(t)
        parameters.set_parameters(parz) #restore all parameters, check that no difference
        self.update()
        assert abs(fzero-self.log_like())<1e-2
        return hess 
       
       
    
