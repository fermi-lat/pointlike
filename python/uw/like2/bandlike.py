"""
Manage spectral and angular models for an energy band to calculate the likelihood, gradient
   Currently delegates some computation to classes in modules like.roi_diffuse, like.roi_extended
   
classes:
    BandSource -- superclass for basic Band/Source association
        BandDiffuse -- diffuse
        BandExtended  -- extended
        BandPoint  -- point source
    BandLike -- manage likelihood calculation, using a list of BandSource objects
functions:
    factory -- create a list of BandLike objects from bands and sources

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/bandlike.py,v 1.7 2011/10/20 21:41:27 burnett Exp $
Author: T.Burnett <tburnett@uw.edu> (based on pioneering work by M. Kerr)
"""

import sys
import numpy as np

class BandSource(object):
    """ superclass for point or diffuse band models, used to implement printout
    subclasses implement code to compute prediction of the source model for the pixels in the band
    """
    def __init__(self, band, source ): 
        """
            band : a like.roi_bands.ROIBand object reference
                depends on: emin,emax,e,radius_in_rad,wsdl,
                    solid_angle,pixelArea(), 
            source : generalized Source 
        """
        self.source= source
        self.band = band
        self.initialize()
        self.update()

    def __str__(self):
        return '%s for source %s (%s), band (%.0f-%.0f, %s)' \
            % (self.__class__.__name__, self.source.name, self.spectral_model.name, self.band.emin, self.band.emax, 
                ('front back'.split()[self.band.b.event_class()]),)

    def dump(self, out=None):
        """ useful summary dump """
        b = self.band
        print >>out, self.__str__()
        if hasattr(self,'exposure_ratio'):
            print >>out, '\texposure ratio, overlap: %.3f %.3f'%( self.exposure_ratio, self.overlap)
        pc = self.pix_counts
        print >>out, ('\tpixel counts: min, max, sum: '+3*'%8.1f') % ( pc.min(), pc.max(), pc.sum())
        print >>out, '\ttotal counts %8.1f'%  self.counts 

    @property 
    def spectral_model(self):  return self.source.spectral_model
    
    def update(self, fixed=False):
        """ use current spectral model to update counts, pix_counts """
        b = self.band
        expected = b.exposure_integral(self.spectral_model) * self.exposure_ratio
        self.pix_counts = self.pixel_values * expected
        self.counts = expected*self.overlap

    def grad(self, weights, phase_factor=1): 
        """ contribution to the overall gradient
        """
        b = self.band
        model = self.spectral_model
        if np.sum(model.free)==0 : return []
        # Calculate the gradient of a spectral model (wrt its parameters) integrated over the exposure.
        g = b.exposure_integral(model.gradient, axis=1)[model.free] * self.exposure_ratio
        apterm = phase_factor* self.overlap
        pixterm = (weights*self.pixel_values).sum() if b.has_pixels else 0
        return g * (apterm - pixterm)


class BandDiffuse(BandSource):
    """  Apply diffuse model to an ROIband
    
        Use a ROIDiffuseModel to compute the expected 
        distribution of pixel counts  for ROIBand
        Convolving is done with 
        
    """
    def __init__(self,  band, source): 
        """
            band : a like.roi_bands.ROIBand object reference
                depends on: emin,emax,e,radius_in_rad,wsdl,
                    solid_angle,pixelArea(), 
            source : diffuse source with a make_grid method
        """
        
        self.source = source
        self.band = band
        self.energy = band.e
        self.setup=False
        self.initialize()
        self.update()
        
    @property 
    def spectral_model(self):  return self.source.diffuse_source.smodel

    def initialize(self):
        """ establish a state from
            which the model counts can be calculated.  E.g., evaluating
            the model on a series of energy subplanes and storing the
            results for a later Simpson's rule.
        """
        assert self.setup==False, 'attempt to re-initialize global source'
        self.setup=True
        band = self.band
        delta_e = band.emax - band.emin
        self.grid = self.source.make_grid(band.emin, band.ct)
        self.ap_evals = self.grid.ap_average(band.radius_in_rad) * band.solid_angle * delta_e

        self.pi_evals = self.grid(band.wsdl, self.grid.cvals) if band.has_pixels else []
        self.pi_evals *= band.b.pixelArea() * delta_e

        
    def update(self, fixed=False):
        """ sets the following attributes:
            counts -- the total counts expected for the model in the aperture
            pix_counts-- if the band has pixels, an npixel vector with the 
                        expected counts from the model for each data pixel
            actually use the ROIDiffuseModel object to modify the band, then we 
            copy them here. 
        """
        self.counts = self.ap_evals * self.spectral_model(self.energy)
        self.pix_counts = self.pi_evals * self.spectral_model(self.energy) if self.band.has_pixels else []

        
    def grad(self, weights, phase_factor=1):
        """ contribution to the gradient
        """
        model = self.spectral_model 
        nfree =np.sum(model.free)
        if nfree==0: return []
        if (nfree==1) and model.free[0]:
            # special case -- only normalization free 
            apterm  =  phase_factor*self.counts
            pixterm =  (self.pix_counts*weights).sum() 
            g = [(apterm - pixterm)/model.getp(0)]
        else:
            # this is only rarely needed and not carefully tested
            self.source.smodel = self.spectral_model
            self.band.pix_weights = weights
            self.band.phase_factor = phase_factor
            g = self.source.gradient([self.band], self.source_index)
        return g

 
class BandPoint(BandSource):
    """ Compute the expected distribution of pixel counts from a point source in a ROIBand
    """
       
    def initialize(self):
        """ setup values that depend on source position:
            exposure_ratio, overlap, pixel values
        """
        band = self.band
        self.exposure_ratio = band.exposure(self.source.skydir)/band.exposure(self.band.sd)
        self.overlap = band.psf_overlap(self.source.skydir)  
        
        #unnormalized PSF evaluated at each pixel
        self.pixel_values = band.pixels_from_psf(self.source.skydir)

class BandExtended(BandPoint):
    """  Apply extended model to an ROIband
    
        Use a ROIExtendedModel to compute the expected distribution of 
        pixel counts 
        Convolving is done with like.roi_extended.ExtendedSource 
    """
    def initialize(self):
        self.source.quiet = True # turn off convolvling messages
        self.source.initialize_counts([self.band])
        myband = self.source.bands[0]
        self.exposure_ratio, self.overlap = myband.er, myband.overlaps 
        self.pixel_values = myband.es_pix_counts
 
    
class BandLike(object):
    """ manage the likelihood calculation for a band 
    """
    def __init__(self, band, bandsources, free):
        """
           band : ROIband object
                only used here to get data and phase_factor
           bandsources : list of BandSource objects associated with this band
           free : array of bool to select models
        """
        self.bandsources = bandsources
        self.band = band # for reference: not used after this 
        self.energy = band.e # also for convenience: center of band
        self.phase_factor = band.phase_factor #this should be global
        self.data = band.pix_counts  if band.has_pixels else []# data from the band
        self.pixels=len(self.data)
        self.initialize(free)
        self.update()
        
    def __str__(self):
        b = self.bandsources[0].band
        return 'BandLike: %d models (%d free) applied to band %.0f-%.0f, %s with %d pixels, %d photons'\
                % (len(self.bandsources), sum(self.free), b.emin, b.emax, 
                 ('front back'.split()[b.b.event_class()]), self.pixels, sum(self.data))
                 
    def __getitem__(self, i): return self.bandsources[i]
        
    def initialize(self, free):
        """ should only call if free array changes.
            Saves the combined prediction from the models with fixed parameters
        """
        self.free = free
        self.free_sources = self.bandsources[free]
        self.fixed_pixels = np.zeros(self.pixels)
        for m in self.bandsources[ ~ free]:
            self.fixed_pixels += m.pix_counts
        self.counts=self.fixed_counts = sum([b.counts for b in self.bandsources[ ~ free]])
        self.model_pixels = self.fixed_pixels.copy()
        for m in self.free_sources:
            self.model_pixels += m.pix_counts
        
    def update(self, reset=False, fixed=False):
        """ assume that parameters have changed. Update only contributions 
        from models with free parameters. *must* be called before evaluating likelihood.
        reset: bool
            if True, need to reinitialize variable source(s), for change of position or shape
        fixed : bool
            if True, will not update prediction, for band_ts use (not implemented??)
        """
        self.model_pixels[:]=self.fixed_pixels
        self.counts = self.fixed_counts
        for m in self.free_sources:
            if reset: m.initialize()
            m.update(fixed)
            self.model_pixels += m.pix_counts
            self.counts+= m.counts

    def log_like(self):
        """ return the Poisson extended log likelihood """
        pix = np.sum( self.data * np.log(self.model_pixels) ) if self.pixels>0 else 0
        return pix - self.counts*self.phase_factor 

    def chisq(self):  
        """ compute the chi squared  """
        return np.sum((self.model_pixels-self.data)**2/self.model_pixels)
    
    def gradient(self):
        """ gradient of the likelihood with resepect to the free parameters
        """
        weights = self.data/self.model_pixels
        return np.concatenate([m.grad(weights, self.phase_factor) for m in self.bandsources])
       
    def dump(self, **kwargs):
        map(lambda bm: bm.dump(**kwargs), self.bandsources)

 
def factory(bands, sources, quiet=False):
    """ return an array, one per band, of BandLike objects 
        bands : list of ROIBand objects
        sources : sourcelist.SourceList object
            a list of Source objects, which must contain a factory attribute 
    """
    def bandsource_factory(band, source):
        """helper factory that returns a BandModel object appropriate for the source"""
        class_name = source.__class__.__name__
        B = dict(PointSource=BandPoint, 
                DiffuseModelFromCache=BandDiffuse,
                IsotropicModel=BandDiffuse,
                ROIExtendedModelAnalytic=BandExtended,
                ROIExtendedModel=BandExtended)[class_name]
        return B(band, source)
    bandlist = []
    for i,band in enumerate(bands):
        bandlist.append(BandLike(band, 
                    np.array( [bandsource_factory(band, source) for source in sources]),
                    sources.free)
                    )
    if not quiet: print
    return np.array(bandlist)
         
