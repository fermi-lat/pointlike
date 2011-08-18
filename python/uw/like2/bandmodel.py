"""
Manage spectral and angular models for an energy band
   Delegates computation to classes in modules like.roi_diffuse, like.roi_extended, like.roi_bands

$Header$
Author: T.Burnett <tburnett@uw.edu> (based on pioneering work by M. Kerr)
"""

import numpy as np
from uw import like # for pypsf

class BandModel(object):
    """ superclass for point or diffuse band models, used to implement printout"""

    def __str__(self):
        return '%s for source %s (%s), band %.0f-%.0f, %s' \
            % (self.__class__.__name__, self.source.name, self.model.name, self.band.emin, self.band.emax, 
                ('front back'.split()[self.band.b.event_class()]),)

    def dump(self, out=None):
        """ useful summary dump """
        i = self.source_index
        b = self.band
        print >>out, self.__str__()
        print >>out, 'exposure ratio, overlap %.3f %.3f:'%( b.er[i], b.overlaps[i])
        pc = self.pix_counts
        print >>out, ('pixel counts: min, max, sum: '+3*'%8.1f') % ( pc.min(), pc.max(), pc.sum())
        print >>out, 'total counts %8.1f'%  self.counts # b.ps_counts[i]

class BandDiffuse(BandModel):
    """ 
        Use a ROIDiffuseModel or ROIExtendedModel to compute the expected distribution of pixel counts
        for ROIBand
    """
    def __init__(self,  band, source): 
        """
            band : a like.roi_bands.ROIBand object reference
                depends on: emin,emax,e,radius_in_rad,wsdl,has_pixels,
                    solid_angle,pixelArea(), 
                sets:    bg_counts[source_index], bg_pix_counts[:,source_index]
            source : generalized Source with factory and manager_index
        """
        
        self.source_index = source.manager_index 
        self.source = source.factory(self.source_index) 
        self.model  = self.source.smodel
        self.source.quiet = True
        self.band = band
        self.pix_counts = band.bg_pix_counts[:,self.source_index]
        self.initialize()
        self.update()
    
    def initialize(self):
        """ establishs a state from
            which the model counts can be calculated.  E.g., evaluating
            the model on a series of energy subplanes and storing the
            results for a later Simpson's rule.
        """
        self.source.initialize_counts([self.band])
        
    def update(self):
        """ sets the following members of the band:
            bg_counts[source_index] -- the total counts expected for
                    the model in the aperture
            bg_pix_counts[:,source_index] -- if the band has pixels,
                    an npixel vector with the expected counts from the
                    model for each data pixel
        """
        self.source.update_counts([self.band], self.source_index)
        self.pix_counts = self.band.bg_pix_counts[:,self.source_index]
        self.counts = self.band.bg_counts[self.source_index]
        
    def grad(self, weights, phase_factor=1):
        """ contribution to the gradient
        """
        nfree =np.sum(self.model.free)
        if nfree==0: return []
        if (nfree==1) and self.model.free[0]:
            # special case -- only normalization free 
            apterm  =  phase_factor*self.counts
            pixterm =  (self.pix_counts*weights).sum() 
            g = [(apterm - pixterm)/self.model.getp(0)]
        else:
            # this is only rarely needed and not carefully tested
            self.band.pix_weights = weights
            self.band.phase_factor = phase_factor
            g = self.source.gradient([self.band], self.source_index)
        return g
    
class BandPoint(BandModel):
    """ Compute the expected distribution of pixel counts from a point source in a ROIBand
    """
    def __init__(self, band, source ): 
        """
            point_source_factory : function that returns point source, roi_dir
            band : ROIBand object
                reads: e, exp, wsdl, psf, b.pixelArea()
                calls: expected
                
        """
        self.source, self.roi_dir = source.factory(source.manager_index) 
        self.model = self.source.model
        self.band = band
        self.pix_counts = np.zeros(band.npix) 
        self.initialize()
        self.update()
       
    def initialize(self):
        band = self.band
        ps = self.source
        roi_dir = self.roi_dir
 
        energy,exposure = band.e, band.exp.value

        # make a first-order correction for exposure variation
        denom = exposure(roi_dir,energy)
        if denom==0:
            raise Exception('BandPoint: exposure is zero for  energy %f'%en)
        self.er = exposure(ps.skydir,energy)/denom 
        
        #unnormalized PSF evaluated at each pixel
        rvals  = np.empty(len(band.wsdl),dtype=float)
        band.psf.cpsf.wsdl_val(rvals, ps.skydir, band.wsdl)
        self.ps_pix_counts =rvals*band.b.pixelArea()

        # note this will use a numerical integration if the ragged edge is impt.
        overlap_function = like.pypsf.PsfOverlap()
        self.overlap = overlap_function(band,roi_dir,self.source.skydir)  

    def update(self):  
        t= self.band.expected(self.source.model) * self.er
        self.pix_counts = self.ps_pix_counts * t
        self.counts = t*self.overlap

    def grad(self, weights, phase_factor=1): 
        """ contribution to the overall gradient
        """
        b = self.band
        model = self.model
        if np.sum(model.free)==0 : return []
        # Calculate the gradient of a spectral model (wrt its parameters) integrated over the exposure.
        ##TODO: encapsulate this simple Simpsons rule integral
        g = (model.gradient(b.sp_points)*b.sp_vector).sum(axis=1)[model.free] * self.er
        apterm = phase_factor* self.overlap
        pixterm = (weights*self.ps_pix_counts).sum() if b.has_pixels else 0
        return g * (apterm - pixterm)
   
    
class BandModelStat(object):
    """ manage the likelihood calculation for a band 
    """
    def __init__(self, band, bandmodels, free):
        """
           band : ROIband object
           bandmodels : list of BandModel objects associated with this band
           free : array of bool to select models
        """
        self.bandmodels = bandmodels
        self.band = band # for reference: no top-level access
        self.phase_factor = band.phase_factor
        #for b in bandmodels: self.append(b)
        self.data = band.pix_counts
        self.pixels=len(self.data)
        self.initialize(free)
        self.update()
        
    def __str__(self):
        b = self.bandmodels[0].band
        return 'BandModelStat: %d models (%d free) applied to band %.0f-%.0f, %s with %d pixels, %d photons'\
                % (len(self.bandmodels), sum(self.free), b.emin, b.emax, 
                 ('front back'.split()[b.b.event_class()]), self.pixels, sum(self.data))
        
    def initialize(self, free):
        """ only call later if free array changes
        """
        self.free = free
        self.freemodels = self.bandmodels[free]
        self.fixed = np.zeros(self.pixels)
        for m in self.bandmodels[-free]:
            self.fixed += m.pix_counts
        self.counts=self.fixed_counts = sum([b.counts for b in self.bandmodels[-free]])
        self.model = self.fixed.copy()
        for m in self.freemodels:
            self.model += m.pix_counts
        
    def update(self):
        """ assume that parameters have changed
        """
        self.model[:]=self.fixed
        self.counts = self.fixed_counts
        for m in self.freemodels:
            m.update()
            self.model += m.pix_counts
            self.counts+= m.counts

    def log_like(self):
        """ return the Poisson extended log likelihood """
        pix = np.sum( self.data * np.log(self.model) ) if self.pixels>0 else 0
        return pix - self.counts*self.phase_factor

    def chisq(self):  
        """ compute the chi squared
        """
        M = self.model        # array of model predictions 
        D = self.data         # array of data counts
        return np.sum((M-D)**2/D)
    
    def gradient(self):
        """ gradient of the likelihood with resepect to the free parameters
        """
        weights = self.data/self.model
        return np.concatenate([m.grad(weights, self.phase_factor) for m in self.bandmodels])
       
    def dump(self, **kwargs):
        map(lambda bm: bm.dump(**kwargs), self.bandmodels)
    
  

def factory(bands, sources):
    """ return an array, one per band, of BandModelStat objects 
        bands : list of ROIBand objects
        sources : sourcelist.SourceList object
            a list of Source objects, which must contain a factory attribute 
    """
    def bandmodel_factory(band, source):
        """helper factory that returns a BandModel object appropriate for the source"""
        factory_name = source.factory.__class__.__name__
        assert factory_name in ('PointSourceFactory', 'ConvolvedSourceFactory'),\
            'source factory name %s not recognized' % factory_name 
        B = BandPoint if factory_name=='PointSourceFactory' else BandDiffuse
        return B(band, source)
        
    return np.array(
        [ BandModelStat(band, 
                    np.array( [bandmodel_factory(band, source) for source in sources]),
                    sources.free)
            for band in bands])
         