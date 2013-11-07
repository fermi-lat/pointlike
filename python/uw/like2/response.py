"""
Classes to compute response from various sources
 


$Header$
author:  Toby Burnett
"""
import os
import numpy as np
from . import convolution
from uw.utilities import keyword_options

class Response(object):
    def __init__(self, source, band, **kwargs):
        """
        source : PointSource object
            skydir 
            model
        band : ROIband object
            psf, exposure functions for the event type
            sd, radius : location, size of ROI
            emin, emax : energy range
            wsdl : list of pixel positions, perhaps None
            pixel_size : pixel solid angle
        """
        self.band=band
        self.source=source
        self.roicenter = self.band.sd
        self.initialize()
        
    def __repr__(self):
        return '%s.%s: ROI at %s source "%s" at %s' %( self.__module__,self.__class__.__name__,
            self.roicenter, self.source.name, self.source.skydir)

    def exposure_integral(self, skydir=None):
        if skydir is None: skydir = self.roicenter
        return self.band.exposure.model_integral(skydir, self.source.model, self.band.emin, self.band.emax)

class PointResponse(Response):
    """Manage predictions of the response of a point source
    
    Calculates:
        counts : the expected counts in an ROI
        pixel_values : arraay of counts in a set of pixels, presumably corresponding to pixels with data
        grad : gradient
        
    """

    def initialize(self):
        """ Only needs to be done if position changes
        """
        self.overlap = self.band.psf.overlap(self.roicenter, self.band.radius, self.source.skydir)
        self._exposure_ratio = self.band.exposure(self.source.skydir) / self.band.exposure(self.roicenter)
        self.evaluate()
        
    def evaluate(self): #, weights=None, exposure_factor=1):
        """ update values of counts, pixel_values, used for likelihood calculation, derivatives
        Called when source parameters change
        """
        self.counts = self.exposure_integral(self.source.skydir) * self.overlap
        
    #def evaluate_at(self, skydirs):
    #    """ set relative exposure values from a list of skydirs
    #    """
    #    return np.array([self.band.psf(s.difference(self.source.skydir))[0] for s in skydirs]) #* exp
    #    
    #def __call__(self, skydir):
    #    return self.band.psf(skydir.difference(self.source.skydir))# * self.band.exposure(self.source.skydir)
    # 

class IsotropicResponse(Response):

    defaults =(
        ('pixelsize', 0.25, 'Size of pixels to use for convolution grid'),
        ('npix',      61,   'Number of pixels (must be an odd number'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, source, band, **kwargs):
        keyword_options.process(self, kwargs)

        super(IsotropicResponse, self).__init__(source, band, **kwargs)
        
    def initialize(self):

        #set up the diffuse model
        self.dmodel = self.source.dmodel[self.band.event_type]
        self.dmodel.setEnergy(self.band.energy)
        
        # create a grid for evaluating counts integral over ROI, individual pixel predictions
        grid = self.grid= convolution.ConvolvableGrid(center=self.roicenter, 
                npix=self.npix, pixelsize=self.pixelsize)
        flux = self.dmodel(self.roicenter)
        center_exposure = self.band.exposure(self.roicenter) # flux at roicenter
        grid.bg_fill(self.band.exposure, None, cache=1)
        inside = grid.dists< np.radians(self.band.radius)
        mean_exposure = grid.bg_vals[inside].mean()
        grid.bg_vals /= mean_exposure
        self.factor = mean_exposure/center_exposure * flux * self.band.solid_angle
        self.evaluate()
        
    def evaluate(self):
        self.counts = self.exposure_integral() * self.factor
    
    def __call__(self, skydir):
        return super(IsotropicConvolver, self).__call__(skydir, self.bg_vals)

    def evaluate_at(self, skydirs):
        return super(IsotropicConvolver, self).__call__(skydirs, self.bg_vals)
        

    
    def convolve(self, energy=None):
        if energy is not None: self.band.set_energy(energy)
        self.bg_fill(self.band.exposure, None, cache=self.dmodel(self.center))
        inside = self.dists< np.radians(self.band.radius)
        self.ap_average = self.bg_vals[inside].mean()
    
    def __call__(self, skydir):
        return super(IsotropicConvolver, self).__call__(skydir, self.bg_vals)
    def evaluate_at(self, skydirs):
        return super(IsotropicConvolver, self).__call__(skydirs, self.bg_vals)
        
    def model_integrator(self):
        return lambda func : self.band.exposure.model_integral(self.center, func, 
            self.band.emin, self.band.emax)*self(self.center) * self.ap_average
