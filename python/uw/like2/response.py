"""
Classes to compute response from various sources
 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/response.py,v 1.2 2013/11/08 04:30:05 burnett Exp $
author:  Toby Burnett
"""
import os, pickle
import numpy as np
import skymaps
from uw.utilities import keyword_options
from . import convolution

class ResponseException(Exception): pass

class Response(object):
    """ base class for classes that manage the response of a source, in total counts or count density for any position within the ROI    
    """
    def __init__(self, source, band, **kwargs):
        """
        source : Source object, inherit from sources.Source
            skydir : position of source, or None if global
            model : associated specral model
        band : ROIband object
            psf, exposure functions for the event type
            skydir, radius : location, size of ROI
            emin, emax : energy range
            ----- pixelization, data ---
            wsdl : list of pixel positions, perhaps None
            pixel_size : pixel solid angle
        """
        self.band=band
        self.source=source
        self.roicenter = self.band.skydir
        self.initialize()
        
    def __repr__(self):
        return '%s.%s: ROI at %s source "%s" at %s' %( self.__module__,self.__class__.__name__,
            self.roicenter, self.source.name, self.source.skydir)

    def exposure_integral(self, skydir=None):
        if skydir is None: skydir = self.roicenter
        return self.band.exposure.model_integral(skydir, self.source.model, self.band.emin, self.band.emax)
    def __call__(self, skydir):
        """return the counts/sr for the source at the position"""
        raise NotImplemented
    @property
    def spectral_model(self):
        return self.source.model

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
        if self.band.has_pixels:
            wsdl = self.band.wsdl
            rvals  = np.empty(len(wsdl),dtype=float)
            self.band.psf.cpsf.wsdl_val(rvals, self.source.skydir, wsdl) #from C++: sets rvals
            self.pixel_values = rvals * self.band.pixel_area
        self.evaluate()
        
    def evaluate(self): #, weights=None, exposure_factor=1):
        """ update values of counts, pix_counts used for likelihood calculation, derivatives
        Called when source parameters change
        """
        expected = self.exposure_integral(self.source.skydir)
        self.counts =  expected * self.overlap
        if self.band.has_pixels:
            self.pix_counts = self.pixel_values* expected
        
    def __call__(self, skydir):
        return self.band.psf(skydir.difference(self.source.skydir)) \
            * self.exposure_integral(self.source.skydir)
     
diffuse_grid_defaults = (
        ('pixelsize', 0.25, 'Size of pixels to use for convolution grid'),
        ('npix',      61,   'Number of pixels (must be an odd number'),
        )
        
extended_grid_defaults = (
        ('pixelsize', 0.1, 'Size of pixels to use for convolution grid'),
        ('npix',      201,   'Number of pixels (must be an odd number'),
        )

class IsotropicResponse(Response):

    defaults = diffuse_grid_defaults
    @keyword_options.decorate(defaults)    
    def __init__(self, source, band, **kwargs):
        keyword_options.process(self, kwargs)

        super(IsotropicResponse, self).__init__(source, band, **kwargs)
        
    def initialize(self):

        #set up the diffuse model
        self.dmodel = self.source.dmodel[self.band.event_type]
        self.dmodel.setEnergy(self.band.energy)
        
        # create a temporary grid for evaluating counts integral over ROI, individual pixel predictions
        grid = self.grid= convolution.ConvolvableGrid(center=self.roicenter, 
                npix=self.npix, pixelsize=self.pixelsize)
        self.flux = self.dmodel(self.roicenter)
        self.center_exposure = self.band.exposure(self.roicenter) # flux at roicenter
        grid.bg_fill(self.band.exposure, None, cache=1)
        inside = grid.dists< np.radians(self.band.radius)
        self.mean_exposure = grid.bg_vals[inside].mean()
        self.factor = self.mean_exposure/self.center_exposure * self.flux * self.band.solid_angle
        self.evaluate()
        
    def evaluate(self):
        self.counts = self.exposure_integral() * self.factor
    
    def __call__(self, skydir):
        return self.exposure_integral(skydir) * self.flux

    def evaluate_at(self, skydirs):
        return super(IsotropicConvolver, self).__call__(skydirs, self.bg_vals)
    def convolve(self, energy=None):
        if energy is not None: self.band.set_energy(energy)
        self.bg_fill(self.band.exposure, None, cache=self.dmodel(self.center))
        inside = self.dists< np.radians(self.band.radius)
        self.ap_average = self.bg_vals[inside].mean()
    
class DiffuseResponse(Response):
        
    defaults = diffuse_grid_defaults
    @keyword_options.decorate(defaults)
    def __init__(self, source, band, **kwargs):
        keyword_options.process(self, kwargs)
        super(DiffuseResponse, self).__init__(source, band, **kwargs)
        
    def initialize(self):
        #set up the spatial model NOTE THIS NEEDS TO BE SAVED
        self.dmodel = self.source.dmodel[self.band.event_type]
        self.dmodel.load()
        self.energy = self.band.energy
        self.setup_grid()
        self.evaluate()
        
    def setup_grid(self):
        self.dmodel.setEnergy(self.band.energy)
        
        # create a temporary grid for evaluating counts integral over ROI, individual pixel predictions
        grid = self.grid= convolution.ConvolvableGrid(center=self.roicenter, 
                npix=self.npix, pixelsize=self.pixelsize)
                
        grid.psf_fill(self.band.psf)
        grid.bg_fill(self.band.exposure, self.dmodel)
        grid.convolve()
        inside = grid.dists< self.band.radius_in_rad
        self.ap_average = grid.cvals[inside].mean()
        self.delta_e = self.band.emax - self.band.emin
        self.factor = self.ap_average * self.band.solid_angle * self.delta_e


    def evaluate(self):
        self.counts = self.source.model(self.band.energy) * self.factor
        
    def __call__(self, skydir):
        return self.grid(skydir, self.grid.cvals)[0] * self.delta_e

        
class CachedDiffuseResponse(DiffuseResponse):

        
    def setup_grid(self):
        """ set up the grid from the cached files """
        
        roi_index = skymaps.Band(12).index(self.roicenter)
        dfun = dmodel = self.dmodel
        try:
            self.filename = dmodel.files[roi_index]
            self.cached_diffuse = pickle.load(dmodel.opener(self.filename))
        except Exception, msg:
            raise ResponseException( 'Diffuse cache file # %d not found:%s' %(roi_index,msg))
        self.emins = [cd['emin'] for cd in self.cached_diffuse]
        if hasattr(dfun, 'kw') and len(dfun.kw.keys())>1: # check for extra keywords from diffuse spec.
            # Manage keywords found in the 
            if dfun.kw['correction'] is not None:
                if not self.quiet:print '\t%s loading corrections from %s.kw:' % (self.__class__.__name__, dfun.__class__.__name__)
                df = pd.read_csv(dfun.kw['correction'], index_col=0) 
                self.corr = df.ix['HP12_%04d'%roi_index].values
                if not self.quiet:print '\tcorrection file: "%s"' % dfun.kw['correction']
                if not self.quiet:print '\tcorrections: %s' %self.corr.round(3)
            else: self.corr=None
            self.systematic = dfun.kw.get('systematic', None)
            if self.systematic is not None:
                if not self.quiet:print '\tsystematic: %.3f' % self.systematic

        
        # find the appropriate cached grid
        energy = self.band.energy
        for index in range(len(self.emins)):
            if energy>self.emins[index] and (index==len(self.emins)-1 or energy<self.emins[index+1])\
                : break
        emin = self.emins[index]    
        assert energy/emin < 1.8 and energy> emin, 'too large a factor: energy, emin=%.0f,%.0f\nemins=%s' % (energy, emin, self.emins)
        cd = self.cached_diffuse[index]
        # create a convolvable grid from it
        self.energy = energy = cd['energy']
        if cd['center'].difference(self.band.sd)>0.1:
            print '%s, WARNING: Logic error? %s not %s' % (self, cd['center'], self.band.sd)
        self.grid = grid = convolution.ConvolvableGrid(self.roicenter,   npix=cd['npix'], pixelsize=cd['pixelsize'])
            
        self.band.set_energy(energy)
        grid.psf_fill(self.band.psf)
        
        # apply correction factor if any; use last value for all higher (depends only on energy)
        vals = cd['vals']
        if hasattr(self, 'corr') and self.corr is not None:
            c = self.corr[index] if index<len(self.corr) else self.corr[-1]
            vals *= c
        else: c = None
        
        # finally do the convolution on the product of exposure and diffuse map, which is passed in 
        grid.bg_fill(self.band.exposure, None, cache=vals)
        grid.convolve()
        inside = grid.dists< np.radians(self.band.radius)
        self.ap_average = grid.cvals[inside].mean()
        self.ap_center = grid(grid.center, grid.cvals)
        
        # a few things needed by evaluate
        self.delta_e = self.band.emax - self.band.emin
        self.factor = self.ap_average * self.band.solid_angle * self.delta_e


class ExtendedResponse(DiffuseResponse):
    defaults = extended_grid_defaults
    
    @keyword_options.decorate(defaults)
    def __init__(self, source, band, **kwargs):
        keyword_options.process(self, kwargs)
        super(ExtendedResponse, self).__init__(source, band, **kwargs)
      
    def initialize(self):
        #set up the spatial model NOTE THIS NEEDS TO BE SAVED
        self.dmodel = self.source.dmodel[self.band.event_type]
  
        self.center = self.source.skydir
         
        self.cvals = dict()
        self.dm_vals = None
        self.grid = None
        self.convolve()
        self.evaluate()

    def setup_image(self):
        # load the SkyImage as a numpy 2-d array (not used yet?)
        assert hasattr(source.dmodel, 'skyfun'), 'Wrong dmodel? %s' % source.dmodel
        skyfun = source.dmodel.skyfun
        naxis1 = skyfun.naxis1()
        naxis2 = skyfun.naxis2()
        self.image = np.asarray(skyfun.image()).reshape(naxis1,naxis2)
        
    def set_energy(self, energy):
        self.band.set_energy(energy)
    
    @property
    def energy(self):
        return round(self.band.energy)
        
    def create_grid(self):
        """create a grid for all convolutions, evaluate exended model on it
        Note that is is normalized so the integral over solid angle should be 1
        """
        self.grid = convolution.ConvolvableGrid(self.center, pixelsize=self.pixelsize, npix=self.npix)
        self.dm_vals = self.grid.fill(self.source.dmodel)
        self.dm_vals/= (self.dm_vals.sum() * np.radians(self.pixelsize)**2) 

    def overlap_mask(self):
        """ 
        return a npix x npix array of bools for the part of the grid inside the ROI circle
        """
        x,y = self.grid.pix(self.roicenter)
        npix = self.grid.npix
        dx2 =((np.arange(npix)-x)**2).reshape(npix,1)
        dy2 =((np.arange(npix)-y)**2)
        d2  = dx2+dy2 # dx2 *dy2 is expanded to a square matrix
        #return d2.max(), d2.min()
        return d2  <= (self.band.radius/self.grid.pixelsize)**2

    def convolve(self, energy=None):
        """ perform the convolution, creating a grid for the specified energy
        
        return the tuple ap_average, overlap, psf_overlap
        """
        if energy is not None: self.band.set_energy(energy)
        if self.grid is None: self.create_grid()
        
        # chedk the PSF size: if less than grid spacing, do not convolve, keep original grid
        #if self.band.psf.inverse_integral(68)< self.pixelsize:
        #    self.cvals[self.energy]=None
        # Fill the psf in the grid
        self.grid.psf_fill(self.band.psf)
        
        # now look at values, decide if want to convolve
        exp_grid = self.grid.fill(self.band.exposure) # this is expensive, 1.6 s for npix=201
        self.grid.bg_vals = exp_grid * self.dm_vals
        self.grid.convolve()  
        
        # save a copy for this energy
        self.cvals[self.energy] = self.grid.cvals.copy()
        
        # calculations needed which depend on this convolution
        cvals = self.cvals[self.energy]
        inside = self.overlap_mask() #self.grid.dists< np.radians(self.band.radius)
        self.ap_average = cvals[inside].mean()
        self.ap_center = self.grid(self.center, cvals)
        self.overlap = cvals[inside].sum() / cvals.sum()
        pvals = self.grid.psf_vals
        self.psf_overlap = pvals[inside].sum() / pvals.sum() 
        self.exposure_ratio = self.band.exposure(self.source.skydir)/self.band.exposure(self.roicenter)
        self.factor = self.overlap * self.exposure_ratio 
   
    def evaluate(self):
        self.counts = self.exposure_integral() * self.factor
        
    def __call__(self, skydir, force=False):
        """ return value of perhaps convolved grid for the position
        skydir : SkyDir object | [SkyDir]
        """
        # TODO: if not convolved, just evaluate product of exposure and 
        #if force or self.cvals[self.energy] is None:
        #    return self.band.exposure(skydir) * self.grid(skydir, self.dm_vals)
            
        return self.grid(skydir, self.cvals[self.energy])
        
    def __repr__(self):
        return '%s.%s: \n\tsource: %s\n\tband  : %s\n\tpixelsize: %.1f, npix: %d' % (
            self.__module__, self.__class__.__name__, self.source,self.band, self.pixelsize, self.npix)
  
          
    def show(self, title, vals, ax=None, logscale=False,
            roi_radius=None, colorbar=True):
        import pylab as plt
        from matplotlib.colors import LogNorm
        if ax is None:
            fig, ax = plt.subplots(figsize=(4,4))
        norm = LogNorm(vmin=0.001, vmax=1.0) if logscale else None
        roi_radius = self.band.radius if roi_radius is None else roi_radius
        self.grid.show_vals(vals/vals.max(), ax=ax, roi_radius=roi_radius,
            roi_dir=self.roi_dir,norm=norm, colorbar=colorbar)
        ax.set_title(title, size=10)
        return ax.figure
    
    def show_psf(self, ax=None, colorbar=True):
        return self.show('PSF for %.0f MeV' % self.band.energy, 
            self.grid.psf_vals, roi_radius=0, ax=ax, colorbar=colorbar)
    
    def show_cvals(self, ax=None, colorbar=True):
        return self.show(
            'convolved %s at %.0f MeV' % (self.source.name,self.band.energy), 
            self.cvals[self.energy], logscale=True, ax=ax, colorbar=colorbar)
            
    def show_source(self, ax=None, colorbar=True):
        return self.show(
            'Source %s' % (self.source.name),self.dm_vals, 
                logscale=True, ax=ax, colorbar=colorbar)
            
    def show_all(self):
        import pylab as plt
        fig, axx = plt.subplots(1, 3, figsize=(12,4))
        self.show_psf(ax=axx[0], colorbar=False)
        self.show_source(ax=axx[1], colorbar=False)
        self.show_cvals(ax=axx[2], colorbar=False)
        return fig    
