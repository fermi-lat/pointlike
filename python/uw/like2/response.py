"""
Classes to compute response from various sources
 
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/response.py,v 1.7 2013/11/18 00:03:24 burnett Exp $
author:  Toby Burnett
"""
import os, pickle
import numpy as np
import skymaps
from uw.utilities import keyword_options
from . import convolution

class ResponseException(Exception): pass

diffuse_grid_defaults = (
        ('pixelsize', 0.25, 'Size of pixels to use for convolution grid'),
        ('npix',      61,   'Number of pixels (must be an odd number'),
        ('quiet',     False, ''),
        )
        
extended_grid_defaults = [
        ['pixelsize', 0.2, 'Size of pixels to use for convolution grid'],
        ('size',      14,  'Size of grid; npix set to size/pixelsize'),
        ('quiet',     False, ''),
        ]

class Response(object):
    """ base class for classes that manage the response of a source, in total counts 
    or count density for any position within the ROI.   
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
        self.quiet = kwargs.pop('quiet', True)
        self.initialize()
        
    def update(self):
        self.evaluate()
        
    def __repr__(self):
        return '%s.%s: ROI at %s source "%s" at %s' %( self.__module__,self.__class__.__name__,
            self.roicenter, self.source.name, self.source.skydir)

    def exposure_integral(self):
        """Integral of the exposure times the flux at the given position"""
        return self.band.integrator(self.source.model)
    
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
        model = self.spectral_model
        self.expected = self.band.integrator(model)
        assert not np.isinf(self.expected), 'model integration failure'
        self.counts =  self.expected * self.overlap
        self.model_grad = self.band.integrator( model.gradient)[model.free] #* self.exposure_ratio
        if self.band.has_pixels:
            self.pix_counts = self.pixel_values * self.expected
        
    def grad(self, weights, exposure_factor=1): 
        """ contribution to the overall gradient
        weights : arrary of float
            
            weights = self.data / self.model_pixels
        Assume that evaluate has set model_grad
        """
        model = self.spectral_model
        if np.sum(model.free)==0 : return []
        # Calculate the gradient of a spectral model (wrt its parameters) integrated over the exposure.
        #g = self.band.integrator( model.gradient)[model.free] #* self.exposure_ratio
        g = self.model_grad
        apterm = exposure_factor* self.overlap
        pixterm = (weights*self.pixel_values).sum() if self.band.has_pixels else 0
        return g * (apterm - pixterm)

    def __call__(self, skydir):
        return self.band.psf(skydir.difference(self.source.skydir))  * self.expected
     

    
class DiffuseResponse(Response):
        
    defaults = diffuse_grid_defaults
    @keyword_options.decorate(defaults)
    def __init__(self, source, band, **kwargs):
        keyword_options.process(self, kwargs)
        self.setup=False
        self.overlap=1
        super(DiffuseResponse, self).__init__(source, band, **kwargs)
        
    def initialize(self):
        if self.setup: return
        self.setup=True
        #set up the spatial model 
        self.dmodel = self.source.dmodel[self.band.event_type]
        self.dmodel.load()
        self.energy = self.band.energy
        self.dmodel.setEnergy(self.band.energy)
        
        self.create_grid()
        grid = self.grid
        inside = grid.dists< self.band.radius_in_rad
        self.ap_average = grid.cvals[inside].mean()
        self.delta_e = self.band.emax - self.band.emin
        self.factor = self.ap_average * self.band.solid_angle * self.delta_e
        if self.band.has_pixels:
            self.pixel_values = grid(self.band.wsdl, grid.cvals) * self.band.pixel_area * self.delta_e

        self.evaluate()
        
    def create_grid(self):
        # create a grid for evaluating counts integral over ROI, individual pixel predictions
        grid = self.grid= convolution.ConvolvableGrid(center=self.roicenter, 
                npix=self.npix, pixelsize=self.pixelsize)
        # this may be overridden
        self.fill_grid()

            
    def fill_grid(self):
        # fill and convolve the grid
        self.grid.psf_fill(self.band.psf)
        self.grid.bg_fill(self.band.exposure, self.dmodel)
        self.grid.convolve()
        
    def evaluate(self):
        norm = self.source.model(self.band.energy)
        self.counts = norm * self.factor
        if self.band.has_pixels:
            self.pix_counts = self.pixel_values * norm

    def grad(self, weights, exposure_factor=1): 
        """ contribution to the overall gradient
        weights : arrary of float
            weights = self.data / self.model_pixels
        """
        model = self.spectral_model
        if np.sum(model.free)==0 : 
            return []
        pixterm = ( self.pixel_values * weights ).sum() if self.band.has_pixels else 0
        return (self.factor*exposure_factor - pixterm) * model.gradient(self.energy)[model.free] 
        
    def __call__(self, skydir):
        return self.grid(skydir, self.grid.cvals)[0] * self.delta_e
    
    def _keyword_check(self, roi_index):
        # check for extra keywords from diffuse spec.
        import pandas as pd
        if hasattr(self.source,'corr'): 
            self.corr = self.source.corr
            return
        dfun = self.dmodel
        
        if hasattr(dfun, 'kw') and len(dfun.kw.keys())>1: 
            # Manage keywords found in the 
            if dfun.kw['correction'] is not None:
                if not self.quiet:print '\t%s loading corrections for source %s from %s.kw:' \
                    % (self.__class__.__name__, self.source.name,  dfun.__class__.__name__)
                print '\t%s loading corrections for source %s from %s.kw:' \
                    % (self.__class__.__name__, self.source.name,  dfun.__class__.__name__)
                corr_file = os.path.expandvars(dfun.kw['correction'])
                if not os.path.exists(corr_file):
                    corr_file = os.path.expandvars(os.path.join('$FERMI','diffuse', corr_file))
                try:
                    df = pd.read_csv(corr_file, index_col=0) 
                except Exception, msg:
                    raise Exception('Error loading correction file %s: %s'% (corr_file,msg))
                self.corr = df.ix['HP12_%04d'%roi_index].values
                if not self.quiet:print '\tcorrection file: "%s"' % corr_file
                if not self.quiet:print '\tcorrections: %s' %self.corr.round(3)
            else: self.corr=None
            self.systematic = dfun.kw.get('systematic', None)
            if self.systematic is not None:
                if not self.quiet:print '\tsystematic: %.3f' % self.systematic
        else: 
            self.corr =self.systematic = None
            
        #cache this result in the diffuse object
        self.source.corr = self.corr 
        self.source.systematic = self.systematic

        
class CachedDiffuseResponse(DiffuseResponse):
        
    def create_grid(self):
        """ set up the grid from the cached files """
        
        roi_index = skymaps.Band(12).index(self.roicenter)
        dmodel = self.dmodel
        try:
            self.filename = dmodel.files[roi_index]
            self.cached_diffuse = pickle.load(dmodel.opener(self.filename))
        except Exception, msg:
            raise ResponseException( 'Diffuse cache file # %d not found:%s' %(roi_index,msg))
        self.emins = [cd['emin'] for cd in self.cached_diffuse]
        
        self._keyword_check(roi_index)
        
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

class IsotropicResponse(DiffuseResponse):

    defaults = diffuse_grid_defaults
    @keyword_options.decorate(defaults)    
    def __init__(self, source, band, **kwargs):
        keyword_options.process(self, kwargs)
        super(IsotropicResponse, self).__init__(source, band, **kwargs)
        
    def evaluate(self):
        # deal with FrontBackConstant case, which uses different conatants for front/back
        if hasattr(self.source.model, 'ct'): # bit ugly
            self.source.model.ct=self.band.event_type
        super(IsotropicResponse, self).evaluate()
            
    def fill_grid(self):
        # fill the grid for evaluating counts integral over ROI, individual pixel predictions
        grid = self.grid
        #grid.bg_fill(self.band.exposure, None, cache=self.dmodel(band.skydir))
        grid.cvals = grid.fill(self.band.exposure) * self.dmodel(self.band.skydir)
        #grid.cvals = grid.bg_vals
        

class ExtendedResponse(DiffuseResponse):
    
    def __init__(self, source, band, **kwargs):
        defaults = dict( [x[:2] for x in extended_grid_defaults])
        defaults.update(kwargs)
        if 'npix' not in defaults:
            npix = defaults.pop('size')/defaults['pixelsize']
            defaults['npix'] = int(npix) | 1 # make odd
            #print 'setting npix', defaults
        self.quiet=defaults.get('quiet', True)
        self.initialized = False
        super(ExtendedResponse, self).__init__(source, band, **defaults)
            
    def exposure_integral(self):
        """Perform integral of the exposure times the flux at the source position
        (Override base to specify source position)
        """
        return self.band.exposure.integrator(self.source.skydir, 
                self.band.emin, self.band.emax)(self.source.model)

      
    def initialize(self):
        #set up the spatial model NOTE THIS NEEDS TO BE SAVED
        self.dmodel = self.source.dmodel[self.band.event_type]
  
        self.center = self.source.skydir
         
        self.cvals = dict()
        self.dm_vals = None
        self.create_grid()
        self.convolve()
        self.evaluate()
        self.initialized = True

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
        Have the source own the grid; each band has 
        """
        if not hasattr(self.source, 'grid'):
            self.source.grid = convolution.ConvolvableGrid(self.center, pixelsize=self.pixelsize, npix=self.npix)
        self.grid =self.source.grid
        self.dm_vals = self.grid.fill(self.source.dmodel)
        self.dm_vals/= (self.dm_vals.sum() * np.radians(self.pixelsize)**2) 

    def overlap_mask(self):
        """ 
        return a npix x npix array of bools for the part of the grid, which is centered on the source,
        which is inside the ROI circle
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
        assert not self.initialized, 'Already initialized?'
        if energy is not None: self.band.set_energy(energy)
        #if self.grid is None: self.create_grid()
        
        # chedk the PSF size: if less than grid spacing, do not convolve, keep original grid
        #if self.band.psf.inverse_integral(68)< self.pixelsize:
        #    self.cvals[self.energy]=None
        # Fill the psf in the grid
        self.grid.psf_fill(self.band.psf)
        
        # now look at values, decide if want to convolve    
        exp_grid = self.grid.fill(self.band.exposure) # this is expensive, 1.6 s for npix=201
        self.grid.bg_vals = exp_grid * self.dm_vals
        self.grid.convolve()  
        
        # save a copy of cvals for this band
        self.cvals = self.grid.cvals.copy()
        
        # calculations needed which depend on this convolution
        cvals = self.cvals
        inside = self.overlap_mask() #self.grid.dists< np.radians(self.band.radius)
        self.ap_average = cvals[inside].mean()
        self.ap_center = self.grid(self.center, cvals)
        self.overlap = cvals[inside].sum() / cvals.sum()
        pvals = self.grid.psf_vals
        self.psf_overlap = pvals[inside].sum() / pvals.sum() 
        self.exposure_at_center = self.band.exposure(self.source.skydir)
        self.exposure_ratio = self.exposure_at_center / self.band.exposure(self.roicenter)
        self.factor = self.overlap * self.exposure_ratio 
        if self.band.has_pixels:
            self.pixel_values = self.grid(self.band.wsdl, self.cvals)\
                /self.exposure_at_center * self.band.pixel_area
        self.expint = self.exposure_integral() # save initial value
        self.initcounts = self.expint * self.factor # "
        self.initmodel = self.source.model.copy()
        
    def evaluate(self):
        total_counts = self.exposure_integral()
        self.counts = total_counts * self.factor
        if self.band.has_pixels:
            self.pix_counts = self.pixel_values * total_counts
        
    def grad(self, weights, exposure_factor=1): 
        """ contribution to the overall gradient
        weights : arrary of float
            weights = self.data / self.model_pixels
        """
        model = self.spectral_model
        if np.sum(model.free)==0 : 
            return []
        #pixterm = ( self.pixel_values * weights ).sum() * self.exposure_at_center if self.band.has_pixels else 0
        #return (self.factor*exposure_factor - pixterm) * model.gradient(self.energy)[model.free] 
        g = self.band.integrator( model.gradient)[model.free] #* self.exposure_ratio
        apterm = exposure_factor * self.factor #self.overlap
        pixterm = (weights*self.pixel_values).sum() if self.band.has_pixels else 0
        return g * (apterm - pixterm)

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
            roi_dir=self.roicenter, norm=norm, colorbar=colorbar)
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
