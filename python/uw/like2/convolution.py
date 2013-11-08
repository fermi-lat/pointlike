"""
Convolution interface for like2
Extends classes from uw.utilities 


$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/convolution.py,v 1.2 2013/10/29 14:09:17 burnett Exp $
author:  Toby Burnett
"""
import os, pickle, zipfile 
import numpy as np
import pandas as pd
from uw.utilities import keyword_options
from uw.utilities import convolution as utilities_convolution

import skymaps #from Science Tools: for SkyDir 

class FillMixin(object):
    """A Mixin class for like2 convolution, to replace functions in utilities.convolution
    """
    def fill(self, skyfun):
        """ Evaluate skyfun along the internal grid and return the resulting array.
        (Identical to superclass, except skyfun can be either a python functor or a 
        C++ SkySkySpectrum)
        """
        v = np.empty(self.npix*self.npix)
        if isinstance(skyfun, skymaps.SkySpectrum):
            skymaps.PythonUtilities.val_grid(v,self.lons,self.lats,self.center,skyfun)
        else:
            def pyskyfun(u):
                return skyfun(skymaps.SkyDir(skymaps.Hep3Vector(u[0],u[1],u[2])))
            skymaps.PythonUtilities.val_grid(v,self.lons,self.lats,self.center,
                skymaps.PySkyFunction(pyskyfun))
        return v.reshape([self.npix,self.npix])
        
    def bg_fill(self, exp, dm, cache=None, ignore_nan=False):
        """ Evaluate product of exposure and diffuse map on the grid
        exp : SkyFunction for exposure
        dm  : [SkyFuntion for diffuse map | None]
            If None, expect predetermined values in cache, which may be an array or a scalar
        """
        #print 'filling with product of exposure "%s" model "%s"' % (exp, dm)
        
        if dm is None:
            assert cache is not None, 'Logic error'
            self.bg_vals = self.fill(exp) * cache
        else:
            def exp_dm(skydir):
                    return exp(skydir)*dm(skydir)
            self.bg_vals = self.fill(exp_dm)
    
        #self.bg_vals = self.fill(exp) * (self.fill(dm) if cache is None else cache) #product of exposure and map
        #self.dm_vals = self.fill(dm) #temporary
        #self.exp_vals = self.fill(exp)
        # check for nans, replace with zeros if not full ROI
        nans = np.isnan(self.bg_vals)
        if np.all(nans):
            raise eException('Diffuse source %s has no overlap with ROi' % dm.filename)
        if np.any(nans) and ignore_nan:
            self.bg_vals[nans]=0
            
    def psf_fill(self, psf):
        """ Evaluate PSF on the grid
        """
        #print 'filling with psf %s' % psf
        psf_vals = psf(self.dists).reshape([self.npix,self.npix])
        self.psf_vals = psf_vals / psf_vals.sum()
        
    def set_npix(self, psf, edge=0, r_multi=1.2, r_max=20):
        """ modify the npix with
            psf : PSF object
            edge: float --Source size (degrees)
            r_multi   float multiple of r95 to set max dimension of grid
            r_max     float an absolute maximum (half)-size of grid (deg)
            """
        r95 = psf.inverse_integral(95)
        rad = r_multi*r95 + edge
        rad = max(min(r_max,rad),edge+2.5)
        npix = int(round(2*rad/self.pixelsize))
        npix += (npix%2 == 0)
        return npix


class ShowMixin(object):
    """ A mixin class to add or replace show methods
    """
    def show_vals(self, vals=None, ax=None, roi_radius=5, roi_dir=None, colorbar=True, npix=None, **kw):
        """Make a display.
        vals : 2-d array of float
            generated by the fill method; expect to be npix x npix
        npix : [int | None]
            if int, override self.npix to for central npix x npix
        """
        import pylab as plt
        if ax is None: fig,ax=plt.subplots()
        if vals is None: vals = self.cvals
        if npix is not None and npix!=self.npix:
            delta = (self.npix-npix)/2
            assert delta>0, 'npix not >= self.npix'
            tvals = vals[delta:delta+npix, delta:delta+npix]
        else: 
            npix=self.npix; tvals = vals
        if roi_radius is not None:
            if roi_dir is None: roi_dir = self.center
            circle = plt.Circle(self.pix(roi_dir),roi_radius/self.pixelsize, color='grey', lw=2,fill=False)
            ax.add_artist(circle)

        v = ax.imshow( tvals.transpose()[::-1],  interpolation='nearest', **kw)
        marker = float(npix)/2
        ax.axvline(marker,color='k')
        ax.axhline(marker,color='k')
        if colorbar: 
            cb = plt.colorbar(v, shrink=0.8)
        def scale(x, factor=1.0):
            return x*factor/self.pixelsize+self.npix/2.
        r = np.arange(-8,9,4)
        ax.set_xticks(scale(r))
        ax.set_xticklabels(map(lambda x:'%.0f'%x ,r))
        ax.set_yticks(scale(r, -1))
        ax.set_yticklabels(map(lambda x:'%.0f'%x ,r))
        return ax.figure
            
    def show(self, roi_radius=None,roi_dir=None, **kwargs):
        """Three subplots: PSF, raw, convolved"""
        import pylab as plt
        from matplotlib.colors import LogNorm
        title = kwargs.pop('title', None)
        if hasattr(self, 'band'):
            roi_radius = self.band.radius
            roi_dir = self.band.sd
        fig, axx = plt.subplots(1,3, figsize=(10,3), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.05)
        if hasattr(self, 'psf_vals'):
            axx[0].imshow(self.psf_vals,interpolation='nearest')
        vmax = self.bg_vals.max()
        norm = LogNorm(vmax=vmax, vmin=vmax/1e3)
        marker = float(self.npix)/2
        for ax,what in zip(axx[1:], (self.bg_vals, self.cvals)  ):
            what[what==0]=vmax/1e6
            ax.imshow(what.transpose()[::-1], norm=norm, interpolation='nearest')
            ax.axvline(marker,color='grey')
            ax.axhline(marker,color='grey')
            if roi_radius is not None:
                if roi_dir is None: roi_dir = self.center
                circle = plt.Circle(self.pix(roi_dir),roi_radius/self.pixelsize, color='grey', lw=2,fill=False)
                ax.add_artist(circle)
        axx[0].set_aspect(1.0)
        if title is not None:
            plt.suptitle(title,fontsize='small')
            
        return fig

class ConvolvableGrid(FillMixin, ShowMixin, utilities_convolution.BackgroundConvolution):
    """ Convolution used by classes below. This subclass uses the mixin classes defined here to:
    
      1) changes the default for a bounds error (to check)
      2) Replaces fill method with version that works for python class
      3) provides useful show methods
      """
    defaults =(
        ('pixelsize', 0.1, 'Size of pixels to use for convolution grid'),
        ('npix',      201, 'Number of pixels (must be an odd number'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, center, **kwargs):
        """ center -- a SkyDir giving the center of the grid on which to convolve bg
            kwargs are passed to Grid.
        """
        keyword_options.process(self, kwargs)
        defaults=dict(bounds_error=False)
        defaults.update(kwargs)
        # note do not use code in superclass needing psf, diffuse function
        super(ConvolvableGrid, self).__init__(center, None, None, **defaults)
        self.center = center
        
    def __repr__(self):
        return '%s.%s: center %s npix %d pixelsize %.2f' %(
            self.__module__,self.__class__.__name__, self.center, self.npix, self.pixelsize)
    def overlap(self): return 1.0
    def exposure_ratio(self): return 1.0

#class Convolver(ConvolvableGrid):
#    pass
#    
#class GridGenerator(object):
#    """Base class for creation of local grid for global diffuse models
#    Note that it applies to a particular subtype (front, back ....)
#    """
#    defaults =(
#        ('pixelsize',0.25,'Pixel size for convolution grid'),
#        ('npix', 61, 'number of pixels for 5 deg'),
#        ('ignore_nan', True, 'replace nan values with zero (generates warning)'), 
#        )
#
#    @keyword_options.decorate(defaults)
#    def __init__(self, band, diffuse_source, **kwargs): #psf, exposure, roi_dir, diffuse_source,  **kwargs):
#        """
#        Factory class to create grids of convolved diffuse
#        
#        parameters
#        ----------
#        band : ROIBand or equivalent
#            Contains position, psf, exposure, selected for a particular event type
#
#        diffuse_source : DiffuseBase object. It must be for the same event type, if relevant
#        """
#        keyword_options.process(self, kwargs)
#        self.roi_dir = band.sd #roi_dir
#        self.psf, self.exposure = band.psf, band.exposure
#        self.band=band
#        self.diffuse_source = diffuse_source
#        self.name = diffuse_source.name
#        
#    def __repr__(self):
#        return '%s for %s' % (self.__class__.__name__, self.diffuse_source.__repr__() )
#        
#    def set_npix(self,  edge=0, r_multi=1.2, r_max=20):
#        """ modify the npix with
#            edge: float --Source size (degrees)
#            r_multi   float multiple of r95 to set max dimension of grid
#            r_max     float an absolute maximum (half)-size of grid (deg)
#            """
#        r95 = self.psf.inverse_integral(95)
#        rad = r_multi*r95 + edge + self.band.radius
#        rad = max(min(r_max,rad),edge+2.5)
#        npix = int(round(2*rad/self.pixelsize))
#        npix += (npix%2 == 0)
#        return npix
#   
#    def _get_convolver(self):
#        """ load a convolution object with npix adjusted for current energy""" 
#        npix = self.npix if self.npix is not None else self.set_npix()
#        return ConvolvableGrid(self.band.sd, npix=npix, pixelsize=self.pixelsize)
#        
#    def __call__(self, energy, tol=0.5):
#        """ return a convolved grid for a band in an ROI
#        parameters
#        ----------
#        energy : float
#            intermediate energy for the band
#        """
#        # First set energy in appropriated places (perhaps Band should know how to do this)
#        dm = self.diffuse_source
#        if hasattr(dm, 'setEnergy'): # case for extended
#            dm.setEnergy(energy)
#        self.exposure.setEnergy(energy)
#        self.psf.setEnergy(energy)
#        
#        # get the convolution class, with perhaps new npix
#        grid = self._get_convolver()
#        
#        grid.bg_fill(self.exposure,dm)
#         
#        #set the PSF map
#        grid.psf_fill(self.psf)
#
#        # finally convolve it unless fairly flat or psf is very local
#        bgmax, bgmin = grid.bg_vals.max(), grid.bg_vals.min()
#        if bgmax==0 or (bgmax-bgmin)/bgmax < tol:
#            grid.cvals = grid.bg_vals
#        else:
#            grid.convolve()
#        return grid
#
#
#class CachedGridGenerator(GridGenerator):
#    def __init__(self, band, center, dfun, **kwargs):
#        self.center=center
#        self.band = band
#        self.__dict__.update(kwargs)
#        self.quiet = True# kwargs['quiet']
#        self.npix=61
#        self.diffuse_source = dfun
#        roi_index = skymaps.Band(12).index(center)
#        try:
#            self.filename = dfun.files[roi_index]
#            self.cached_diffuse = pickle.load(dfun.opener(self.filename))
#        except Exception, msg:
#            raise Exception( 'Diffuse cache file # %d not found:%s' %(roi_index,msg))
#        self.emins = [cd['emin'] for cd in self.cached_diffuse]
#        if hasattr(dfun, 'kw') and dfun.kw is not None: # check for extra keywords from diffuse spec.
#            # Manage keywords found in the 
#            if dfun.kw['correction'] is not None:
#                if not self.quiet:print '\t%s loading corrections from %s.kw:' % (self.__class__.__name__, dfun.__class__.__name__)
#                df = pd.read_csv(dfun.kw['correction'], index_col=0) 
#                self.corr = df.ix['HP12_%04d'%roi_index].values
#                if not self.quiet:print '\tcorrection file: "%s"' % dfun.kw['correction']
#                if not self.quiet:print '\tcorrections: %s' %self.corr.round(3)
#            else: self.corr=None
#            self.systematic = dfun.kw.get('systematic', None)
#            if self.systematic is not None:
#                if not self.quiet:print '\tsystematic: %.3f' % self.systematic
#
#
#    def __call__(self, energy):
#        # find the appropriate cached grid
#        for index in range(len(self.emins)):
#            if energy>self.emins[index] and (index==len(self.emins)-1 or energy<self.emins[index+1])\
#                : break
#        emin = self.emins[index]    
#        assert energy/emin < 1.8 and energy> emin, 'too large a factor: energy, emin=%.0f,%.0f\nemins=%s' % (energy, emin, self.emins)
#        cd = self.cached_diffuse[index]
#        # create a convolvable grid from it
#        energy = cd['energy']
#        if cd['center'].difference(self.band.sd)>0.1:
#            print '%s, WARNING: Logic error? %s not %s' % (self, cd['center'], self.band.sd)
#        grid = ConvolvableGrid(self.center,   npix=cd['npix'], pixelsize=cd['pixelsize'])
#            
#        self.band.set_energy(energy)
#        grid.psf_fill(self.band.psf)
#        
#        # apply correction factor if any; use last value for all higher (depends only on energy)
#        vals = cd['vals']
#        ds = self.diffuse_source
#        if hasattr(self, 'corr') and self.corr is not None:
#            c = self.corr[index] if index<len(self.corr) else self.corr[-1]
#            vals *= c
#        else: c = None
#        
#        # finally do the convolution on the product of exposure and diffuse map, which is passed in 
#        #print 'applying correction %s: exposure at center %.3g average vals %.3g convolving...' % (
#        #   c, self.band.exposure(self.center), vals.mean())
#        grid.bg_fill(self.band.exposure, None, cache=vals)
#        grid.convolve()
#        return grid
#
#
#class CachedGridConvolver(Convolver):
#
#    def __init__(self, source, band, **kwargs):
#        self.center=band.sd
#        super(CachedGridConvolver, self).__init__(self.center, **kwargs)
#        self.band = band
#        self.quiet = True# kwargs['quiet']
#        self.npix=61
#        self.source = source
#        dfun = self.source.dmodel[band.event_type]
#        roi_index = skymaps.Band(12).index(self.center)
#        try:
#            self.filename = dfun.files[roi_index]
#            self.cached_diffuse = pickle.load(dfun.opener(self.filename))
#        except Exception, msg:
#            raise Exception( 'Diffuse cache file # %d not found:%s' %(roi_index,msg))
#        self.emins = [cd['emin'] for cd in self.cached_diffuse]
#        
#        # find the appropriate cached grid for this index and energy
#        energy = band.energy
#        for index in range(len(self.emins)):
#            if energy>self.emins[index] and (index==len(self.emins)-1 or energy<self.emins[index+1])\
#                : break
#        emin = self.emins[index]    
#        assert energy/emin < 1.8 and energy> emin, 'too large a factor: energy, emin=%.0f,%.0f\nemins=%s' % (energy, emin, self.emins)
#        cd = self.cached_diffuse[index]
#        # create a convolvable grid from it
#        energy = cd['energy']
#        if cd['center'].difference(self.band.sd)>0.1:
#            print '%s, WARNING: Logic error? %s not %s' % (self, cd['center'], self.band.sd)
#        # now set up the convolver
#        super(CachedGridConvolver, self).__init__(self.center, npix=cd['npix'], pixelsize=cd['pixelsize'])
#        self.cd = cd 
#   
#    def convolve(self):
#        self.psf_fill(self.band.psf)
#        self.bg_fill(self.band.exposure, None, cache=self.cd['vals'])
#        super(CachedGridConvolver,self).convolve()
#        inside = self.dists< np.radians(self.band.radius)
#        ap_average = self.cvals[inside].mean()
#        return ap_average, 1.0
#        
#    def __call__(self, skydir):
#        return super(CachedGridConvolver, self).__call__(skydir, self.cvals)
# 
#
#class IsotropicConvolver(Convolver):
#    """ avoid convolution, using only the exposure and the value of the isotropic for the roi band
#    """
#
#    defaults =(
#        ('pixelsize', 0.25, 'Size of pixels to use for convolution grid'),
#        ('npix',      61, 'Number of pixels (must be an odd number'),
#        )
#
#    @keyword_options.decorate(defaults)
#    def __init__(self, source, band, **kwargs):
#        keyword_options.process(self, kwargs)
#        super(IsotropicConvolver, self).__init__(band.sd, **kwargs)
#
#        self.band = band
#        self.source = source
#        self.center = band.sd
#        self.dmodel = source.dmodel[band.event_type]
#    
#    def convolve(self, energy=None):
#        if energy is not None: self.band.set_energy(energy)
#        self.bg_fill(self.band.exposure, None, cache=self.dmodel(self.center))
#        inside = self.dists< np.radians(self.band.radius)
#        self.ap_average = self.bg_vals[inside].mean()
#    
#    def __call__(self, skydir):
#        return super(IsotropicConvolver, self).__call__(skydir, self.bg_vals)
#    def evaluate_at(self, skydirs):
#        return super(IsotropicConvolver, self).__call__(skydirs, self.bg_vals)
#        
#    def model_integrator(self):
#        return lambda func : self.band.exposure.model_integral(self.center, func, 
#            self.band.emin, self.band.emax)*self(self.center) * self.ap_average
#
#
#    
#class MapCubeConvolver(Convolver):
#
#    @keyword_options.decorate(ConvolvableGrid.defaults)
#    def __init__(self, source, band, **kwargs):
#        keyword_options.process(self, kwargs)
#        super(MapCubeConvolver, self).__init__(band.sd, **kwargs)
#
#        self.band = band
#        self.source = source
#        self.center = band.sd
#        self.dmodel = source.dmodel[band.event_type]
#    
#    def convolve(self, energy=None):
#        if energy is not None: self.band.set_energy(energy)
#        self.psf_fill(self.band.psf)
#        self.bg_fill(self.band.exposure, self.dmodel)
#        super(MapCubeConvolver,self).convolve()
#        inside = self.dists< np.radians(self.band.radius)
#        ap_average = self.cvals[inside].mean()
#        return ap_average, 1.0
#    
#    def __call__(self, skydir):
#        return super(MapCubeConvolver, self).__call__(skydir, self.bg_vals)
#
#
#class ExtendedConvolver(Convolver):
#    """Manage the convolution of an extended source with the PSF
#    Maintains separate grids corresponding to different energies
#    """
#    defaults =(
#        ('pixelsize', 0.1, 'Size of pixels to use for convolution grid'),
#        ('npix',      201, 'Number of pixels (must be an odd number'),
#        )
#
#    @keyword_options.decorate(defaults)
#    def __init__(self, source, band, **kwargs):
#        """
#        parameters
#        -----------
#        band  : ROIBand object, which has appropriate psf and exposure for a given event type
#        source: An extended source, must have a map cube description
#        """
#        assert hasattr(band, 'energy'), 'Not a band object? %s' %band
#        assert hasattr(source, 'dmodel'), 'Not an extended source object?' % source
#        keyword_options.process(self, kwargs)
#        self.band=band
#        self.source=source
#        self.center = self.source.skydir
#        self.roi_dir =band.sd
#        
#        self.cvals = dict()
#        self.dm_vals = None
#        self.grid = None
#
#    def setup_image(self):
#        # load the SkyImage as a numpy 2-d array (not used yet?)
#        assert hasattr(source.dmodel, 'skyfun'), 'Wrong dmodel? %s' % source.dmodel
#        skyfun = source.dmodel.skyfun
#        naxis1 = skyfun.naxis1()
#        naxis2 = skyfun.naxis2()
#        self.image = np.asarray(skyfun.image()).reshape(naxis1,naxis2)
#        
#    def set_energy(self, energy):
#        self.band.set_energy(energy)
#    
#    @property
#    def energy(self):
#        return round(self.band.energy)
#        
#    def create_grid(self):
#        """create a grid for all convolutions, evaluate exended model on it
#        Note that is is normalized so the integral over solid angle should be 1
#        """
#        self.grid = ConvolvableGrid(self.center, pixelsize=self.pixelsize, npix=self.npix)
#        self.dm_vals = self.grid.fill(self.source.dmodel)
#        self.dm_vals/= (self.dm_vals.sum() * np.radians(self.pixelsize)**2) 
#
#    def overlap_mask(self):
#        """ 
#        return a npix x npix array of bools for the part of the grid inside the ROI circle
#        """
#        x,y = self.grid.pix(self.roi_dir)
#        npix = self.grid.npix
#        dx2 =((np.arange(npix)-x)**2).reshape(npix,1)
#        dy2 =((np.arange(npix)-y)**2)
#        d2  = dx2+dy2 # dx2 *dy2 is expanded to a square matrix
#        #return d2.max(), d2.min()
#        return d2  <= (self.band.radius/self.grid.pixelsize)**2
#
#    def convolve(self, energy=None):
#        """ perform the convolution, creating a grid for the specified energy
#        
#        return the tuple ap_average, overlap, psf_overlap
#        """
#        if energy is not None: self.band.set_energy(energy)
#        if self.grid is None: self.create_grid()
#        
#        # chedk the PSF size: if less than grid spacing, do not convolve, keep original grid
#        #if self.band.psf.inverse_integral(68)< self.pixelsize:
#        #    self.cvals[self.energy]=None
#        # Fill the psf in the grid
#        self.grid.psf_fill(self.band.psf)
#        
#        # now look at values, decide if want to convolve
#        exp_grid = self.grid.fill(self.band.exposure) # this is expensive, 1.6 s for npix=201
#        self.grid.bg_vals = exp_grid * self.dm_vals
#        self.grid.convolve()  
#        
#        # save a copy for this energy
#        self.cvals[self.energy]=self.grid.cvals.copy()
#        
#        # calculations needed which depend on this convolution
#        cvals = self.cvals[self.energy]
#        inside = self.overlap_mask() #self.grid.dists< np.radians(self.band.radius)
#        ap_average = cvals[inside].mean()
#        overlap = cvals[inside].sum() / cvals.sum()
#        pvals = self.grid.psf_vals
#        psf_overlap = pvals[inside].sum() / pvals.sum() 
#        return ap_average, overlap,  psf_overlap
#   
#    def __call__(self, skydir, force=False):
#        """ return value of perhaps convolved grid for the position
#        skydir : SkyDir object | [SkyDir]
#        """
#        # TODO: if not convolved, just evaluate product of exposure and 
#        #if force or self.cvals[self.energy] is None:
#        #    return self.band.exposure(skydir) * self.grid(skydir, self.dm_vals)
#            
#        return self.grid(skydir, self.cvals[self.energy])
#        
#    def __repr__(self):
#        return '%s.%s: \n\tsource: %s\n\tband  : %s\n\tpixelsize: %.1f, npix: %d' % (
#            self.__module__, self.__class__.__name__, self.source,self.band, self.pixelsize, self.npix)
#  
#        
#    def exposure_ratio(self):
#        """ the ratio of the exposure at the source, to that at the center of the ROI)
#        """
#        return self.band.exposure(self.source.skydir)/self.band.exposure(self.band.sd)
#          
#    def show(self, title, vals, ax=None, logscale=False,
#            roi_radius=None, colorbar=True):
#        import pylab as plt
#        from matplotlib.colors import LogNorm
#        if ax is None:
#            fig, ax = plt.subplots(figsize=(4,4))
#        norm = LogNorm(vmin=0.001, vmax=1.0) if logscale else None
#        roi_radius = self.band.radius if roi_radius is None else roi_radius
#        self.grid.show_vals(vals/vals.max(), ax=ax, roi_radius=roi_radius,
#            roi_dir=self.roi_dir,norm=norm, colorbar=colorbar)
#        ax.set_title(title, size=10)
#        return ax.figure
#    
#    def show_psf(self, ax=None, colorbar=True):
#        return self.show('PSF for %.0f MeV' % self.band.energy, 
#            self.grid.psf_vals, roi_radius=0, ax=ax, colorbar=colorbar)
#    
#    def show_cvals(self, ax=None, colorbar=True):
#        return self.show(
#            'convolved %s at %.0f MeV' % (self.source.name,self.band.energy), 
#            self.cvals[self.energy], logscale=True, ax=ax, colorbar=colorbar)
#            
#    def show_source(self, ax=None, colorbar=True):
#        return self.show(
#            'Source %s' % (self.source.name),self.dm_vals, 
#                logscale=True, ax=ax, colorbar=colorbar)
#            
#    def show_all(self):
#        import pylab as plt
#        fig, axx = plt.subplots(1, 3, figsize=(12,4))
#        self.show_psf(ax=axx[0], colorbar=False)
#        self.show_source(ax=axx[1], colorbar=False)
#        self.show_cvals(ax=axx[2], colorbar=False)
#        return fig