"""
Provides classes to encapsulate and manipulate diffuse sources.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/diffuse.py,v 1.11 2012/01/11 14:10:22 burnett Exp $

author: Matthew Kerr, Toby Burnett
"""
import sys, os, types, pickle, glob, copy
import numpy as np
from . import models
from uw.utilities import keyword_options, convolution
import skymaps #from Science Tools
from scipy import optimize # for fsolve

class DiffuseException(Exception):pass

class Simpson(object):
    """ implement quick simpson integral with log scale """
    ## not used in this module yet.
    def __init__(self, emin, emax, nsimps=4):

        # use a higher nsimps at low energy where effective area is jagged
        ns= (2 if emin<200 else 1)*nsimps
        if ns > 0:
            self.points = sp = np.logspace(np.log10(emin),np.log10(emax),ns+1)
            self.vector = sp * (np.log(sp[-1]/sp[0])/(3.*ns)) * \
                                     np.asarray([1.] + ([4.,2.]*(ns/2))[:-1] + [1.])
        else:
            self.points = sp = np.asarray([np.sqrt(emin*emax)])
            self.vector = sp * np.log(emax/emin)
        
        self.delta = emax-emin
            
    def __call__(self, fn):
        vals = np.array([fn(x) for x in self.points]) # this might be inefficient
        return np.sum( self.vector * vals)
 
    def _mean_energy(self, roi_dir, dmodel, exp):
        fcn = lambda e : dmodel.value(roi_dir, e) * exp(roi_dir,e)
        w = self(fcn)/self.delta
        print w, fcn(self.points[0]), fcn(self.points[-1]) 
        return optimize.brentq(lambda x : fcn(x)-w, self.points[0], self.points[-1]) 
       

###=========================================================================###
class DiffuseModel(object):
    """ Base class to implement diffuse angular distributions for an ROI
    """

    defaults = (
        ('quiet',   False,     'Set True to quiet'),
        ('efactor', 10**0.125, 'Energy factor for bin center'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, psf, exposure_manager, roi_dir, diffuse_source, 
                **kwargs):
        """
        exposure_manager : ExposureManager object
            provides exposure for front and back
        psf: object
        
        roi_dir : SkyDir object
            center of the ROI 
        diffuse_source :
        """
        keyword_options.process(self, kwargs)
        self.roi_dir = roi_dir
        self.skydir = None # (flag that this is global)
        self.psf, self.exposure = psf, exposure_manager.exposure
        #self.exposure_correction = exposure.correction
        self.diffuse_source = diffuse_source
        self.name = diffuse_source.name
        self.setup()
        
    def setup(self): pass
    
    def make_grid(self, emin, conversion_type): pass
    
    @property
    def spectral_model(self):
        return self.diffuse_source.smodel


class ConvolvableGrid(convolution.BackgroundConvolution):
    """ subclass of the basic convolution used by all classes below.
    It
      1) changes the default for a bounds error (to check)
      2) provides useful show method
      """
    def __init__(self, *args, **kwargs):
        defaults=dict(bounds_error=False)
        defaults.update(kwargs)
        super(ConvolvableGrid, self).__init__(*args, **defaults)
        
    def show(self, fignum=None, **kwargs):
        import pylab as plt
        title = kwargs.pop('title', None)
        if fignum is not None:
            plt.close(fignum)
            fig = plt.figure(fignum, figsize=(6,3))
        super(ConvolvableGrid, self).show(**kwargs)
        if title is not None:
            plt.suptitle(title,fontsize='small')
            
    
###====================================================================================================###
class DiffuseModelFromCache(DiffuseModel):
    """ define the diffuse model by reading a grid of the flux values, then
        applying exposure and convolution
        """
    defaults = DiffuseModel.defaults 
    
    @keyword_options.decorate(defaults)
    def __init__(self, *args, **kwargs):
        super(DiffuseModelFromCache, self).__init__(*args, **kwargs)
        keyword_options.process(self, kwargs)
    
    def setup(self):
        filename = self.diffuse_source.dmodel[0].filename
        cache_path = os.path.splitext(filename)[0]
        assert os.path.exists(filename), 'oops, %s not found' %filename
        if not os.path.exists(cache_path):
            raise DiffuseException('cache folder %s not found' %cache_path)

        test = 'HP12_%04d*.pickle'% self.diffuse_source.index
        globbed = glob.glob(os.path.join(cache_path, test))
        if len(globbed)!=1:
            raise DiffuseException( 'Diffuse cache file %s not found' %test)
        self.filename= globbed[0]
        self.cached_diffuse = pickle.load(open(self.filename))
        if not self.quiet: print 'Using cached diffuse in %s'%self.filename
        self.emins = [cd['emin'] for cd in self.cached_diffuse]

        
    def make_grid(self, energy, conversion_type):
        """ return a convovled grid
        
        parameters
        ----------
        energy : float
            intermediate energy for the band
        conversion_type : int
            0 or 1 for front or back
        
        This needs some more care to deal with 8 bands/decade
        """
        for index in range(len(self.emins)):
            if energy>self.emins[index] and (index==len(self.emins)-1 or energy<self.emins[index+1])\
                : break
        emin = self.emins[index]    
        assert energy/emin < 1.8 and energy> emin, 'too large a factor: energy, emin=%.0f,%.0f\nemins=%s' % (energy, emin, self.emins)
        #try:
        #    index = self.emins.index(emin)
        #except:
        #    raise DiffuseException(
        #        'Specified emin, %.1f, not in list of cached values'%emin)
        cd = self.cached_diffuse[index]
        
        energy = cd['energy']
        grid = ConvolvableGrid(cd['center'], None, self.psf, 
            npix=cd['npix'], pixelsize=cd['pixelsize'])
            
        # determine the exposure for this energy and conversion type
        exp = self.exposure[conversion_type]
        exp.setEnergy(energy)
        expgrid = grid.fill(exp) 

        # correction factor
        #expgrid *= self.exposure_correction[conversion_type](energy)
        
        # finally to the convolution on the product of exposure and diffuse map, 
        #  using the appropriate PSF
        grid.do_convolution(energy, conversion_type, override_vals=expgrid*cd['vals'] )
        return grid

    
    def show(self, iband=0):
        for ct, title in ((0,'front'),(1,'back')):
            grid = self.make_grid(self.emins[iband],ct)
            for t in (grid.bg_vals, grid.cvals):
                print '%.3e, %.3f' % (t.mean(), t.std()/t.mean())
                grid.show(fignum=ct+1, 
                    title='Emin=%.0f, %s'%(self.emins[iband], title), origin='upper')
  
  
class IsotropicModel(DiffuseModel):
    """ processing appropriate for Isotropic: no convolution needed
    
    """

    def __init__(self, *pars, **kwargs):
        super(IsotropicModel,self).__init__(*pars, **kwargs)
        
    def setup(self):
        if not self.quiet:
            print 'set up isotropic_model(s):'
            for dm in self.diffuse_source.dmodel:
                print '\t%s'% dm.name()
                
            
    def make_grid(self, energy, conversion_type, npix=61, pixelsize=0.25):
        dmodels = self.diffuse_source.dmodel
        dm = dmodels[conversion_type if len(dmodels)>1 else 0]
        dm.setEnergy(energy)
        exp = self.exposure[conversion_type]
        exp.setEnergy(energy)
        grid = ConvolvableGrid(self.roi_dir, None, self.psf, 
            npix=npix, pixelsize=pixelsize)
        
        grid.cvals = grid.fill(exp) * dm(self.roi_dir, energy) 
        assert not np.any(np.isnan(grid.cvals)), \
            'Grid for %s has nan values' %dm.filename
        return grid

class DiffuseModelFromFits( DiffuseModel):
    """ load pattern from a FITS file """
    def __init__(self, *pars, **kwargs):
        super(DiffuseModelFromFits,self).__init__(*pars, **kwargs)
        
    def setup(self):
        if not self.quiet:
            print 'set up diffuse_model(s):'
            for dm in self.diffuse_source.dmodel:
                print '\t%s'% dm.name()
                
    def make_grid(self, energy, conversion_type, npix=61, pixelsize=0.25):
        dmodels = self.diffuse_source.dmodel
        dm = dmodels[conversion_type if len(dmodels)>1 else 0]
        dm.setEnergy(energy)
        exp = self.exposure[conversion_type]
        exp.setEnergy(energy)

        grid = ConvolvableGrid(self.roi_dir, None, self.psf, 
            npix=npix, pixelsize=pixelsize)
        
        grid.cvals = grid.fill(exp) * grid.fill(dm) #product of exposure and map
        # check for nans, replace with zeros if not full ROI
        nans = np.isnan(grid.cvals)
        if np.all(nans):
            raise DiffuseException('Diffuse cube %s has no overlap with ROi' % dm.filename)
        if np.any(nans):
            grid.cvals[nans]=0

        return grid
     
    def copy(self):
        t = copy.copy(self)
        t.diffuse_source.smodel = self.diffuse_source.smodel.copy()
        return t
    
    @property
    def spectral_model(self):
        return self.diffuse_source.smodel

class DiffuseModelFB(DiffuseModelFromFits): 
    def __init__(self, *pars, **kwargs):
        super(DiffuseModelFB,self).__init__(*pars, **kwargs)


class CacheableDiffuse(convolution.Grid):
    def __init__(self, *args, **kwargs):
        super(CacheableDiffuse,self).__init__(*args, **kwargs)
    
    def make_dict(self):
        return dict(center=self.center, npix=self.npix, pixelsize=self.pixelsize,
                vals = np.array(self.vals, np.float32),
        )
    def dump(self,filename):
        pickle.dump(self.make_dict(), open(filename, 'w'))
    
    def show(self, fignum=1, ax=None, nocolorbar=False, notitle=False, **kwargs):
        """ A simple plot for a sanity check """
        import pylab as plt
        if ax is None:
            plt.close(fignum);
            fig = plt.figure(fignum, figsize=(5,5)); ax = plt.gca() 
        t = np.log10(self.vals)
        imshow_kw = dict(interpolation='nearest', origin='upper'); dict.update(kwargs)
        mappable = ax.imshow(t.transpose()[::-1],**imshow_kw)
        marker = float(self.npix)/2
        ax.axvline(marker,color='k')
        ax.axhline(marker,color='k')
        if not notitle:
            ax.set_title('l, b= (%06.2f,%+06.2f)'%(self.center.l(), self.center.b()),fontsize='medium')
        if not nocolorbar: plt.colorbar(mappable,ax=ax)
        
    def load(self, filename):
        t = pickle.load(open(filename))
        self.__dict__.update(t)
 
class CacheDiffuseConvolution(object):
    """ 
    Create convolved cache
    Deprecated, remove at some point
    
    """
    defaults = DiffuseModel.defaults + (
        ('pixelsize',0.25,'Pixel size for convolution grid'),
        ('r_multi',1.0,"Multiple of r95 to set max dimension of grid"),
        ('r_max',20,"An absolute maximum (half)-size of grid (deg)"),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, exposure, diffuse_function, **kwargs):
        """ 
            exposure : skymaps.Exposure object
            diffuse_filename :string
        """
        keyword_options.process(self, kwargs)
        self.df = diffuse_function
        self.bg = skymaps.Background(self.df,exposure[0],exposure[1])
        
    def convolve(self, band):
        """ band : an ROIBand object
        returns a convolution.Grid object
        """
        psf = band.psf.parent_psf # use this for compatibilty
        roi_dir = band.sd
        energy = band.e #only at geometric average energy
        conversion_type = band.ct 
        
        multi = 1 + 0.01*(energy==band.emin) -0.01*(energy==band.emax)
        r95 = psf.inverse_integral(energy*multi,conversion_type,95)
        #r95 = psf.inverse_integral_on_axis(0.95)
        rad = self.r_multi*r95 + np.degrees(band.radius_in_rad)
        rad = max(min(self.r_max,rad),np.degrees(band.radius_in_rad)+2.5)
        npix = int(round(2*rad/self.pixelsize))
        npix += (npix%2 == 0)
        bgc = CacheableBackgroundConvolution(roi_dir,self.bg, psf,
                npix=npix, pixelsize=self.pixelsize)
        bgc.setup_grid(npix,self.pixelsize)
        bgc.do_convolution(energy,conversion_type)
        return bgc

    def process_bands(self, bands, dumpto):
        dicts = []
        for band in bands:
            d = dict(energy=band.e, conversion_type=band.ct)
            d.update( self.convolve(band).make_dict())
            dicts.append(d)
        pickle.dump(dicts, open(dumpto, 'w'))
        print 'wrote pickle file %s' %dumpto       

class CacheDiffuseModel(CacheDiffuseConvolution):
    """
    create a cache of the diffuse flux values, for later multiplication by exposure and convolution
    Note that the default npix correspond to M. Kerr's code for 4 bands/decade, back.
    This will be used for front as well.
    """
    defaults =CacheDiffuseConvolution.defaults+ (
        ('npix',      None, ""),
        ('emin',   100., 'minimum enery'),
        ('bpd',      4,  'bands per decade'),
        ('nbands',  14, 'number of bands'),
        ('efactor', 10**(0.125), 'energy factor for evaluation'),
        ('npix_list',    (161,  141,  103,   81,   67,   61,), 'number of pixels: array for bands'),
        )

    @keyword_options.decorate(defaults)
    def __init__(self, diffuse_function, **kwargs):
        """
        parameters
        ----------
        diffuse_function : skymaps.SkySpectrum object
            presumably a 
        """
        keyword_options.process(self, kwargs)
        self.df = diffuse_function
        logemin = np.log10(self.emin)
        self.energies = np.logspace(logemin, logemin+float(self.nbands)/self.bpd, 
            self.nbands+1)  
        
    def fill(self, sdir, energy,  npix):
        """ sdir : Skydir
            energy : energy to evaluate 
            npix : number of pixels in x and y
        returns a convolution.Grid object
        """
        grid = CacheableDiffuse(sdir, npix=npix, pixelsize=self.pixelsize)
        grid.setup_grid(npix,self.pixelsize)
        self.df.setEnergy(energy) 
        grid.vals = grid.fill(self.df)
        if np.any(np.isnan(grid.vals)):
            raise DiffuseException('grid had nan(s)')
        return grid
        
    def process_bands(self, center, dumpto):
        """  
        create and write a list of dictionaries for the energy bands,
        using npix and energies from the 
        """
        dicts = []
        for iband, emin in enumerate(self.energies[:-1]): 
            npix = self.npix_list[min(iband, len(self.npix_list)-1)]
            energy = emin*self.efactor
            d = dict(emin=emin, emax=self.energies[iband+1], energy=energy)
            d.update( self.fill(center, energy, npix).make_dict())
            dicts.append(d)
        pickle.dump(dicts, open(dumpto, 'w'))
        print 'wrote pickle file %s' %dumpto


def create_diffuse_cache(diffuse, outdir, name='ring', **kwargs):
    """
    Create a file for each healpixel containing grids for each energy
    
    """
    diffuse_file = os.path.expandvars(diffuse)
    if not os.path.exists(diffuse_file):
        DiffuseException('Diffuse file %s nof found'%diffuse_file)
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    diffusemodel = skymaps.DiffuseFunction(diffuse_file, **kwargs)
    cdm = CacheDiffuseModel(diffusemodel)
    
    def makeone(hp12):
        skydir = skymaps.Band(12).dir(hp12)
        filename = 'HP12_%04d_%s.pickle'%(hp12,name)
        cdm. process_bands(skydir, os.path.join(outdir,filename))
        
    map(makeone, range(1728))
  

def mapper(roi_factory, roiname, skydir, source, **kwargs): 
    """
    return a DiffuseModel appropriate for the source
    """
    psf = roi_factory.psf
    exposure = roi_factory.exposure
    if source.name.startswith('iso'):
        return IsotropicModel(psf, exposure,skydir, source, **kwargs)
    elif source.name.startswith('limb'):
        for dmodel in source.dmodel:
            if not getattr(dmodel,'loaded', False): dmodel.load()
        if source.smodel.name=='Constant':
            # limb model must be separate front and back!
            source.smodel = models.FrontBackConstant(0.001, 2.0)
        return DiffuseModelFB(psf, exposure,skydir, source, **kwargs)
    return DiffuseModelFromCache(psf, exposure,skydir, source, **kwargs)
        
   