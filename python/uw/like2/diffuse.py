"""
Provides classes to encapsulate and manipulate diffuse sources.

$Header$

author: Matthew Kerr, Toby Burnett
"""
import sys, os, types, pickle
import numpy as np
from uw.utilities import keyword_options, convolution
from uw.like import roi_diffuse #until not needed
import skymaps #from Science Tools


###=========================================================================###
class DiffuseModel(object):
    """ Associate a SkySpectrum with a spectral model and an energy band
    
        Provide the interface that DiffuseModel et al. should satisfy."""

    defaults = (
        ('quiet', False, 'Set True to quiet'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, band, 
                diffuse_source, 
                roi_dir,
                **kwargs):
        
        self.band = band
        self.roi_dir = roi_dir
        self.diffuse_source = diffuse_source
        keyword_options.process(self, kwargs)
        self.setup()
   
    def setup(self): pass
    
    def __str__(self):
        return '%s scaled with %s\n'%(self.name,self.smodel.pretty_name)+self.smodel.__str__()

    def initialize_counts(self):
        """ This method is responsible for establishing a state from
            which the model counts can be calculated.  E.g., evaluating
            the model on a series of energy subplanes and storing the
            results for a later Simpson's rule.
        """
        raise NotImplementedError,'Classes must implement this method!'

    def update_counts(self ):
        """ This method *must* set the following members of each band:
            band.bg_counts[model_index] -- the total counts expected for
                    the model in the aperture
            band.bg_pix_counts[:,model_index] -- if the band has pixels,
                    an npixel vector with the expected counts from the
                    model for each data pixel
        """
        raise NotImplementedError,'Classes must implement this method!'

    def gradient(self):
        """ This method should return the gradient with respect to the
            free parameters of the spectral model.
        """
        raise NotImplementedError,'Classes must implement this method!'



###====================================================================================================###

class DiffuseModelFromCache(DiffuseModel):
    cache = None
    defaults = DiffuseModel.defaults + (
        ('pixelsize',0.25,'Pixel size for convolution grid'),
        ('npix',101,"Note -- can be overridden at the band level"),
        ('nsimps',4,"Note -- some energies use a multiple of this"),
        ('r_multi',1.0,"Multiple of r95 to set max dimension of grid"),
        ('r_max',20,"An absolute maximum (half)-size of grid (deg)"),
        )

    class GridFromDict(convolution.BackgroundConvolution):
        """ subclass of BackgroundConvolution that fills the grid from a dict"""
        
        def __init__(self, cdict, bounds_error=True,fill_value=0):
            self.set_center(cdict['center'])
            self.setup_grid(cdict['npix'],cdict['pixelsize'])
            if False: #cdict['energy']<188 and cdict['conversion_type']==0:
                factor=0.91
                print '\n*** correcting front low energy exposure by %.2f factor\n' %factor
                self.cvals = cdict['cvals'] * factor
                self.bg_vals=cdict['bg_vals'] * factor
            else:
                self.cvals = cdict['cvals'] 
                self.bg_vals=cdict['bg_vals'] 
        
            self.bounds_error=bounds_error
            self.fill_value=fill_value
            
    def __init__(self, roi_factory, source, skydir, **kwargs):
        roiname = kwargs.pop('name')
        self.filename = os.path.join(DiffuseModelFromCache.cache, roiname+'_ring.pickle')
        assert os.path.exists(self.filename), 'Diffuse cache file %s not found' %filename
        keyword_options.process(self, kwargs)

        self.sa = roi_factory
        self.diffuse_source = source
        self.smodel = source.smodel
        self.dmodel = source.dmodel
        self.grid = None
        self.name = source.name
        
    def _ap_value(self,center,radius):
        """ Return the integral of the model over the aperture    """
        return self.grid.ap_average(radius)

    def _pix_value(self,pixlist):
        """ Return the model evaluated at each data pixel """
        return self.grid(pixlist,self.grid.cvals)
    
    def initialize_counts(self, band): 
        self.band=band
        energy = self.energy = band.e
        delta_e = band.emax - band.emin
        dicts=pickle.load(open(self.filename))
        for d in dicts:
            if d['energy']!=energy or d['conversion_type']!=self.band.ct: continue
            self.grid = DiffuseModelFromCache.GridFromDict(d)
            break
        if self.grid is None: 
            raise Exception('Cached band convolution not found for energy %f type %d' %(energy,band.ct))
 
        self.ap_evals = self._ap_value(band.sd,band.radius_in_rad) * band.solid_angle * delta_e
        if not band.has_pixels: return
        self.pi_evals = self._pix_value(band.wsdl) if band.has_pixels else []
        self.pi_evals *= band.b.pixelArea() * delta_e

    def update_counts(self):
        # calculate integral counts 
        self.counts = self.ap_evals * self.smodel(self.energy)
        self.pix_counts = self.pi_evals * self.smodel(self.energy) if self.band.has_pixels else []
        

###====================================================================================================###
class CacheableBackgroundConvolution(convolution.BackgroundConvolution):
    """
    
    """
    def __init__(self, *args, **kwargs):
        if len(args)==1 and type(args[0])==types.StringType:
            t = pickle.load(open(args[0]))
            self.set_center(t['center'])
            self.setup_grid(t['npix'],t['pixelsize'])
            self.bg_vals = t['bg_vals']
            self.cvals = t['cvals']
            self.psf_vals = t['psf_vals']
        else:
            super(CacheableBackgroundConvolution,self).__init__(*args, **kwargs)
    
    def make_dict(self):
        return  dict(center=self.center, npix=self.npix, pixelsize=self.pixelsize,
                    bg_vals= np.array(self.bg_vals, dtype=np.float32), 
                    cvals=   np.array(self.cvals, dtype=np.float32),
                    psf_vals = np.array(self.psf_vals, dtype=np.float32))

    def dump(self,filename):
        pickle.dump( self.makedict(), open(filename, 'w') )
        
class CacheDiffuseConvolution(object):
    """ 
    
    """
    defaults = DiffuseModel.defaults + (
        ('pixelsize',0.25,'Pixel size for convolution grid'),
        ('npix',101,"Note -- can be overridden at the band level"),
        ('nsimps',4,"Note -- some energies use a multiple of this"),
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

class IsoDiffuseModel(roi_diffuse.ROIDiffuseModel_PC):
    def initialize_counts(self, band):
        self.band = band
        super(IsoDiffuseModel,self).initialize_counts([band], self.skydir)

    def update_counts(self, index=1):
        super(IsoDiffuseModel,self).update_counts([self.band], index)
        self.counts = self.band.bg_counts[index]
        if self.band.has_pixels:
            self.pix_counts = self.band.bg_pix_counts[:,index]

class GalDiffuseModel(roi_diffuse.ROIDiffuseModel_OTF):
    def initialize_counts(self, band):
        self.band = band
        super(GalDiffuseModel,self).initialize_counts([band], self.skydir)

    def update_counts(self, index=0):
        super(GalDiffuseModel,self).update_counts([self.band], index)
        self.counts = self.band.bg_counts[index]
        if self.band.has_pixels:
            self.pix_counts = self.band.bg_pix_counts[:,index]

def mapper(roi_factory, name, skydir, source, **kwargs): 
    """
    return a DiffuseModel appropriate for the source
    """
    cache = kwargs.pop('cache', '/phys/groups/tev/scratch1/users/Fermi/cached_diffuse/P72Y')
    
    if source.name.startswith('iso'):
        return IsoDiffuseModel(roi_factory, source, skydir, **kwargs)
    
    elif cache is not None and source.name.startswith('ring'):
        # cache only for ring
        DiffuseModelFromCache.cache = cache
        kwargs.update(name=name)
        return DiffuseModelFromCache(roi_factory,  source, skydir, **kwargs)
    else:
        # here if not using cache or not ring (limb perhaps)
        for dmodel in source.dmodel:
            if not getattr(dmodel,'loaded', False): dmodel.load()
        return GalDiffuseModel(roi_factory, source, skydir, **kwargs)

        
   