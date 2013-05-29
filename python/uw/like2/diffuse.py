"""
Provides classes to encapsulate and manipulate diffuse sources.

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/diffuse.py,v 1.23 2013/02/12 15:20:50 burnett Exp $

author: Matthew Kerr, Toby Burnett
"""
import sys, os, types, pickle, glob, copy, zipfile
import numpy as np
import pandas as pd
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
        ('binsperdec',    4,           'bins per decade'),
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
        
    def show(self, **kwargs):
        import pylab as plt
        title = kwargs.pop('title', None)
        fig, axx = plt.subplots(1,3, figsize=(8,3), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.05)
        axx[0].imshow(self.psf_vals,interpolation='nearest')
        axx[0].set_aspect(1.0)
        norm = plt.normalize(np.log10(self.bg_vals.min()),np.log10(self.bg_vals.max()))
        marker = float(self.npix)/2
        for ax,what in zip(axx[1:], (self.bg_vals, self.cvals)  ):
            ax.imshow(np.log10(what).transpose()[::-1], norm=norm, interpolation='nearest')
            ax.axvline(marker,color='k')
            ax.axhline(marker,color='k')

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
        self.efactor = 10**(0.5/self.binsperdec)
    
    def setup(self):
        dfun = self.diffuse_source.dmodel[0]
        try:
            filename = dfun.filename
        except AttributeError:
            filename = dfun.name()
        cache_path = os.path.splitext(filename)[0]+'_%dbpd'%self.binsperdec
        cache_path_alt = os.path.splitext(filename)[0]
        assert os.path.exists(filename), 'oops, %s not found' %filename
        if os.path.exists(cache_path+'.zip'):
            z = zipfile.ZipFile(cache_path+'.zip')
            t = z.namelist()
            if len(t)==1729: # if ziped with foldername
                t = t[1:]
            files = sorted(t) 
            opener = z.open
        elif os.path.exists(cache_path_alt+'.zip'):
            z = zipfile.ZipFile(cache_path_alt+'.zip')
            t = z.namelist()
            if len(t)==1729: # if ziped with foldername
                t = t[1:]
            files = sorted(t) 
            opener = z.open
        else:
            if not os.path.exists(cache_path):
                cache_path = cache_path_alt
                if not os.path.exists(cache_path):
                    raise DiffuseException('cache folder or zip %s not found' %cache_path)
            files = sorted(glob.glob(cache_path+'/*'))
            opener=open
        assert len(files)==1728, 'wrong number of files: expected 1728, found %d' % len(files)
        try:
            self.filename = files[self.diffuse_source.index]
            self.cached_diffuse = pickle.load(opener(self.filename))
        except:
            raise DiffuseException( 'Diffuse cache file # %d not found' %self.diffuse_source.index)

        if not self.quiet: print 'Using cached diffuse in %s'%self.filename
        self.emins = [cd['emin'] for cd in self.cached_diffuse]
        if hasattr(dfun, 'kw') and dfun.kw is not None: # check for extra keywords
            self.corr = pd.read_csv(dfun.kw['correction'])
            self.systematic = dfun.kw['systematic']
            print '\tusing file "%s" for corrections' % dfun.kw['correction']
            print '\tsystematic factor:%.3f' % dfun.kw['systematic']

        
    def make_grid(self, energy, conversion_type):
        """ return a convovlved grid
        
        parameters
        ----------
        energy : float
            intermediate energy for the band
        conversion_type : int
            0 or 1 for front or back
      
        """
        # find the appropriate cached grid
        for index in range(len(self.emins)):
            if energy>self.emins[index] and (index==len(self.emins)-1 or energy<self.emins[index+1])\
                : break
        emin = self.emins[index]    
        assert energy/emin < 1.8 and energy> emin, 'too large a factor: energy, emin=%.0f,%.0f\nemins=%s' % (energy, emin, self.emins)
        cd = self.cached_diffuse[index]
        
        # create a convolvable grid from it, with current PSF
        energy = cd['energy']
        grid = ConvolvableGrid(cd['center'], None, self.psf,  npix=cd['npix'], pixelsize=cd['pixelsize'])
            
        # determine the current exposure for this energy and conversion type
        exp = self.exposure[conversion_type]
        exp.setEnergy(energy)
        
        # get values of exposure on the grid, multiply by saved diffuse
        vals = grid.fill(exp) * cd['vals'] 

        # finally do the convolution on the product of exposure and diffuse map, which is passed in 
        grid.do_convolution(energy, conversion_type, override_vals=vals)
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
        cflux = dm(self.roi_dir, energy)
        if np.isnan(cflux):
            print 'WARNING: nan flux from diffuse map: using zero'
            cflux=0
        grid.cvals = grid.fill(exp) * cflux 
        nnan = np.sum(np.isnan(grid.cvals))
        if nnan>0: print 'Grid for %s has %d nan values ' %( dm.name(), nnan)
        assert nnan<5, \
            'Grid for %s, %.0f Mev< has %d >5 nan values ' %( dm.name(), energy, nnan)
        return grid

class DiffuseModelFromFits( DiffuseModel):

    defaults =DiffuseModel.defaults+ (
        ('pixelsize',0.25,'Pixel size for convolution grid'),
        ('npix_list',    (161,  141,  103,   81,   67,   61,), 'number of pixels: array for bands'),
        ('ignore_nan', True, 'replace nan values with zero (generates warning)'), 
        )

    """ load pattern from a FITS file """
    def __init__(self, *pars, **kwargs):
        keyword_options.process(self, kwargs)
        super(DiffuseModelFromFits,self).__init__(*pars, **kwargs)
        
    def setup(self):
        if not self.quiet:
            print 'set up %d diffuse_models from Fits files:' %len(self.diffuse_source.dmodel)
            for dm in self.diffuse_source.dmodel:
                print '\t%s'% dm.name()
                
    def make_grid(self, energy, conversion_type):
        """ return a convovlved grid
        
        parameters
        ----------
        energy : float
            intermediate energy for the band (perhaps derive it here)
        conversion_type : int
            0 or 1 for front or back
       
        """
        # get the diffuse map, set energy
        dmodels = self.diffuse_source.dmodel
        dm = dmodels[conversion_type if len(dmodels)>1 else 0]
        dm.setEnergy(energy)
        
        # get the exposure
        exp = self.exposure[conversion_type]
        exp.setEnergy(energy)
        
        # create a grid with size determine by npix_list, appropriate PSF function
        iband = min(int(np.log10(energy/100.)*4), len(self.npix_list)-1)
        npix = self.npix_list[iband] 
        
        grid = ConvolvableGrid(self.roi_dir, None, self.psf,  npix=npix, pixelsize=self.pixelsize)
        # set it up with product of exposure and diffuse map
        grid.cvals = grid.fill(exp) * grid.fill(dm) #product of exposure and map
        # check for nans, replace with zeros if not full ROI
        nans = np.isnan(grid.cvals)
        if np.all(nans):
            raise DiffuseException('Diffuse cube %s has no overlap with ROi' % dm.filename)
        if np.any(nans) and self.ignore_nan:
            grid.cvals[nans]=0
        
        # finally convolve it
        grid.do_convolution(energy, conversion_type, override_vals=grid.cvals) 
        return grid
     
    def copy(self):
        t = copy.copy(self)
        t.diffuse_source.smodel = self.diffuse_source.smodel.copy()
        return t
    
    @property
    def spectral_model(self):
        return self.diffuse_source.smodel

class DiffuseModelFB(IsotropicModel): #DiffuseModelFromFits): 
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
        ('binsperdec',      4,  'bands per decade'),
        ('decades',     3.5, 'number of decades'),
        ('npix_list',    (161,  141,  103,   81,   67,   61,), 'number of pixels: array for bands'),
        ('ignore_nan', True, 'replace nan values with zero (generates warning)'), 
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
        if self.binsperdec==8:
            self.npix_list = (161,161,141,141,103,103,81,81,67,67, 61,61)
            self.nbands=28
        elif self.binsperdec!=4:
            raise DiffuseException('binsperdec parameter=%d: must be either 4 or 8'%self.binsperdec)
        self.efactor = 10**(0.5/self.binsperdec)
        self.nbands = int(self.binsperdec*self.decades)
        self.df = diffuse_function
        logemin = np.log10(self.emin)
        self.energies = np.logspace(logemin, logemin+float(self.nbands)/self.binsperdec, 
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
            vals = grid.vals.flatten()
            ok = ~np.isnan(vals)
            print '%d/%d nans for energy %.0f:  min %.2g, mean %.2g, max %.2g:'%(sum(np.isnan(vals)), len(vals), energy,
                    vals[ok].min(), vals[ok].mean(), vals[ok].max())
            if self.ignore_nan:
                grid.vals[np.isnan(grid.vals)]=0
            else:
                raise DiffuseException('Convolution problem?: energy=%.1f, %d/%d nan(s)'\
                    % (energy, np.sum(np.isnan(grid.vals)), len(grid.vals)) )
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
            try:
                d.update( self.fill(center, energy, npix).make_dict())
            except DiffuseException, msg:
                print 'exception in process_bands, emin=%.0f'%emin
                raise
            
            dicts.append(d)
        pickle.dump(dicts, open(dumpto, 'w'))
        print 'wrote pickle file %s with %d bands' % (dumpto, self.nbands)


def create_diffuse_cache(name, **kwargs):
    """
    Create a file for each healpixel containing grids for each energy
    
    name : string
        name of the global diffuse. Unless specific file specified, expect to find a FITS file name_*fit*
    
    optional parameters
        outdir :string
            actual path to folder in which to create files (eventually zip?)
            default determined from name of fits file found
        idlist : int
            number of HEALPix rois, default 1728
        binsperdec : int
            either 4 or 8, default 4
        infile : string
            name file to analyze. Default None, meaning expect to find with glob
       others passed on to skymaps.DiffuseFunction
    """
    
    diffuse_dir = kwargs.pop('outdir', os.path.expandvars('$FERMI/diffuse')) 
    assert os.path.exists(diffuse_dir), 'diffuse directory %s not found' %diffuse_Dir
    binsperdec = kwargs.pop('binsperdec',4)
    idlist = kwargs.pop('idlist', 1728) 

    diffuse_file = kwargs.pop('infile', None)
    if diffuse_file is None:
        pattern = os.path.join(diffuse_dir, name+'_*.fit*')
        diffuse_file = glob.glob(pattern)
        if len(diffuse_file)==0:
            raise DiffuseException('no files found with pattern %s'%pattern)
        elif len(diffuse_file)>1:
            raise DIffuseException("""more than one input file found, %s:
                specify which one with infile"""%diffuse_file)
        else: diffuse_file = diffuse_file[0]
    if not os.path.exists(diffuse_file):
        DiffuseException('Diffuse file %s nof found'%diffuse_file)
    
    outdir = kwargs.pop('outdir', None)
    if outdir is None:
        outdir = os.path.join(diffuse_dir,
                os.path.split(diffuse_file)[-1].split('.')[0])+'_%dbpd'%binsperdec
        if not os.path.exists(outdir):
            os.mkdir(outdir)

    print 'will create cache from: %s\n\t\t    to: %s' % (diffuse_file, outdir)
    if kwargs.pop('test', False): return

    # load the FITS diffuse file
    diffusemodel = skymaps.DiffuseFunction(diffuse_file, **kwargs)
    
    cdm = CacheDiffuseModel(diffusemodel, binsperdec=binsperdec)
    
    def makeone(hp12):
        skydir = skymaps.Band(12).dir(hp12)
        filename = 'HP12_%04d_%s.pickle'%(hp12,name)
        cdm. process_bands(skydir, os.path.join(outdir,filename))
        
    map(makeone, range(idlist))
  

def mapper(roi_factory, roiname, skydir, source, **kwargs): 
    """
    return a convolved DiffuseModel appropriate for the source
    """
    psf = roi_factory.psf
    exposure = roi_factory.exposure
    try:
        print '>>>>Loading', source.name, source.dmodel[0].name()
        if source.name.lower().startswith('iso'):
            return IsotropicModel(psf, exposure,skydir, source, **kwargs)
        elif source.name.startswith('limb'):
            for dmodel in source.dmodel:
                if not getattr(dmodel,'loaded', False): dmodel.load()
            if source.smodel.name=='Constant':
                # limb model must be separate front and back!
                source.smodel = models.FrontBackConstant(0.001, 2.0)
            return DiffuseModelFB(psf, exposure,skydir, source, **kwargs)
        elif source.name=='ring':
            for dmodel in source.dmodel:
                if not getattr(dmodel,'loaded', False): dmodel.load()
            if os.path.splitext(source.dmodel[0].name())[-1]=='.fits':
                return DiffuseModelFromFits(psf, exposure, skydir, source, **kwargs)
        ### cache is default, it seems
        return DiffuseModelFromCache(psf, exposure,skydir, source, **kwargs)
    except Exception, msg:
        print 'Failed to map diffuse source "%s" : %s' % (source.name, msg)
        raise
        
   
