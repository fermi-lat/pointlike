"""
Basic ROI analysis
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pipeline/skyanalysis.py,v 1.5 2011/01/25 03:24:14 kerrm Exp $
"""
import os, pickle, glob, types
import numpy as np
from skymaps import SkyDir, PySkyFunction
from . import skymodel, dataspec
from uw.utilities import keyword_options
from uw.like import pointspec, roi_analysis, roi_managers, roi_diffuse, roi_localize
from uw.like import sed_plotter, tsmap_plotter, counts_plotter
from uw.like import roi_extended #ExtendedSource,ROIExtendedModel

# special function to replace or extend a docstring from that of another function
def decorate_with(other_func, append=False):
    def decorator(func):
        if append: func.__doc__ += other_func.__doc__ 
        else:      func.__doc__  = other_func.__doc__ 
        return func
    return decorator


class SkyAnalysis(pointspec.SpectralAnalysis):

    config = dict(
        fit_emin     = 100,
        fit_emax     = 600000.,
        minROI       = 5,
        maxROI       = 5,
        radius       =10, 
        irf          = 'P6_v11_diff',
        convo_reso   = 0.25,
        fit_kw = dict(fit_bg_first = False,
            use_gradient = True, ),
        log          = None,
        )

    def __init__(self, sky, dataset, **kwargs):
        """
        sky:  skymodel object
        dataset:
        """
        # extract updates to local kw
        
        for kw in self.config:
            if kw in kwargs: self.config[kw]=kwargs.pop(kw)
        default_keys = [x[0] for x in self.defaults]
        # add keys to modify to kwargs
        for kw in default_keys: 
            if kw in self.config: kwargs[kw]= self.config[kw]
        super(SkyAnalysis,self).__init__( dataspec.DataSpec(dataset), **kwargs)
        # now add waht is left
        self.__dict__.update(self.config)
        self.skymodel = sky
        if not self.quiet: 
            print >>self.log, self
            if self.log is not None: self.log.close()
    
    def __str__(self):
        s = '%s configuration:\n'% self.__class__.__name__
        show = """CALDB irf skymodel dataspec fit_emin fit_emax fit_kw 
               minROI maxROI process_kw""".split()
        for key in show:
            s += '\t%-20s: %s\n' %(key, 
                self.__dict__[key] if key in self.__dict__.keys() else 'not in self.__dict__!')
        return s

        s = 'SkyAnalysis analysis environment:\n'
        s += super(SkyAnalysis,self).__str__()
        return s
    
    def _diffuse_sources(self, src_sel):
        """ return a source manager for the diffuse,  global and extended sources
        """
        skydir = src_sel.skydir()
        # get all diffuse models appropriate for this ROI
        globals, extended = self.skymodel.get_diffuse_sources(src_sel)
       
        # perform OTF convolutions with PSFs
        diffuse_mapper = lambda x: roi_diffuse.ROIDiffuseModel_OTF(self, x, skydir, pixelsize=self.convo_reso)
        global_models = map(diffuse_mapper, globals)
        
        extended_mapper =lambda x: roi_extended.ROIExtendedModel.factory(self,x,skydir)
        extended_models = map(extended_mapper, extended)
        
        # create and return the manager
        return roi_managers.ROIDiffuseManager(global_models+extended_models, skydir, quiet=self.quiet)
        
    def _local_sources(self, src_sel):
        """ return a manager for the local sources with significant overlap with the ROI
        """
        ps = self.skymodel.get_point_sources(src_sel)
        skydir = src_sel.skydir()
        return roi_managers.ROIPointSourceManager(ps, skydir,quiet=self.quiet)
        
    def roi(self, index, src_sel=None):
        """ return a roi_analysis.ROIAnalysis object based on the roi index
        """
        if src_sel is None:
            src_sel = HEALPixSourceSelector(index,self.skymodel,max_radius=self.radius)
        else:
            if (src_sel.max_radius != self.radius):
                print 'Warning, overriding max radius of %.2f with %.2f'%(src_sel.max_radius,self.radius)
                src_sel.max_radius = self.radius

        ps_manager = self._local_sources(src_sel)
        bg_manager = self._diffuse_sources(src_sel)
        
        def iterable_check(x):
            return x if hasattr(x,'__iter__') else (x,x)

        r = PipelineROI(ps_manager.roi_dir, 
                    ps_manager, bg_manager, 
                    self, 
                    name = 'HP12_%04d' % index,
                    fit_emin=iterable_check(self.fit_emin), 
                    fit_emax=iterable_check(self.fit_emax),
                    quiet=self.quiet, 
                    fit_kw = self.fit_kw)
        return r

  
class SourceSelector(object):
    """ Manage inclusion of sources in an ROI."""
    
    defaults = (
        ('max_radius',10,'Maximum radius (deg.) within which sources will be selected.'),
        ('free_radius',3,'Radius (deg.) in which sources will have free parameters'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, skydir, **kwargs):
        self.mskydir = skydir
        keyword_options.process(self,kwargs)

    def include(self,source):
        """ source -- an instance of Source """
        return source.near(self.mskydir,self.max_radius)

    def free(self,source):
        """ source -- an instance of Source """
        return source.near(self.mskydir,self.free_radius)

    def frozen(self,source): return not self.free(source)

    def skydir(self): return self.mskydir
        
class HEALPixSourceSelector(SourceSelector):
    """ Manage inclusion of sources in an ROI based on HEALPix.
    Overrides the free method to define HEALpix-based free regions
    """

    @keyword_options.decorate(SourceSelector.defaults)
    def __init__(self, index, skymodel,  **kwargs):
        """ index : int
                HEALpix index for the ROI
            skymodel : an instance of SkyModel
        """
        keyword_options.process(self,kwargs)
        self.index = index
        self.skymodel = skymodel
        self.mskydir =  skymodel.skydir(index)

    def free(self,source):
        """
        source : instance of skymodel.Source
        -> bool, if this source in in the region where fit parameters are free
        """
        return self.skymodel.index(source.skydir) == self.index

    

class PipelineROI(roi_analysis.ROIAnalysis):
    """ sub class of the standard ROIAnalysis class to cusomize the fit, add convenience functions
    """

    def __init__(self, *pars, **kwargs):
        self.fit_kw = kwargs.pop('fit_kw', dict())
        self.likelihood_count=0
        self.prior = lambda x : 0 # default, no prior
        self.name = kwargs.pop('name', None)
        if self.name is None:
            if len(self.psm.point_sources)>0:
                self.name=self.psm.point_sources[0].name
            else: self.name='(not set)'
        self.center= pars[0] #roi_dir
        super(PipelineROI, self).__init__(*pars, **kwargs)
        
    def logLikelihood(self, parameters, *args):
        """ the total likelihood, according to model
            parameters parameters to pass to model
        """
        
        self.likelihood_count +=1
        if np.any(np.isnan(parameters)):
            # pretty ridiculous that this check must be made, but fitter passes NaNs...
            return 1e6
            # not sure if should "set parameters" in this case

        self.update_counts(parameters)

        ll = sum(band.logLikelihood(phase_factor=self.phase_factor) for band in self.bands)
        if np.isnan(ll) : return 1e6
        return ll -self.prior(self.psm.models)

    def fit(self, **kwargs):
        """ invoke base class fitter, but insert defaults first 
        """
        ignore_exception = kwargs.pop('ignore_exception', True)
        fit_kw = self.fit_kw
        fit_kw.update(kwargs)
        ts = 0
        initial_count = self.likelihood_count
        initialL = self.logl
        try:
            super(PipelineROI, self).fit( **fit_kw)
            ts = self.TS()
        except Exception, msg:
            if not self.quiet: print 'Fit failed: %s' % msg
            if not ignore_exception: raise
        if not self.quiet:
            print 'logLikelihood called %d times, change: %.1f' % (self.likelihood_count - initial_count, initialL-self.logl )
        return ts
        
    def dump(self, sdir=None, galactic=False, maxdist=5, title=''):
        """ formatted table point sources positions and parameter in the ROI"""
        self.print_summary(sdir, galactic, maxdist, title)
    
    @decorate_with(sed_plotter.plot_sed)
    def plot_sed(self, **kwargs):
        return sed_plotter.plot_sed(self,**kwargs)

    @decorate_with(counts_plotter.plot_counts)
    def plot_counts(self, **kwargs):
        return counts_plotter.plot_counts(self, **kwargs)
    
    @decorate_with(tsmap_plotter.plot_tsmap)
    def plot_tsmap(self, *pars, **kwargs):
        return tsmap_plotter.plot_tsmap(self, *pars, **kwargs)
     
    def band_ts(self, which=0):
        """ return the sum of the individual band ts values
        """
        self.setup_energy_bands()
        ts = 0
        for eb in self.energy_bands:
            eb.bandFit(which)
            ts += eb.ts
        return ts

    def localize(self,which=0, tolerance=1e-3,update=False, verbose=False, bandfits=True, seedpos=None):
        """Localize a source using an elliptic approximation to the likelihood surface.

          which     -- index of point source; default to central 
                      **if localizing non-central, ensure ROI is large enough!**
          tolerance -- maximum difference in degrees between two successive best fit positions
          update    -- if True, update localization internally, i.e., recalculate point source contribution
          bandfits  -- if True, use a band-by-band (model independent) spectral fit; otherwise, use broabband fit
          seedpos   -- use for a modified position (pass to superclass)

         return fit position, change in TS
        """
        try:
            quiet, self.quiet = self.quiet, not verbose # turn off details of fitting
            loc, i, delta, deltaTS= super(PipelineROI,self).localize(which=which,bandfits=bandfits,
                            tolerance=tolerance,update=update,verbose=verbose, seedpos=seedpos)
            self.quiet = quiet
            if not self.quiet: 
                name = self.psm.point_sources[which].name if type(which)==types.IntType else which
                print 'Localization of %s: %d iterations, moved %.3f deg, deltaTS: %.1f' % \
                    (name,i, delta, deltaTS)
                self.print_ellipse()
        except Exception, e:
            print 'Localization failed! %s' % e
            self.qform=None
            loc, deltaTS = None, 99 
        #self.find_tsmax()
        return loc, deltaTS
    
    def get_model(self,which):
        """ return a reference to the model
            which : integer or string
            
        """
        psm, index = self.mapper(which) #raise exception if wrong.
        return psm.models[index]
  
    def tsmap(self, which=0, bandfits=True):
        """ return function of likelihood in neighborhood of given source
            tsm = roi.tsmap(which)
            size=0.25
            tsp = image.TSplot(tsm, center, size, pixelsize =size/20, axes=plt.gca())
            tsp.plot(center, label=name)
            tsp.show()
        """
        self.localizer = roi_localize.localizer(self, which, bandfits=bandfits)
        return PySkyFunction(self.localizer)
        
    def fit_ts_list(self, which=0):
        """ return breakdown of final fit ts per band """
        self.zero_ps(which)
        self.update_counts(self.get_parameters())
        w0 = np.array([band.logLikelihood() for band in self.bands])
        self.unzero_ps(which)
        self.update_counts(self.get_parameters())
        w1 = np.array([band.logLikelihood() for band in self.bands])
        return 2*(w0-w1)
        
