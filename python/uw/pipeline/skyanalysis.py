"""
Basic ROI analysis
$Header$
"""
import os, pickle, glob, types
import numpy as np
from skymaps import SkyDir, PySkyFunction
from uw.pipeline import skymodel, dataspec
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
        fit_emin     = 175.,
        fit_emax     = 600000.,
        minROI       = 5,
        maxROI       = 5,
        radius       =10, 
        irf          = 'P6_v11_diff',
        fit_kw = dict(fit_bg_first = False,
            use_gradient = True, ),
        log          = None,
        )

    def __init__(self, sky, dataset='24MP7source', **kwargs):
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
            s += '\t%-20s: %s\n' %(key, self.__dict__[key])
        return s

        s = 'SkyAnalysis analysis environment:\n'
        s += super(SkyAnalysis,self).__str__()
        return s
    
    def _diffuse_sources(self, index ):
        """ return a source manager for the diffuse,  global and extended sources
        """
        skydir = self.skymodel.skydir(index)
        # get all diffuse models appropriate for this ROI
        globals, extended = self.skymodel.get_diffuse_sources(index, radius=self.radius)
       
        # perform OTF convolutions with PSFs
        diffuse_mapper = lambda x: roi_diffuse.ROIDiffuseModel_OTF(self, x, skydir)
        global_models = map(diffuse_mapper, globals)
        
        extended_mapper =lambda x: roi_extended.ROIExtendedModel.factory(self,x,skydir)
        extended_models = map(extended_mapper, extended)
        
        # create and return the manager
        return roi_managers.ROIDiffuseManager(global_models+extended_models, skydir, quiet=self.quiet)
        
    def _local_sources(self, index ):
        """ return a manager for the local sources with significant overlap with the ROI
        """
        ps = self.skymodel.get_point_sources(index, self.radius)
        skydir = self.skymodel.skydir(index)
        return roi_managers.ROIPointSourceManager(ps, skydir,quiet=self.quiet)
        
    def roi(self, index):
        """ return a roi_analysis.ROIAnalysis object based on the roi index
        """
        ps_manager = self._local_sources( index)
        bg_manager = self._diffuse_sources(index)
        
        def iterable_check(x):
            return x if hasattr(x,'__iter__') else (x,x)

        r = ROIanalysis(ps_manager.roi_dir, 
                    ps_manager, bg_manager, 
                    self, 
                    name = 'HP12_%04d' % index,
                    fit_emin=iterable_check(self.fit_emin), 
                    fit_emax=iterable_check(self.fit_emax),
                    quiet=self.quiet, 
                    fit_kw = self.fit_kw)
        return r

class ROIanalysis(roi_analysis.ROIAnalysis):
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
        super(ROIanalysis, self).__init__(*pars, **kwargs)
        
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
            super(ROIanalysis, self).fit( **fit_kw)
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
            loc, i, delta, deltaTS= super(ROIanalysis,self).localize(which=which,bandfits=bandfits,
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
    
    def get_model(which):
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
        self.localizer = roi_localize.ROILocalizer(self, which, bandfits=bandfits)
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
        
