"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.91 2018/01/27 15:37:17 burnett Exp $

"""
import types, time, glob
import numpy as np
import pandas as pd
from uw.utilities import keyword_options
from skymaps import SkyDir, Band
from . import (views,  configuration, extended,  roimodel, from_xml, from_healpix,
                bands,  localization, sedfuns, tools,
                plotting, associate, printing, to_healpix
        )

import warnings
warnings.resetwarnings()
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)
import astropy
warnings.filterwarnings('ignore', category=astropy.wcs.FITSFixedWarning)
warnings.filterwarnings('ignore', category=astropy.io.fits.verify.VerifyWarning,)

class ROI(views.LikelihoodViews):
    """ROI analysis
    This is a list of the properties and functions appropriate for user analysis of an ROI
    
    properties
    ----------
    all_energies -- list of all band energies
    bands -- BandSet object, a list of all thebands
    emin, emax -- range of energies used by currently selected bands
    energies -- list of central energies
    quiet  -- flag to suppress most output
    roi_dir -- center of the ROI
    sources -- ROImodel object, managing the source list
    selected -- list of selected bands

    likelihood-related functions
    ----------------------------
    fit  -- perform a fit
    gradient -- return the gradient
    hessian -- return the hessian, a n x n matrix
    log_like -- value of the log likelihood
    summarize -- a table of free parameters, values, errors, gradient
    delta_loglike -- estimate of possible improvement.


    source-related functions
    ------------------------
      All of these but add_source take the name of a source, which is an optional parameter: 
      if a source has already been selected, it will be used
    TS -- TS for current source
    add_source -- add a new source
    band_ts -- band TS 
    del_source -- delete a source
    find_associations -- find associations for current source
    get_model -- the model
    get_sed -- return SED information for the current source
    get_counts -- return array of counts per energy band
    Npred -- return predicted total counts for current source
    get_source -- select and return a source
    freeze -- freeze a parameter
    set_model -- change the spectral model for a source
    thaw -- thaw (unfreeze) a parameter


    plot/print
    ----------
    plot_counts
    plot_sed
    plot_tsmap
    plot_roi_position
    print_summary
    print_sed

    change bands
    ------------
    select -- select a subset of the bands, which will be used subsequently to define likelihood functions
    """
    
    defaults = (
        ('quiet', True, 'set to suppress output'),
        ('load_kw', {'rings':2, 'tsmin':0}, 'a dict specific for the loading'),
        ('postpone', False, 'Set True to not load data until requested'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, config_dir, roi_spec=None, xml_file=None, **kwargs):
        """Start pointlike v2 (like2) in the specified ROI
        
        parameters
        ----------
        config_dir : string
            file path to a folder containing a file config.txt
            see configuration.Configuration
        roi_spec : [None |integer | (ra,dec) tuple ]
            If None, require that the input_model dict has a key 'xml_file'
            if an integer, it must be <1728, the ROI number
            if a string, assume a source name and load the ROI containing it
            
        """
        keyword_options.process(self, kwargs)
        self.config=config = configuration.Configuration(config_dir, quiet=self.quiet, postpone=self.postpone)
        ecat = extended.ExtendedCatalog(config.extended)
        
        if isinstance(roi_spec, str):
            # try:
            #     roi_sources =from_xml.ROImodelFromXML(config, roi_spec)
            #     roi_index = roi_sources.index
            # except:
            #     print 'No ROI specification (an index) or presence of an xml file'
            #     raise
            # Change to just expecting the name of a source
            sourcelist=glob.glob('sources_*.csv')[0]
            df = pd.read_csv(sourcelist, index_col=3 if roi_spec[0]=='J' else 0)
            if roi_spec not in df.index:
                print 'Source name "{}" not found '.format(roi_spec)
                raise Exception
            roi_index = int(df.loc[roi_spec]['roiname'][-4:]) 
            print 'Loading ROI #{}, containing source "{}"'.format(roi_index, roi_spec)
        elif isinstance(roi_spec, int):
            roi_index = roi_spec
        elif type(roi_spec)==tuple and len(roi_spec)==2:
            roi_index = Band(12).index(SkyDir(*roi_spec))

        else:
            raise Exception('Did not recoginze roi_spec: %s' %(roi_spec))             
            
        roi_sources = from_healpix.ROImodelFromHealpix(config, roi_index, ecat=ecat,)
        config.roi_spec = configuration.ROIspec(healpix_index=roi_index)

        self.name = config.roi_spec.name if config.roi_spec is not None else roi_spec
        
        roi_bands = bands.BandSet(config, roi_index)
        roi_bands.load_data()
        super(ROI, self).__init__( roi_bands, roi_sources)
    
    def __repr__(self):
        if hasattr(self, 'sources'):
            return '%s.%s :\n\t%s\n\t%s' % (self.__module__, self.__class__.__name__, self.config, self.sources)
        else:
            return '%s.%s :\n\t %s\n\t%s' % (self.__module__, self.__class__.__name__, self.config, 'No sources')
        
    def fit(self, select=None, exclude=None,  summarize=True, setpars=None, **kwargs):
        """ Perform fit, return fitter object to examine errors, or refit
        
        Parameters
        ----------
        select : None, item or list of items, where item is an int or a string
            if not None, it defines a subset of the parameter numbers to select
                to define a projected function to fit
            int:  select the corresponding parameter number
            string: select parameters according to matching rules
                The name of a source (with possible wild cards) to select for fitting
                If initial character is '_', match the rest with parameter names
                if initial character is '_' and last character is '*', treat as wild card
        
        exclude : None, int, or list of int 
            if specified, will remove parameter numbers from selection

        summarize : bool
            if True (default) call summary after succesful fit

        setpars : dict | None
            set a set of parameters by index: the dict has keys that are either the index, 
            or the name of the varialbe, and float values,
            e.g. {1:1e-14, 2:2.1, 'Source_Index': 2.0}
            Note that this uses *internal* variables

        kwargs 
        ------
        ignore_exception : bool
                if set, run the fit in a try block and return None
        update_by : float
            set to zero to not change parameters, or a number between 0 and 1 to make a partial update
        tolerance : float, default 0.0
            If current fit quality, an estimate of potential improvent of the log likelihood, which is
            based on gradient and hessian is less than this, do not fit
        plot : bool
            if set to True, create plots of the parameter fits

            others passed to the fitter minimizer command. defaults are
                estimate_errors = True
                use_gradient = True
                
        """
        if len(self.sources.parameters)==0:
            print 'No parameters to fit'
            return
        ignore_exception = kwargs.pop('ignore_exception', False)
        update_by = kwargs.pop('update_by', 1.0)
        tolerance = kwargs.pop('tolerance', 0.0)
        plot = kwargs.pop('plot', False)
       
        if setpars is not None: 
            self.sources.parameters.setitems(setpars)
            
        fit_kw = dict(use_gradient=True, estimate_errors=True)
        fit_kw.update(kwargs)

        with self.fitter_view(select, exclude=exclude) as fv:
            qual = fv.delta_loglike()
            if qual < tolerance and qual>0:
                if summarize:
                    print 'Not fitting, estimated improvement, %.2f, is less than tolerance= %.1f' % (qual, tolerance)
                    return
            try:
                fv.maximize(**fit_kw)
                w = fv.log_like()
		self.fmin_ret = fv.fmin_ret
                if summarize:
                    print '%d calls, function value, improvement, quality: %.1f, %.2f, %.2f'\
                        % (fv.calls, w, w - fv.initial_likelihood, fv.delta_loglike())
                self.fit_info = dict(
                    loglike = fv.log_like(),
                    pars = fv.parameters[:], 
                    covariance  = fv.covariance,
                    mask = fv.mask,
                    qual = qual,)
                fv.modify(update_by)
                if fit_kw['estimate_errors']: fv.save_covariance()
                if summarize: fv.summary()
                if plot: fv.plot_all()
                
            except Exception, msg:
                print 'Fit Failure %s: quality: %.2f' % (msg, qual)
                fv.summary() # 
                if not ignore_exception: raise
        return 
    
    
    
    def summarize(self, select=None, exclude=None):
        """construct a summary of the parameters, or subset thereof
        
        """
        if len(self.sources.parameters)==0:
            print 'No parameters to fit'
            return
        with self.fitter_view(select, exclude=exclude) as fv:
            print 'current likelihood, est. diff to peak: %.1f, %.2f' % (fv.log_like(), fv.delta_loglike())
            fv.summary()
    
    def profile(self, source_name=None, set_normalization=False):
        """Return profile info as a dict. Add profile info to source object
        
            source_name : str | 'all'
            set_normalization : bool
                if True, use peak likelihood to set normalization of the model
        """
        if source_name=='all':
            for s in self.free_sources:
                if s.isglobal: continue
                self.profile(s.name)
            return
        source = self.sources.find_source(source_name)
        with sedfuns.SED(self, source.name) as sf:
            try:
                p,maxdev= sf.full()
            except Exception, msg:
                print 'Failed profile fit for {}: {}'.format(source.name, msg)
                source.profile=None
                return 
            err=p.errors
            t= dict(peak=p.flux, low=err[0], high=err[1], ts=p.ts, zf=p.zero_fraction())
            err= np.array(p.errors)/p.flux-1 if p.flux>0 else 0
            
            if not self.quiet:
                str = '{:20} flux'.format(source.name)
                if p.flux>0 and p.errors[0]>0:
                    print '{} = {:6.3f} (1 + {:.3f} - {:.3f}) eV/cm^2/s, TS={:.1f} '.format(
                        str, p.flux,err[1],-err[0], p.ts)
                else:
                    print '{} < {:6.3f}'.format(str, t['high'])

            source.profile = t
        if set_normalization:
            m = source.model
            m[0]=t['peak']/(m.e0**2 * 1e06) if t['peak']>0 else 1e-15
        return t


    def localize(self, source_name=None, update=False, ignore_exception=True, **kwargs):
        """ localize the source, return elliptical parameters 
        """
        if source_name=='all':
            localization.localize_all(self, ignore_exception=ignore_exception, **kwargs)
            return
        source = self.sources.find_source(source_name)
        with self.tsmap_view(source.name) as tsm:
            loc = localization.Localization(tsm, **kwargs)
            try: 
                loc.localize()
                t = loc.ellipse
            except Exception, e:
                print 'Failed localization for source %s: %s' % (tsm.source.name, e)
                if ignore_exception: return None
                raise 
        if update:
            tsm.source.skydir = SkyDir(t['ra'], t['dec'])
    
    def get_model(self, source_name=None):
        return self.sources.find_source(source_name).spectral_model
        
    def get_source(self, source_name=None):
        return self.sources.find_source(source_name)
        
    def TS(self, source_name=None, quick=True):
        """ measure the TS for the given source and set the ts property of the source
        """
        source = self.sources.find_source(source_name)
        with self.fitter_view(source.name) as fv:
             source.ts  = fv.ts()
        return source.ts

    def get_count_dict(self, event_type=None):
        """ return the count analysis dict
        """
        return plotting.counts.get_counts(self, event_type=event_type)
        
    def get_sed(self, source_name=None, event_type=None, update=False, tol=0.1):
        """ return the SED recarray for the source, including npred info
        source_name : string
            Name of a source in the ROI, with possible wildcards
        event_type : None, or integer, 0/1 for front/back, 2-5 for psf0-3
        update : bool
            set True to force recalculation of sed recarray
        """
        source = self.sources.find_source(source_name)
        if not hasattr(source, 'sedrec') or source.sedrec is None\
                 or (update and np.any(source.model.free)):
            with sedfuns.SED(self, source.name) as sf:
                source.sedrec = sf.sed_rec(event_type=event_type, tol=tol)
        
        return source.sedrec

    def band_ts(self, source_name=None, update=False):
        sedrec = self.get_sed(source_name, update=update)
        return sedrec.ts.sum()
    
    def get_counts(self, source_name=None, event_type=None):
        """ return the array of predicted counts, or npred, for the specified, or all sources if 'all'
        """
        if source_name=='all': 
            return plotting.counts.get_counts(self, event_type=event_type)['models']

        source = self.sources.find_source(source_name)
        return plotting.counts.get_npred(self, source.name)
        
    def Npred(self, source_name=None):
        """ Return the total predicted counts for the specified, or default source
        """
        return sum(self.get_counts(source_name))
    
    @tools.decorate_with(plotting.sed.stacked_plots)    
    def plot_sed(self, source_name=None, **kwargs):
        if source_name=='all':
            #flag to do them all
            sedfuns.makesed_all(self, **kwargs)
            return
        source = self.sources.find_source(source_name)
        showts = kwargs.pop('showts', True)
        if kwargs.pop('update', False) or not hasattr(self,'sedrec') or self.sedrec is None:
            self.get_sed(update=True)
        annotation =(0.04,0.88, 'TS=%.0f' % source.ts ) if showts and hasattr(source, 'ts') else None 
        kwargs.update(galmap=self.roi_dir, annotate=annotation)
        with sedfuns.SED(self, source.name) as sf:
            t = plotting.sed.stacked_plots(sf, **kwargs)
        return t

    @tools.decorate_with(plotting.counts.stacked_plots)
    def plot_counts(self, **kwargs):
        return plotting.counts.stacked_plots(self, **kwargs)
        
    @tools.decorate_with(plotting.tsmap.plot)
    def plot_tsmap(self, source_name=None, tsplot=False, factor=1.0, refit=False, **kwargs):
        """ create a TS map showing the source localization
        """
        source = self.sources.find_source(source_name)
        plot_kw = dict(size=0.25, pixelsize=0.25/15, outdir=None, 
            assoc=getattr(source, 'associations', None) ) 
        plot_kw.update(kwargs)
        with self.tsmap_view(source.name) as tsm:

            loc = localization.Localization(tsm, factor=factor)
            try: 
                if refit or not hasattr(source,'ellipse') or source.ellipse is None:
                    loc.localize()
                    loc.summary()
                tsize = kwargs.pop('size', source.ellipse[2]*15.) # scale according to major axis s
                plot_kw.update(size=tsize, pixelsize=kwargs.pop('pixelsize', tsize/15.))
            except Exception, e:
                print 'Failed localization for source %s: %s' % (source.name, e)
            tsp = plotting.tsmap.plot(loc, **plot_kw)
        return tsp if tsplot else tsp.axes.figure # might want access to TSplot.
    
    def plot_roi_position(self, ax=None):
        """ define an Axes with a grid showing the position of this ROI """
        from uw.utilities import image
        from matplotlib import pyplot as plt
        if ax is None:
            fig, ax = plt.subplots(figsize=(2,1))
        else:
            fig  = ax.figure
        ag = image.AIT_grid(axes=ax, labels=False)
        ag.plot([self.roi_dir], c='r', s=30)
        return fig
    
    @tools.decorate_with(printing.print_summary)
    def print_summary(self, **kwargs):
        selected_source= self.sources.selected_source
        t = printing.print_summary(self, **kwargs)
        self.sources.selected_source=selected_source

    def print_sed(self, source_name=None):
        """print a formated table of SED info"""
        sedfuns.print_sed(self, source_name)
        
        
    def find_associations(self, source_name=None, classes='all_but_gammas', srcid='srcid', quiet=True):
        """ find associations, using srcid object.
        If source was not localized, run that first
        Expect to find files at $FERMI/catalog.
        """
        if not hasattr(self, 'srcid'):
            self.srcid=associate.SrcId('$FERMI/catalog',classes=classes, srcid=srcid, quiet=quiet)

        def find_one(source):
            try:
                if not hasattr(source, 'ellipse'):
                    self.localize(source.name)
                with self.tsmap_view(source.name) as tsv:
                    associate.make_association(source, tsv, self.srcid)
            except Exception, msg:
                print 'Exception raised while associating souurce %s: %s' %(source.name, msg)

        if source_name=='all':
            sources = filter(lambda s: np.any(s.model.free) and not s.isglobal and not s.isextended, self.sources)
            for source in sources:
                find_one(source)
        else:
            source = self.sources.find_source(source_name)
            find_one(source)
            
                
    def to_healpix(self, pickle_dir, dampen, **kwargs):
        to_healpix.pickle_dump(self, pickle_dir, dampen=dampen, **kwargs)
        
    
    def to_xml(self, filename, **kwargs):
        """Save the current ROI to an XML file. 
        Note that it will include info on the diffuse components not consistent with gtlike
        """
        return self.sources.to_xml(filename, **kwargs)
        
    def move(self, new_location, source_name=None):
        """Move the position of a source
        new_location : SkyDir object | (ra,dec) tuple
        source_name : name of source | None
            if None, use currently selected source
        """
        source = self.get_source(source_name)
        tsv = self.tsmap_view(source.name)
        old_loc = tsv.saved_skydir
        loc = new_location if isinstance( new_location,SkyDir,) else SkyDir(*new_location)
        print 'moved position of source %s from %s to %s, change in TS: %.2f'\
                % (source.name, old_loc, loc, tsv(loc) )
        tsv.set_dir(loc)
        
    def ts_beta(self, source_name=None, ignore_exception=True): 
        """evaluate ts_beta for a Log Parabola source
        
            returns TS(beta_fit)-TS(beta=0), or None if the model is not LogParabola
        """
        fit_pars = dict(tolerance=0, ignore_exception=ignore_exception)
        source = self.sources.find_source(source_name)
        if source.model.name != 'LogParabola': return None
        model_copy = source.model.copy()
        ts_saved=source.ts
        if source.model.free[2]:
            self.fit(source.name, **fit_pars)
            ts2 = self.TS()
            fit_beta=source.model['beta']
            self.freeze('beta', source.name, 0)#reeze at zero
            self.fit(source.name, **fit_pars)
            ts1=self.TS()
            self.thaw('beta', source.name)
            source.model['beta']=fit_beta
        else: 
            #frozen, 
            self.fit(source.name, **fit_pars)
            ts1=self.TS()
            self.thaw('beta', source.name)
            self.fit(source.name, **fit_pars)
            ts2 = self.TS()
            self.freeze('beta', source.name, 0)
        source.model=model_copy #make sure no change
        source.ts=ts_saved
        return ts2-ts1
        

class MultiROI(ROI):
    """ROI subclass that will perform a fixed analysis on multiple ROIs
    Intended for subclasses
    """
    
    def __init__(self, config_dir,  quiet=False, postpone=False, **kwargs):
        """
        """
        self.config = configuration.Configuration(config_dir, quiet=quiet, postpone=postpone,
             **self.config_kw)
        self.ecat = extended.ExtendedCatalog(self.config.extended)

    def setup_roi(self, roi_index):
        roi_bands = bands.BandSet(self.config, roi_index)
        roi_bands.load_data()
        if self.config.modeldir is not None:
            roi_sources = from_healpix.ROImodelFromHealpix(self.config, roi_index, ecat=self.ecat,)
        else:
            roi_sources = from_xml.ROImodelFromXML(self.config, roi_index, ecat=self.ecat)
            
        self.name = 'HP12_%04d' % roi_index
        self.setup( roi_bands, roi_sources)
