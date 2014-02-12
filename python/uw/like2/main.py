"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.65 2014/02/09 19:22:16 burnett Exp $

"""
import types, time
import numpy as np
from uw.utilities import keyword_options
from skymaps import SkyDir
from . import (views,  configuration, extended,  roimodel, from_xml, from_healpix,
                bands,  localization, sedfuns, tools,
                plotting, associate, printing, to_healpix
        )


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
        ('load_kw', {'rings':2}, 'a dict specific for the loading'),
        ('postpone', False, 'Set True to not load data until requested'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, config_dir, roi_spec, **kwargs):
        """Start pointlike v2 (like2) in the specified ROI
        
        parameters
        ----------
        config_dir : string
            file path to a folder containing a file config.txt
            see configuration.Configuration
        roi_spec : integer or TODO (ra,dec) tuple, name of an xml file
            
        """
        keyword_options.process(self, kwargs)
        self.config=config = configuration.Configuration(config_dir, quiet=self.quiet, postpone=self.postpone)
        ecat = extended.ExtendedCatalog(config.extended)
        if isinstance(roi_spec, int):
            roi_sources = from_healpix.ROImodelFromHealpix(config, roi_spec, ecat=ecat,load_kw=self.load_kw)
            roi_index = roi_spec
            self.name = 'HP12_%04d' % roi_index
        elif isinstance(roi_spec, str):
            roi_sources =from_xml.ROImodelFromXML(config, roi_spec)
            roi_index = roi_sources.index
            self.name = roi_spec
        else:
            raise Exception('Did not recoginze roi_spec: %s' %roi_spec)
        
        roi_bands = bands.BandSet(config, roi_index)
        roi_bands.load_data()
        super(ROI, self).__init__( roi_bands, roi_sources)
    
    def __repr__(self):
        if hasattr(self, 'sources'):
            return '%s.%s :\n\t%s\n\t%s' % (self.__module__, self.__class__.__name__, self.config, self.sources)
        else:
            return '%s.%s :\n\t %s\n\t%s' % (self.__module__, self.__class__.__name__, self.config, 'No sources')
        
    def fit(self, select=None, exclude=None,  summarize=True,  **kwargs):
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

        kwargs 
        ------
        ignore_exception : bool
                if set, run the fit in a try block and return None
        update_by : float
            set to zero to not change parameters, or a number between 0 and 1 to make a partial update
        tolerance : float, default 0.1
            If current fit quality, an estimate of potential improvent of the log likelihood, which is
            based on gradient and hessian is less than this, do not fit
            
            others passed to the fitter minimizer command. defaults are
                estimate_errors = True
                use_gradient = True
                
        """
        ignore_exception = kwargs.pop('ignore_exception', False)
        update_by = kwargs.pop('update_by', 1.0)
        tolerance = kwargs.pop('tolerance', 0.1)
        
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
                if summarize:
                    print '%d calls, function value, improvement, quality: %.1f, %.2f, %.2f'\
                        % (fv.calls, w, w - fv.initial_likelihood, fv.delta_loglike())
                self.fit_info = dict(
                    loglike = fv.log_like(),
                    pars = fv.parameters[:], 
                    covariance  = fv.covariance,
                    qual = qual,)
                fv.modify(update_by)
                fv.save_covariance()
                if summarize: fv.summary()
                
            except Exception, msg:
                print 'Fit Failure %s: quality: %.2f' % (msg, qual)
                fv.summary() # 
                if not ignore_exception: raise
        return 
    
    def summarize(self, select=None, exclude=None):
        """construct a summary of the parameters, or subset thereof
        
        """
        with self.fitter_view(select, exclude=exclude) as fv:
            print 'current likelihood, est. diff to peak: %.1f, %.2f' % (fv.log_like(), fv.delta_loglike())
            fv.summary()
            
    def localize(self, source_name=None, update=False, **kwargs):
        """ localize the source, return elliptical parameters 
        """
        if source_name=='all':
            localization.localize_all(self, **kwargs)
            return
        source = self.sources.find_source(source_name)
        with self.tsmap_view(source.name) as tsm:
            loc = localization.Localization(tsm, **kwargs)
            try: 
                loc.localize()
                t = loc.ellipse
            except Exception, e:
                print 'Failed localization for source %s: %s' % (tsm.source.name, e)
                return None
        if update:
            tsm.source.skydir = skymaps.SkyDir(t['ra'], t['dec'])
    
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

    def get_counts(self):
        return plotting.counts.get_counts(self)
        
    def get_sed(self, source_name=None, event_type=None, update=False, tol=0.1):
        """ return the SED recarray for the source
        source_name : string
            Name of a source in the ROI, with possible wildcards
        event_type : None, or integer, 0/1 for front/back
        update : bool
            set True to force recalculation of sed recarray
        """
        source = self.sources.find_source(source_name)
        if not hasattr(source, 'sedrec') or source.sedrec is None or update:
            with sedfuns.SED(self, source.name) as sf:
                source.sedrec = sf.sed_rec(event_type=event_type, tol=tol)
        return source.sedrec

    def band_ts(self, source_name=None, update=False):
        sedrec = self.get_sed(source_name, update=update)
        return sedrec.ts.sum()
    
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
    def plot_tsmap(self, source_name=None, **kwargs):
        """ create a TS map showing the source localization
        """
        source = self.sources.find_source(source_name)
        plot_kw = dict(size=0.25, pixelsize=0.25/15, outdir=None, 
            assoc=getattr(source, 'associations', None) ) 
        plot_kw.update(kwargs)
        with self.tsmap_view(source.name) as tsm:

            loc = localization.Localization(tsm)
            try: 
                if not hasattr(source,'ellipse'):
                    loc.localize()
                    loc.summary()
                tsize = kwargs.pop('size', source.ellipse[2]*15.) # scale according to major axis s
                plot_kw.update(size=tsize, pixelsize=kwargs.pop('pixelsize', tsize/15.))
            except Exception, e:
                print 'Failed localization for source %s: %s' % (source.name, e)
            tsp = plotting.tsmap.plot(loc, **plot_kw)
        return tsp.axes.figure # might want access to TSplot.
    
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
        
        
    def find_associations(self, source_name=None, classes='all_but_gammas', quiet=True):
        """ find associations, using srcid object.
        If source was not localized, run that first
        Expect to find files at $FERMI/catalog.
        """
        if not hasattr(self, 'srcid'):
            self.srcid=associate.SrcId('$FERMI/catalog',classes, quiet=quiet)

        def find_one(source):
            if not hasattr(source, 'ellipse'):
                self.localize(source.name)
            with self.tsmap_view(source.name) as tsv:
                associate.make_association(source, tsv, self.srcid)

        if source_name=='all':
            sources = filter(lambda s: np.any(s.model.free) and not s.isglobal and not s.isextended, self.sources)
            for source in sources:
                find_one(source)
        else:
            source = self.sources.find_source(source_name)
            find_one(source)
            
                
    def to_healpix(self, pickle_dir, dampen, **kwargs):
        to_healpix.pickle_dump(self, pickle_dir, dampen=dampen, **kwargs)
        
    
    def to_xml(self, filename):
        """Save the current ROI to an XML file
        """
        return self.sources.to_xml(filename)

class MultiROI(ROI):
    """ROI subclass that will perform a fixed analysis on multiple ROIs
    Intended for subclasses to override the proc function
    """
    
    def __init__(self, config_dir,  quiet=False):
        """
        """
        self.config = configuration.Configuration(config_dir, quiet=quiet, postpone=False)
        self.ecat = extended.ExtendedCatalog(self.config.extended)
            
    def setup_roi(self, roi_index):
        roi_bands = bands.BandSet(self.config, roi_index)
        roi_bands.load_data()
        roi_sources = from_healpix.ROImodelFromHealpix(self.config, roi_index, ecat=self.ecat,)
        self.name = 'HP12_%04d' % roi_index
        self.setup( roi_bands, roi_sources)
