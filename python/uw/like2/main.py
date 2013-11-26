"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.40 2013/11/26 04:25:47 burnett Exp $

"""
import numpy as np
from . import (views,  configuration, extended,  roimodel, bands,  
                localization, tools, plotting,
        )


class ROI_user(views.LikelihoodViews):
    """ subclass of ROIstat that adds user-interface tools: fitting, localization, plotting sources, counts
    methods
    =======
    fit
    localize
    get_source
    get_model
    summary
    TS, band_ts
    get_sed
    plot_tsmap
    plot_sed
    plot_counts
    print_summary
    find_associations
    add_source
    del_source
    zero_source
    unzero_source
    freeze
    thaw
    set_model
    """
    def __init__(self, config_dir, roi_index, **kwargs):
        config = configuration.Configuration(config_dir, quiet=True, postpone=True)
        ecat = extended.ExtendedCatalog(config.extended)
        roi_sources = roimodel.ROImodel(config, ecat=ecat, roi_index=roi_index)
        roi_bands = bands.BandSet(config, roi_index)
        roi_bands.load_data()
        super(ROI_user, self).__init__( roi_bands, roi_sources)

       
    def fit(self, select=None, exclude=None,  summarize=True, quiet=True, **kwargs):
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
            ignore_exceptions : bool
                if set, run the fit in a try block and return None
            call_limit : int
                if set, modify default limit on number of calls
            others passed to the fitter minimizer command. defaults are
                estimate_errors = True
                use_gradient = True
        """
        ignore_exception = kwargs.pop('ignore_exception', False)
        #self.call_limit = kwargs.pop('call_limit', self.call_limit)
        fit_kw = dict(use_gradient=True, estimate_errors=True)
        fit_kw.update(kwargs)

        wfit = None; pfit=None
        with self.fitter_view(select, exclude=exclude) as fv:
            try:
                fv.maximize(**fit_kw)
                w = fv.log_like()
                if summarize:
                    print '%d calls, function value, improvement: %.0f, %.1f'\
                        % (fv.calls, w, w - fv.initial_likelihood)
                if fit_kw['estimate_errors'] :
                    pass
                    #self.sources.set_covariance_matrix(mm.cov_matrix, selected)
                if summarize: fv.summary()
                wfit = fv.log_like()
                pfit = fv.parameters[:] 
            except FloatingPointError, e:
                if not ignore_exception: raise
            except Exception:
                if ignore_exception: return 
                else: raise
        return wfit, pfit
    
    def localize(self, source_name=None, update=False, **kwargs):
        """ localize the source, return elliptical parameters 
        """

        with self.tsmap_view(source_name, **kwargs) as tsm:
            loc = localization.Localization(tsm, **kwargs)
            try: 
                loc.localize()
                t = loc.ellipse
            except Exception, e:
                print 'Failed localization for source %s: %s' % (tsm.source.name, e)
                return None
        if update:
            tsm.source.skydir = skymaps.SkyDir(t['ra'], t['dec'])
        return t
    
    def get_model(self, source_name=None):
        return self.sources.find_source(source_name).spectral_model
    def get_source(self, source_name=None):
        return self.sources.find_source(source_name)
        
    def TS(self, source_name=None, quick=True):
        """ measure the TS for the given source and set the ts property of the source
        """
        source = self.sources.find_source(source_name)
        with self.fitter_view(source_name) as fv:
             source.ts  = fv.ts()
        return source.ts
#
#    def band_ts(self, source_name=None):
#        sed = self.get_sed(source_name)
#        return np.sum(sed.ts)
#        
    def get_sed(self, source_name=None, event_type=None, update=False, **kwargs):
        """ return the SED recarray for the source
        source_name : string
            Name of a source in the ROI, with possible wildcards
        event_type : None, or integer, 0/1 for front/back
        """
        source = self.sources.find_source(source_name)
        if not hasattr(source, 'sedrec') and not update:
            with sedfuns.SED(self, source.name, **kwargs) as sf:
                source.sedrec = sf.sed_rec(event_type=event_type)
        return source.sedrec

    @tools.decorate_with(plotting.sed.stacked_plots)    
    def plot_sed(self, source_name=None, **kwargs):
        source = self.sources.find_source(source_name)
        showts = kwargs.pop('showts', True)
        if kwargs.get('update', True):
            self.get_sed(update=True)
        annotation =(0.04,0.88, 'TS=%.0f' % source.ts ) if showts else None 
        kwargs.update(galmap=self.roi_dir, annotate=annotation)
        return plotting.sed.stacked_plots(self, **kwargs); 

        
    @tools.decorate_with(plotting.counts.stacked_plots)
    def plot_counts(self, **kwargs):
        return plotting.counts.stacked_plots(self, **kwargs)
        
#    @decorate_with(pointlike_plotting.tsmap.plot)
#    def plot_tsmap(self, source_name=None, **kwargs):
#        """ create a TS map showing the source localization
#        """
#        source = self.sources.find_source(source_name)
#        plot_kw = dict(size=0.25, pixelsize=0.25/15, outdir=None, 
#            assoc=getattr(source, 'adict', None) ) 
#        plot_kw.update(kwargs)
#        with localization.Localization(self, source.name, quiet=True) as loc:
#            try: 
#                if not hasattr(source,'ellipse'):
#                    loc.localize()
#                    loc.summary()
#                tsize = kwargs.pop('size', source.ellipse[2]*15.) # scale according to major axis size
#                plot_kw.update(size=tsize, pixelsize=tsize/15.)
#            except Exception, e:
#                print 'Failed localization for source %s: %s' % (source.name, e)
#            tsp = pointlike_plotting.tsmap.plot(loc, **plot_kw)
#        return tsp
#    @decorate_with(printing.print_summary)
#    def print_summary(self, **kwargs):
#        selected_source= self.sources.selected_source
#        t = printing.print_summary(self, **kwargs)
#        self.sources.selected_source=selected_source
#
#    def find_associations(self, source_name=None, classes='all', quiet=True):
#        """ find associations, using srcid object.
#        Expect to find files at $FERMI/catalog.
#        """
#        try:
#            if not hasattr(self,'srcid'):
#                from uw.like2.pipeline import associate
#                self.srcid=associate.SrcId('$FERMI/catalog',classes, quiet=quiet)
#            source = self.sources.find_source(source_name)
#            with localization.Localization(self, source.name, quiet=quiet) as loc:
#                if not hasattr(source, 'ellipse'):
#                    loc.localize()
#                localization.make_association(source, loc.TSmap, self.srcid)
#        except Exception, msg:
#            print 'Failed to find associations for %s: %s' % (self.get_source().name, msg)
#            
#    @property
#    def bounds(self):
#        return self.sources.bounds
#    @property
#    def cov_matrix(self):
#        """ the current covariance matrix, determined from the gradient 
#            assumes that there is a consistent fit for for the full set of parameters
#        """
#        return fitter.Minimizer.mycov(self.gradient, self.get_parameters())
#        
#    @property
#    def correlations(self):
#        """Return the linear correlation coefficients for the estimated covariance matrix.
#        """
#        cm = self.cov_matrix
#        diag = cm.diagonal()
#        diag[diag<=0]=np.nan # protect against zero or negative
#        s = np.sqrt(diag)
#        return cm / np.outer(s,s)
#        
#    def add_source(self, **kwargs):
#        """ add a source to the ROI
#        keywords:
#            name : string
#            model
#            skydir
#        """
#        newsource = sources.PointSource(**kwargs)
#        self.sources.add_source(newsource)
#        for band in self.all_bands:
#            band.add_source(newsource)
#        self.initialize()
#        return self.get_source(newsource.name)
#        
#    def del_source(self, source_name):
#        """ delete the specifiec source (which can be expressed with wildcards """
#        source = self.sources.del_source(source_name)
#        for band in self.all_bands:
#            band.del_source(source)
#        self.initialize()
#        
#    def zero_source(self, source_name=None, **kwargs):
#        model = self.get_model(source_name)
#        model.zero()
#        self.initialize()
#        
#    def unzero_source(self, source_name=None, **kwargs):
#        model = self.get_model(source_name)
#        model.unzero()
#        self.initialize()
#    
#    def poisson_likelihood(self, source_name=None, **kwargs):
#        """ return a functor that corresponds to the likelihood function for the given source
#        """
#        source = self.sources.find_source(source_name)
#
#        with sedfuns.SourceFlux(self, source.name) as sf:
#            sf.select_band(None)
#            pf = loglikelihood.PoissonFitter(sf, **kwargs)
#            p = pf.poiss
#        return p 
#        
#    def freeze(self, parname, source_name=None, value=None):
#        """ freeze the parameter, optionally set the value
#        
#        parname : name or index
#        source_name: None or string
#            if None, use currently selected source
#        value : float or None
#            if float, set the value
#        """
#        model = self.get_model(source_name)
#        model.freeze(parname)
#        if value is not None: model.setp(parname, value)
#        self.initialize()
#        
#    def thaw(self, parname, source_name=None):
#        """ thaw the parameter
#        
#        parname : name or index
#            if a string with an underscore, interpret as source_parname
#        source_name: None or string
#            if None, use currently selected source
#        """
#        if parname.find('_')>0 and source_name is None:
#            source_name, parname = parname.split('_')
#        model = self.get_model(source_name)
#        model.freeze(parname, freeze=False)
#        self.initialize()
#        
#    def set_model(self, model, source_name=None):
#        """ replace the current model, return reference to previous
#        
#        model : string, or like.Models.Model object
#            if string, evaluate. Note that 'PowerLaw(1e-11,2.0)' will work. Also supported:
#            ExpCutoff, PLSuperExpCutoff, LogParabola, each with all parameters required.
#        source_name: None or string
#            if None, use currently selected source
#        """
#        from sources import ExpCutoff, PLSuperExpCutoff, LogParabola, PowerLaw
#        src = self.get_source(source_name)
#        old_model = src.spectral_model
#        if type(model)==types.StringType:
#            model = eval(model) 
#        #assert model.isinstance(Models.Model), 'model must inherit from Model class'
#        sourcelist.set_default_bounds(model)
#        src.spectral_model = model
#        self.initialize()
#        return old_model
        
