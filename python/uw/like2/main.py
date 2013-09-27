"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.36 2013/09/04 12:34:58 burnett Exp $

"""
import types
import numpy as np
from . import roistat, localization, printing, roisetup, sedfuns, sources, loglikelihood, sourcelist
from . import plotting as pointlike_plotting
from uw.utilities import fitter
import skymaps

# special function to replace or extend a docstring from that of another function
def decorate_with(other_func, append=False, append_init=False):
    """ append_init: If decorating with an object (which has an __init__ function),
                     then append the doc for the __init__ after the doc for the
                     overall class. """
    def decorator(func):
        if append: func.__doc__ += other_func.__doc__ 
        else:      func.__doc__  = other_func.__doc__ 

        if append_init and hasattr(other_func,'__init__'):
                func.__doc__ += other_func.__init__.__doc__
        return func
    return decorator


#

class ROI_user(roistat.ROIstat, fitter.Fitted):
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
    def __init__(self, *pars, **kwargs):
        source_name= kwargs.pop('source_name', None)
        super(ROI_user, self).__init__(*pars, **kwargs)
        self.selected_source = self.sources.find_source(source_name) if source_name is not None else None

    def _selector(self, select=None, exclude=None):
        # select a list of parameter numbers, or None for all free parameters
        selected= set()
        npars = len(self.get_parameters())
        
        if select is not None:
            selectpar = select
            if not hasattr(select, '__iter__'): select = [select]
            for item in select:
                if type(item)==types.IntType or type(item)==np.int64:
                    selected.add(item)
                    if item>=npars:
                        raise Exception('Selected parameter number, %d, not in range [0,%d)' %(item, npars))
                elif type(item)==types.StringType:
                    if item.startswith('_'):
                        # look for parameters
                        if item[-1] != '*':
                            toadd = filter( lambda i: self.parameter_names[i].endswith(item), range(npars) )
                        else:
                            def filt(i):
                                return self.parameter_names[i].find(item[:-1])!=-1
                            toadd = filter( filt, range(npars) )
                    else:
                        src = self.get_source(item)
                        toadd = filter(lambda i: self.parameter_names[i].startswith(src.name), range(npars))
                    selected = selected.union(toadd )
                else:
                    raise Exception('fit parameter select list item %s, type %s, must be either an integer or a string' %(item, type(item)))
            select = sorted(list(selected))
            if len(select)==0:
                raise Exception('nothing selected for fit from selection "%s"' % selectpar)
        
        if exclude is not None:
            if not hasattr(exclude, '__iter__'): exclude = [exclude]
            all = set(range(npars)) if select is None else set(select)
            select = list( all.difference(exclude))
                
        return select
       
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
        self.call_limit = kwargs.pop('call_limit', self.call_limit)
        fit_kw = dict(use_gradient=True, estimate_errors=True)
        fit_kw.update(kwargs)
        self.update()
        initial_value, self.calls = self.log_like(), 0
        saved_pars = self.get_parameters()
        
        ## process select: return None or a list of parameter numbers, create the function to fit
        selected = self._selector(select, exclude)
        fn = self if selected is None  else fitter.Projector(self, select=selected)
        try:
            mm = fitter.Minimizer(fn, quiet=quiet)
            mm(**fit_kw)
            w = self.log_like()
            if summarize:
                print '%d calls, function value, improvement: %.0f, %.1f'\
                    % (self.calls, w, w - initial_value)
            if fit_kw['estimate_errors'] :
                self.sources.set_covariance_matrix(mm.cov_matrix, selected)
            if summarize: self.summary(selected)
            return mm
        except FloatingPointError, e:
            print 'Fit error: restoring parameters since  "%s"' %e 
            self.set_parameters(saved_pars)
            if not ignore_exception: raise
            return mm
        
        except Exception:
            if ignore_exception: return None
            else: raise
    
    def cov_matrix(self, par=None):
        return fitter.Minimizer.mycov(self, self.get_parameters() if par is None else par)
        
    def localize(self, source_name=None, **kwargs):
        """ localize the source, return elliptical parameters 
        """
        source = self.sources.find_source(source_name)
        with localization.Localization(self, source.name, **kwargs) as loc:
            try: 
                loc.localize()
                t = loc.ellipse
            except Exception, e:
                print 'Failed localization for source %s: %s' % (source.name, e)
                return None
        if kwargs.get('update', True):
            source.skydir = skymaps.SkyDir(t['ra'], t['dec'])
        return t
    
    def get_sources(self):
        return [ s for s in self.sources if s.skydir is not None]
    def get_model(self, source_name=None):
        return self.get_source(source_name).spectral_model
    def get_source(self, source_name=None):
        return self.sources.find_source(source_name)
        
    def summary(self, select=None, exclude=None, out=None, title=None, gradient=True):
        """ summary table of free parameters, values uncertainties gradient
        
        Parameters:
        ----------
        select : list of integers or string
            integers are indices of parameters
            string is the wildcarded name of a source
        out : open file or None
        title: None or string
        gradient: bool
            set False to not print gradient
            
        """
        if title is not None:
            print >>out, title
        fmt, tup = '%-21s%6s%10s%10s', tuple('Name index value error(%)'.split())
        if gradient:
            grad = self.gradient()
            fmt +='%10s'; tup += ('gradient',)
        print >>out, fmt %tup
        prev=''
        selected = self._selector(select, exclude)

        for index, (name, value, rsig) in enumerate(zip(self.parameter_names, self.model_parameters, self.sources.uncertainties)):
            if selected is not None and index not in selected: continue
            t = name.split('_')
            pname = t[-1]
            sname = '_'.join(t[:-1])
            #sname,pname = name.split('_',1)
            if sname==prev: name = len(sname)*' '+'_'+pname
            prev = sname
            fmt = '%-21s%6d%10.4g%10s'
            psig = '%.1f'%(rsig*100) if rsig>0 and not np.isnan(rsig) else '***'
            tup = (name,index, value,psig)
            if gradient:
                fmt +='%10.1f'; tup += (grad[index],)
            print >>out,  fmt % tup
            
    def TS(self, source_name=None, quick=True, zero=1e-20):
        """ measure the TS for the given source
        """
        source = self.sources.find_source(source_name)
        model = source.spectral_model
        norm =model[0]
        model[0]=zero
        self.update()
        llzero = self.log_like()
        model[0]=norm; self.update()
        ts= 2*(self.log_like()-llzero)
        source.ts=max(ts,0)
        return source.ts

    def band_ts(self, source_name=None):
        sed = self.get_sed(source_name)
        return np.sum(sed.ts)
        
    def get_sed(self, source_name=None, event_class=None, **kwargs):
        """ return the SED recarray for the source
        source_name : string
            Name of a source in the ROI, with possible wildcards
        event_class : None, or integer, 0/1 for front/back
        update : if present in kwargs and True, force a new evaluation; otherwise return present one
        """
        source = self.sources.find_source(source_name)
        update = kwargs.pop('update', False)
        if hasattr(source, 'sedrec') and not update:
            return source.sedrec
        with sedfuns.SourceFlux(self, source.name, **kwargs) as sf:
            source.sedrec = sedfuns.SED(sf, event_class=event_class).rec
            return source.sedrec

    @decorate_with(pointlike_plotting.tsmap.plot)
    def plot_tsmap(self, source_name=None, **kwargs):
        """ create a TS map showing the source localization
        """
        source = self.sources.find_source(source_name)
        plot_kw = dict(size=0.25, pixelsize=0.25/15, outdir=None, 
            assoc=getattr(source, 'adict', None) ) 
        plot_kw.update(kwargs)
        with localization.Localization(self, source.name, quiet=True) as loc:
            try: 
                if not hasattr(source,'ellipse'):
                    loc.localize()
                    loc.summary()
                tsize = kwargs.pop('size', source.ellipse[2]*15.) # scale according to major axis size
                plot_kw.update(size=tsize, pixelsize=tsize/15.)
            except Exception, e:
                print 'Failed localization for source %s: %s' % (source.name, e)
            tsp = pointlike_plotting.tsmap.plot(loc, **plot_kw)
        return tsp
        
    @decorate_with(pointlike_plotting.sed.stacked_plots)    
    def plot_sed(self, source_name=None, **kwargs):
        source = self.sources.find_source(source_name)
        showts = kwargs.pop('showts', True)
        if kwargs.get('update', True):
            self.get_sed(update=True)
        annotation =(0.04,0.88, 'TS=%.0f' % source.ts ) if showts else None 
        kwargs.update(galmap=self.roi_dir, annotate=annotation)
        return pointlike_plotting.sed.stacked_plots(self, **kwargs); 

    @decorate_with(pointlike_plotting.counts.stacked_plots)
    def plot_counts(self, **kwargs):
        return pointlike_plotting.counts.stacked_plots(self, **kwargs)
        
    @decorate_with(printing.print_summary)
    def print_summary(self, **kwargs):
        selected_source= self.sources.selected_source
        t = printing.print_summary(self, **kwargs)
        self.sources.selected_source=selected_source

    def find_associations(self, source_name=None, classes='all', quiet=True):
        """ find associations, using srcid object.
        Expect to find files at $FERMI/catalog.
        """
        try:
            if not hasattr(self,'srcid'):
                from uw.like2.pipeline import associate
                self.srcid=associate.SrcId('$FERMI/catalog',classes, quiet=quiet)
            source = self.sources.find_source(source_name)
            with localization.Localization(self, source.name, quiet=quiet) as loc:
                if not hasattr(source, 'ellipse'):
                    loc.localize()
                localization.make_association(source, loc.TSmap, self.srcid)
        except Exception, msg:
            print 'Failed to find associations for %s: %s' % (self.get_source().name, msg)
            
    @property
    def bounds(self):
        return self.sources.bounds
    @property
    def cov_matrix(self):
        """ the current covariance matrix, determined from the gradient 
            assumes that there is a consistent fit for for the full set of parameters
        """
        return fitter.Minimizer.mycov(self.gradient, self.get_parameters())
        
    @property
    def correlations(self):
        """Return the linear correlation coefficients for the estimated covariance matrix.
        """
        cm = self.cov_matrix
        diag = cm.diagonal()
        diag[diag<=0]=np.nan # protect against zero or negative
        s = np.sqrt(diag)
        return cm / np.outer(s,s)
        
    def add_source(self, **kwargs):
        """ add a source to the ROI
        keywords:
            name : string
            model
            skydir
        """
        newsource = sources.PointSource(**kwargs)
        self.sources.add_source(newsource)
        for band in self.all_bands:
            band.add_source(newsource)
        self.initialize()
        return self.get_source(newsource.name)
        
    def del_source(self, source_name):
        """ delete the specifiec source (which can be expressed with wildcards """
        source = self.sources.del_source(source_name)
        for band in self.all_bands:
            band.del_source(source)
        self.initialize()
        
    def zero_source(self, source_name=None, **kwargs):
        model = self.get_model(source_name)
        model.zero()
        self.initialize()
        
    def unzero_source(self, source_name=None, **kwargs):
        model = self.get_model(source_name)
        model.unzero()
        self.initialize()
    
    def poisson_likelihood(self, source_name=None, **kwargs):
        """ return a functor that corresponds to the likelihood function for the given source
        """
        source = self.sources.find_source(source_name)

        with sedfuns.SourceFlux(self, source.name) as sf:
            sf.select_band(None)
            pf = loglikelihood.PoissonFitter(sf, **kwargs)
            p = pf.poiss
        return p 
        
    def freeze(self, parname, source_name=None, value=None):
        """ freeze the parameter, optionally set the value
        
        parname : name or index
        source_name: None or string
            if None, use currently selected source
        value : float or None
            if float, set the value
        """
        model = self.get_model(source_name)
        model.freeze(parname)
        if value is not None: model.setp(parname, value)
        self.initialize()
        
    def thaw(self, parname, source_name=None):
        """ thaw the parameter
        
        parname : name or index
            if a string with an underscore, interpret as source_parname
        source_name: None or string
            if None, use currently selected source
        """
        if parname.find('_')>0 and source_name is None:
            source_name, parname = parname.split('_')
        model = self.get_model(source_name)
        model.freeze(parname, freeze=False)
        self.initialize()
        
    def set_model(self, model, source_name=None):
        """ replace the current model, return reference to previous
        
        model : string, or like.Models.Model object
            if string, evaluate. Note that 'PowerLaw(1e-11,2.0)' will work. Also supported:
            ExpCutoff, PLSuperExpCutoff, LogParabola, each with all parameters required.
        source_name: None or string
            if None, use currently selected source
        """
        from sources import ExpCutoff, PLSuperExpCutoff, LogParabola, PowerLaw
        src = self.get_source(source_name)
        old_model = src.spectral_model
        if type(model)==types.StringType:
            model = eval(model) 
        #assert model.isinstance(Models.Model), 'model must inherit from Model class'
        sourcelist.set_default_bounds(model)
        src.spectral_model = model
        self.initialize()
        return old_model
        
    def diffuse_correction(self, corr=None):
        """ set or return the diffuse correction factors
        corr : array of float or None
             if None, return the current correction factors
             
        """
        bands = self.selected_bands
        bgal = [b[0] for b in bands]
        dc = [x.diffuse_correction for x in bgal]
        if corr is None:
            return np.array(dc)
        for i,c in enumerate(corr):
            dc[2*i]=dc[2*i+1] = c
            for x,y in zip(dc,bgal):
                y.diffuse_correction = x
 
    def weights(self, source_name, data, emin=None):
        """ return an array of weights, signal/total for the source
        
        data : dict-like 
            expect to have keys ENERGY, CONVERSION_TYPE, RA, DEC
            can just be a pyfits.FITS_rec object
        emin : [float | None ]
            minimum energy, if specified. Assign zero weight for events with energy<emin
        """
        self.get_source(source_name)
        weight = np.empty(len(data))
        source_index = self.sources.selected_source_index
        for i, (energy, ct, ra, dec) in enumerate(zip(
                data['ENERGY'],data['CONVERSION_TYPE'], data['RA'], data['DEC'])):
            if emin is not None and energy<emin: 
                weight[i]=0 
            else:
                bl = self.select_band(energy, ct)
                signal,back= bl.counts_in_pixel(source_index, skymaps.SkyDir(float(ra),float(dec)) )
                weight[i]=signal/(signal+back)
        return weight
        
class Factory(roisetup.ROIfactory):
    """ subclass of ROIfactory that sets up a ROI_user analysis object"""
    def __call__(self, sel):
        """
        Select and setup an ROI.
        sel : string or SkyDir
            if string, the name of a source, or coordinate in form J123.9+30.1
        """
        source_name=None
        if type(sel)==types.IntType:
            index = sel
        elif type(sel)==skymaps.SkyDir:
            index = self.skymodel.hpindex(sel)
        elif type(sel)==types.StringType:
            source = self.skymodel.find_source(sel)
            if source is not None:
                skydir = source.skydir
            elif sel[0]=='J':
                # starts with 'J': try coordinate like J123.6-60.1
                t = sel[1:]
                i = max(t.find('+'), t.find('-'))
                skydir = skymaps.SkyDir(float(t[0:i]),float(t[i:]))
            else:
                raise Exception('Source %s not found in model' % sel)
            index = self.skymodel.hpindex(skydir)
            source_name=sel
        else:
            raise Exception( 'factory argument "%s" not recognized.' %sel)
        roi = ROI_user(super(Factory,self).__call__(index))
        ## preselect the given source after setting up the ROI
        if source_name is not None: 
            roi.sources.find_source(source_name)
            assert source_name==roi.sources.selected_source.name
        return roi


@decorate_with(roisetup.ROIfactory, append_init=True)
def factory(  modeldir='.', dataspec=None,     **kwargs    ):
    """ will then return a ROI_user object 
    """
    return Factory(modeldir,  dataspec,  **kwargs)
