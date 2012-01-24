"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.14 2011/12/21 16:05:59 burnett Exp $

"""
import types
import numpy as np
from . import roistat, localization, printing, roisetup, sedfuns, sources, plotting
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
    """

    def fit(self, select=None, exclude=None,  summarize=True, quiet=True, **kwargs):
        """ Perform fit, return fitter object to examine errors, or refit
        
        Parameters
        ----------
        select : None, item or list of items, where item is an int or a string
            if not None, it defines a subset of the parameter numbers to select
                to define a projected function to fit
            int:  select the corresponding parameter number
            string: select parameters according to maching rules
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
        npars = len(saved_pars)
        
        # process select
        selected= set() 
        if select is not None:
            selectpar = select
            if not hasattr(select, '__iter__'): select = [select]
            for item in select:
                if type(item)==types.IntType:
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
                    raise Exception('fit parameter select list item %s must be either an integer or a string' %item)
            select = sorted(list(selected))
            if len(select)==0:
                raise Exception('nothing selected for fit from selection "%s"' % selectpar)
        
        if exclude is not None:
            if not hasattr(exclude, '__iter__'): exclude = [exclude]
            all = set(range(npars)) if select is None else set(select)
            select = list( all.difference(exclude))
                
        fn = self if select is None  else fitter.Projector(self, select=select)
        try:
            mm = fitter.Minimizer(fn, quiet=quiet)
            mm(**fit_kw)
            print '%d calls, likelihood improvement: %.1f'\
                % (self.calls, self.log_like() - initial_value)
            if fit_kw['estimate_errors'] :
                self.sources.set_covariance_matrix(mm.cov_matrix, select)
            if summarize: self.summary(select)
            return mm
        except FloatingPointError, e:
            print 'Fit error: restoring parameters since  "%s"' %e 
            self.set_parameters(saved_pars)
            return mm
        
        except Exception:
            if ignore_exception: return None
            else: raise
    
    def cov_matrix(self, par=None):
        return fitter.Minimizer.mycov(self, self.get_parameters() if par is None else par)
        
    def localize(self, source_name, **kwargs):
        """ localize the source, return elliptical parameters 
        """
        source = self.sources.find_source(source_name)
        with localization.Localization(self, source_name, **kwargs) as loc:
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
    def get_model(self, name):
        return self.get_source(name).spectral_model
    def get_source(self, name):
        return self.sources.find_source(name)
        
    def summary(self, select=None, out=None, title=None, gradient=True):
        """ summary table of free parameters, values uncertainties gradient
        
        Parameters:
        ----------
        select : list of integers or string
            integers are indices of parameters
            string is the wildcarded name of a source
        out : open file or None
            
        """
        if title is not None:
            print >>out, title
        fmt, tup = '%-21s%6s%10s%10s', tuple('Name index value error(%)'.split())
        if gradient:
            grad = self.gradient()
            fmt +='%10s'; tup += ('gradient',)
        print >>out, fmt %tup
        prev=''
        if type(select)==types.StringType:
            src = self.get_source(select)
            select = [i for i in range(len(saved_pars)) if self.parameter_names[i].startswith(src.name)]

        for index, (name, value, rsig) in enumerate(zip(self.parameter_names, self.parameters, self.sources.uncertainties)):
            if select is not None and index not in select: continue
            t = name.split('_')
            pname = t[-1]
            sname = '_'.join(t[:-1])
            #sname,pname = name.split('_',1)
            if sname==prev: name = len(sname)*' '+'_'+pname
            prev = sname
            fmt = '%-21s%6d%10.4g%10.1f'
            tup = (name,index, value,rsig*100)
            if gradient:
                fmt +='%10.1f'; tup += (grad[index],)
            print >>out,  fmt % tup
            
    def TS(self, which, quick=True):
        """ measure the TS for the given source
        """
        source_name=which
        source = self.sources.find_source(source_name)
        model = source.spectral_model
        norm =model[0]
        model[0]=1e-15
        self.update()
        llzero = self.log_like()
        model[0]=norm; self.update()
        ts= 2*(self.log_like()-llzero)
        return max(ts, 0)

    def band_ts(self, source_name):
        sed = self.get_sed(source_name)
        return np.sum(sed.ts)
        
    def get_sed(self, source_name, event_class=None, **kwargs):
        """ return the SED recarray for the source
            event_class : None, or integer
        """
        source = self.sources.find_source(source_name)
        update = kwargs.pop('update', False)
        if hasattr(source, 'sedrec') and not update:
            return source.sedrec
        sf = sedfuns.SourceFlux(self, source_name, **kwargs)
        source.sedrec = sedfuns.SED(sf, event_class=event_class).rec
        return source.sedrec

    @decorate_with(plotting.tsmap.plot)
    def plot_tsmap(self, source_name, **kwargs):
        """ create a TS map showing the source localization
        """
        plot_kw = dict(size=0.25, pixelsize=0.25/15, outdir='plots')
        plot_kw.update(kwargs)
        loc= self.localize( source_name)
        tsp = plotting.tsmap.plot(loc, **plot_kw)
        return tsp
        
    @decorate_with(plotting.sed.Plot, append_init=True)    
    def plot_sed(self, source_name, **kwargs):
        source = self.sources.find_source(source_name)
        source.sedrec = self.get_sed(source_name, 
            event_class=kwargs.pop('event_class', None), update=kwargs.pop('update',True))
        ps = plotting.sed.Plot(source)
        ps(**kwargs)
        return ps

    @decorate_with(plotting.counts.stacked_plots)
    def plot_counts(self, **kwargs):
        return plotting.counts.stacked_plots(self, **kwargs)
        
    @decorate_with(printing.print_summary)
    def print_summary(self, **kwargs):
        return printing.print_summary(self, **kwargs)

    @property
    def bounds(self):
        return self.sources.bounds
    @property
    def cov_matrix(self):
        """ the current covariance matrix, determined from the gradient """
        return fitter.Minimizer.mycov(self.gradient, self.get_parameters())
        
    @property
    def correlations(self):
        """Return the linear correlation coefficients for the estimated covariance matrix."""
        cm = self.cov_matrix
        s = np.sqrt(cm.diagonal())
        return cm / np.outer(s,s)
        
    def add_source(self, **kwargs):

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
        
    def zero_source(self, source, **kwargs):
        raise Exception('not implemented')
    def unzero_source(self, source, **kwargs):
        raise Exception('not implemented')
        
class Factory(roisetup.ROIfactory):
    def __call__(self, *pars, **kwargs):
        return ROI_user(super(Factory,self).__call__(*pars, **kwargs))

@decorate_with(roisetup.ROIfactory, append_init=True)
def factory(dataname='P7_V4_SOURCE_4bpd', indir='uw26', 
        analysis_kw={'minROI': 7, 'maxROI': 7, 'irf': 'P7SOURCE_V6', 'emin':100,},):
    """ will then return a ROI_user object """
    return Factory(dataname, indir, analysis_kw=analysis_kw)
