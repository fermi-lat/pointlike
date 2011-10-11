"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.6 2011/10/03 22:04:11 burnett Exp $

"""

import numpy as np
from . import roistat, tools, printing, roisetup
from . import plotting 
from uw.utilities import fitter

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

class ROI_user(roistat.ROIstat):
    """ subclass of ROIstat that adds user-interface tools: fitting, localization, plotting sources, counts
    """

    def fit(self, select=None, **kwargs):
        """ Perform fit, return fitter object to examine errors, or refit
        parameters:
            select : list type of int or None
                if not None, the list is a subset of the parameter numbers to select
                to define a projected function to fit
        kwargs :
            ignore_exceptions : bool
                if set, run the fit in a try block and return None
            call_limit : int
                if set, modify default limit on number of calls
            others passed to the fitter minimizer command. defaults are
                estimate_errors = True
                use_gradient = True
        
        notes:
            if select is None, it will set the cov_matrix attributes of the corresponding Model
            objects
        """
        ignore_exception = kwargs.pop('ignore_exception', False)
        self.call_limit = kwargs.pop('call_limit', self.call_limit)
        fit_kw = dict(use_gradient=True, estimate_errors=True)
        fit_kw.update(kwargs)
        self.update()
        initial_value, self.calls = self.log_like(), 0
        fn = self if select is None else fitter.Projector(self, select=select)
        try:
            mm = fitter.Minimizer(fn)
            mm(**fit_kw)
            print '%d calls, likelihood improvement: %.1f'\
                % (self.calls, self.log_like() - initial_value)
            if fit_kw['estimate_errors'] and select is None:
                # tricky to do if fitting subset, punt for now
                self.sources.set_covariance_matrix(mm.cov_matrix)
            return mm
        except Exception:
            if ignore_exception: return None
            else: raise
        
    def localize(self, source_name, **kwargs):
        """ localize the source: adding the localization info to the source object
        If already done, just return the info (the tools.Localization object)
        """
        source = self.sources.find_source(source_name)
        update = kwargs.pop('update', False)
        if hasattr(source, 'loc') and not update:
            return source.loc
        loc= tools.Localization(self, source_name, **kwargs)
        try: 
            loc.localize()
        except Exception, e:
            print 'Failed localization for source %s: %s' % (source.name, e)
        source.loc = loc
        return loc
    
    def get_sources(self):
        return [ s for s in self.sources if s.skydir is not None]
    def get_model(self, name):
        return self.get_source(name).spectral_model
    def get_source(self, name):
        return self.sources.find_source(name)
        
    def summary(self, out=None, title=None):
        """ summary table of free parameters, values uncertainties if any"""
        if title is not None:
            print >>out, title
        print >>out, '%-20s%10s%10s'% tuple('Name value error(%)'.split())
        prev=''
        for (name, value, rsig) in zip(self.parameter_names, self.parameters, self.sources.uncertainties):
            sname,pname = name.split('_',1)
            if sname==prev: name = len(sname)*' '+'_'+pname
            prev = sname
            print >>out, '%-20s%10.4g%10.1f' %(name,value,rsig*100)
            
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
        
    def get_sed(self, source_name, **kwargs):
        """ return the SED recarray for the source
        """
        source = self.sources.find_source(source_name)
        update = kwargs.pop('update', False)
        if hasattr(source, 'sed_rec') and not update:
            return source.sed_rec
        sf =tools.SourceFlux(self, source_name, **kwargs)
        source.sed_rec = tools.SED(sf).rec
        return source.sed_rec

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
        self.get_sed(source_name)
        ps = plotting.sed.Plot(source)
        ps(**kwargs)
        return ps

    @decorate_with(plotting.counts.stacked_plots)
    def plot_counts(self, **kwargs):
        return plotting.counts.stacked_plots(self, **kwargs)
        
    @decorate_with(printing.print_summary)
    def print_summary(self, **kwargs):
        return printing.print_summary(self, **kwargs)

    def tsmap(self, source_name, **kwargs):
        """ provided for the like2.processor
            return function of likelihood in neighborhood of given source
            tsm = roi.tsmap(which)
            size=0.25
            tsp = image.TSplot(tsm, center, size, pixelsize =size/20, axes=plt.gca())
            tsp.plot(center, label=name)
            tsp.show()
        """
        loc = self.localize(source_name,**kwargs)
        return PySkyFunction(loc)
        
class Factory(roisetup.ROIfactory):
    def __call__(self, *pars, **kwargs):
        return ROI_user(super(Factory,self).__call__(*pars, **kwargs))

@decorate_with(roisetup.ROIfactory, append_init=True)
def factory(dataname='P7_V4_SOURCE_4bpd', indir='uw26', 
        analysis_kw={'minROI': 7, 'maxROI': 7, 'irf': 'P7SOURCE_V6', 'emin':100,},):
    """ will then return a ROI_user object """
    return Factory(dataname, indir, analysis_kw=analysis_kw)
