"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.2 2011/09/05 18:59:26 burnett Exp $

"""

import numpy as np
from . import roistat, tools, plotting
from uw.utilities import fitter


class Main(roistat.ROIstat):
    """ subclass of ROIstat that adds user-interface tools
    """

    def fit(self, select=None, **kwargs):
        """ Perform fit, return fitter object to examine errors, or refit
        parameters:
            select : list type of int or None
                if not None, the list is a subset of the paramter numbers to select
                to define a projected function to fit
        kwargs :
            passed to the fitter minimizer command. defaults are
                estimate_errors = True
                use_gradient = True
        
        notes:
            if select is None, it will set the cov_matrix attributes of the corresponding Model
            objects
        """
        fit_kw = dict(use_gradient=True, estimate_errors=True)
        fit_kw.update(kwargs)
        initial_value, self.calls = self.log_like(), 0
        fn = self if select is None else fitter.Projector(self, select=select)
        mm = fitter.Minimizer(fn)
        mm(**fit_kw)
        print '%d calls, likelihood improvement: %.1f'\
            % (self.calls, self.log_like() - initial_value)
        if fit_kw['estimate_errors'] and select is None:
            # tricky to do if fitting subset, punt for now
            self.sources.set_covariance_matrix(mm.cov_matrix)
        return mm
        
    def localize(self, source_name, **kwargs):
        """ localize the source: adding the localization info to the source object
        If already done, just return the info (the tools.Localization object)
        """
        source = self.sources.find_source(source_name)
        update = kwargs.pop('update', False)
        if hasattr(source, 'loc') and not update:
            return source.loc
        loc= tools.Localization(self, source_name, **kwargs)
        loc.localize()
        source.loc = loc 
        return loc
        
    def summary(self, out=None):
        """ summary table of free parameters, values"""
        print >>out, '%-20s%10s'% tuple('Name value'.split())
        prev=''
        for (name, value) in zip(self.parameter_names, self.parameters):
            sname,pname = name.split('_',1)
            if sname==prev: name = len(sname)*' '+'_'+pname
            prev = sname
            print >>out, '%-20s%10.4g' %(name,value)
            
    def TS(self, source_name, quick=True):
        """ measure the TS for the given source
        """
        source = self.sources.find_source(source_name)
        model = source.spectral_model
        norm =model[0]
        model[0]=1e-15
        self.update()
        llzero = self.log_like()
        model[0]=norm; self.update()
        ts= 2*(self.log_like()-llzero)
        return ts

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

    def plot_tsmap(self, source_name, **kwargs):
        """ create a TS map showing the source localization
        """
        plot_kw = dict(size=0.25, pixelsize=0.25/15, outdir='plots')
        plot_kw.update(kwargs)
        loc= self.localize( source_name)
        tsp = plotting.tsmap.plot(loc, **plot_kw)
        return tsp
        
    def plot_sed(self, source_name, **kwargs):
        source = self.sources.find_source(source_name)
        self.get_sed(source_name)
        ps = plotting.sed.Plot(source)
        ps(**kwargs)
        return ps
