"""
Top-level code for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/main.py,v 1.1 2011/09/02 16:30:43 burnett Exp $

"""

import numpy as np
from . import roistat, tools
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
        """ for now, return a Localization object, which see """
        return tools.Localization(self, source_name, **kwargs)
        
           
    def summary(self, out=None):
        """ summary table of free parameters"""
        print >>out, '%-20s%10s'% tuple('Name value'.split())
        for (name, value) in zip(self.parameter_names, self.parameters):
            print >>out, '%-20s %10.4g' %(name,value)
            
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

    def source_flux(self, source_name, **kwargs):
        """ return a SourceFlux object for the source, which see"""
        return tools.SourceFlux(self, source_name, **kwargs)