"""
Top-level code for ROI analysis

$Header$

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
        if fit_kw['estimate_errors'] and select is not None:
            # tricky to do if fitting subset, punt for now
            self.sources.set_covariance_matrix(mm.cov_matrix)
        return mm
        
    def localize(self, **kwargs):
        return tools.Localize(self, **kwargs)
        
    def pardump(self, eps=1e-1):
        print 'variable parameters'
        par = self.get_parameters()
        grad = self.gradient(par)
        f0 = self(par)
        for i, name in enumerate(self.parameter_names):
            dpar = par.copy()
            dpar[i] += eps/grad[i]
            print '%-20s%10.2e%10.2e %10.2e'% (name, par[i], grad[i], (self(dpar)-f0)*grad[i]/eps)
            self.set_parameters(par)
            
    def summary(self, out=None):
        """ summary table of free parameters"""
        print >>out, '%-20s%10s'% tuple('Name value'.split())
        for (name, value) in zip(self.parameter_names, self.parameters):
            print >>out, '%-20s %10.4g' %(name,value)
            
    