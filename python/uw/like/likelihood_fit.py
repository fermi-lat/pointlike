"""A module providing functionality for parametrizing a likelihood curve
by a simple function.
Classes:
    LogLikelihood: a representation of a likelihood curve
    
Authors: Eric Wallace, Matthew Kerr
"""

__version__ = "$Revision: 1.1 $"
#$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like/likelihood_fit.py,v 1.1 2011/09/27 19:54:39 wallacee Exp $

import numpy as np
import scipy.optimize as opt


class LogLikelihood(object):
    """A representation of a log likelihood curve by a Poisson-like function
    
    The representation used here follows the approach of Nolan, et al. The
    likelihood is represented by a three-parameter function of a form similar
    to the Poisson distribution PMF. The peak is found by maximizing the
    provided log likelihood function. The parameters are then found by a least
    squares fit using the peak, two points equispaced around it, and zero.
    The parametrizing function is
       f(s) = logL - logL_max = n*np.log(e*(s+b)) - e*(s+b) - n*np.log(n) + n
    with n = e*(s_peak+b)
    """

    def __init__(self,loglike,initial_value=1e-10,fit_max = True,pars=None):
        """Create a LogLikelihood instance.
        loglike: The log likelihood function to be parametrized. Should be a
            callable with one argument.
        initial_value: An initial guess at the maximum of the provided
            function. The default of 1e-10 should be a reasonable guess for
            the normalization of a PowerLaw model.
        fit_max: Whether to use fmin to maximize the log likelihood. If False
            initial_value will be taken as the position of the maximum of the
            log likelihood, so this should only be set to False if the value
            passed as initial_value is the result of a previous  maximization
            of the provided function.
        pars: A length three sequence providing the values for the parameters
            of the fit function: s_peak,e, and b. If provided, these values
            will be used and the loglike argument will be ignored.
        
        """
        self.function = self._setup_function(loglike)
        self.saved_points = np.array([])
        self.saved_values = np.array([])
        if pars is not None:
            try:
                assert(hasattr(pars,'__iter__') and len(pars)==3)
                self.pars = pars
            except AssertionError:
                print('Keyword argument pars must be a sequence of length 3.')
                print('Will attempt to derive parameters from provided function')
                self.pars = self._find_pars(initial_value,fit_max = fit_max)
        else:
            self.pars = self._find_pars(initial_value,fit_max = fit_max)
        self._check_agreement()

    def _setup_function(self,function):
        """Setup caching of values passed to the log likelihood function."""
        def _function(x):
            if x in self.saved_points:
                ll = self.saved_values[self.saved_points==x][0]
            else:
                ll = function(x) 
                self.saved_points = np.append(self.saved_points,x)
                self.saved_values = np.append(self.saved_values,ll)
            return ll
        return _function

    def _find_pars(self,initial_value,fit_max = False):
        """Find the best fit parameters for the fit function"""
        if fit_max:
            self.mode = opt.fmin(lambda x: -self.function(x),initial_value)[0]
        else:
            self.mode = initial_value
        self.max = self.function(self.mode)
        #xs = np.array([0,max/2,max,max*2])
        #ys = np.array([self.function(x) for x in xs]) 
        xs = self.saved_points.copy()
        ys = self.saved_values.copy()
        ys = ys - ys.max()
        return opt.leastsq(lambda x:self._poisson(x,xs)-ys,np.array([self.mode,10/self.mode,xs[-1]]),maxfev=5000)[0]

    def _poisson(self,pars,s):
        """Calculate the value of the parametrizing function for some parameters.
        pars: A sequence of length 3 providing the parameters s_peak, e, and b.
        s: The point at which to evaluate the function. Can be a numpy array.
        """
        
        if pars[0]<0: return -1e10
        s_peak,e,b = pars[0],pars[1],pars[2];n = e*(s_peak+b)
        #logL - logL_max = n*np.log(e*(s+b))-e*(s+b) - n*np.log(e*(s_peak+b))+e*(s_peak+b)
        #simplified:
        return n*np.log((s+b)/(s_peak+b)) + e*(s_peak-s)

    def __call__(self,x):
        """Return the value of the parametrizing function at point x."""
        return self._poisson(self.pars,x) + self.max

    def find_logl_change(self,initial_value,delta_logl):
        """Find the points where the likelihood has decreased by delta_logl.

           Returns a tuple of the (low, high) values. If the likelihood at zero
           differs from the max by less than the specified change, return zero
           for the lower value.
        """
        #First, find lower value
        lo = 1e-20 #basically zero
        hi = initial_value
        ll_0 = self.function(hi)
        if ll_0-self.function(lo)>delta_logl:
            for i in range(20):
                avg = .5*(hi+lo)
                ll = self.function(avg)
                if ll_0-ll<delta_logl: hi = avg
                else: lo = avg
                if abs(ll_0-ll-delta_logl)<.01: break
            lo_val = avg
        else: lo_val = lo
        #Now the upper value
        lo = initial_value
        hi = initial_value*10
        while ll_0-self.function(hi)<delta_logl: hi+=1
        for i in range(20):
            avg = .5*(lo+hi)
            ll = self.function(avg)
            if ll_0-ll<delta_logl: lo = avg
            else: hi = avg
            if abs(ll_0-ll-delta_logl)<.01: break
        hi_val = avg
        return (lo_val,hi_val)

    def _check_agreement(self):
        lo,hi = self.find_logl_change(self.mode,10)
        lo_ll,hi_ll = self.function(lo),self.function(hi)
        lo_val,hi_val = self(lo),self(hi)
        if abs(1-lo_ll/lo_val) > .05:
            print("Warning: fit function differs from log likelihood by {0:.02}\% in the low tail".format((1-lo_ll/lo_val)*100))
        if abs(1-hi_ll/hi_val) > .05:
            print("Warning: fit function differs from log likelihood by {0:.02}\% in the high tail".format((1-lo_ll/lo_val)*100))

    def ts(self):
        return self(self.mode)-self(0)
