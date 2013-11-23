"""
Tools for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/tools.py,v 1.13 2012/06/24 04:52:29 burnett Exp $

"""
import os
import numpy as np
#from scipy import optimize, integrate
#from skymaps import SkyDir
#from uw.like import quadform, srcid
#from uw.utilities import keyword_options
#from . plotting import tsmap
#

def ufunc_decorator(f): # this adapts a bound function
    def new_ufunc(self, par):
        return np.array(map(lambda x: f(self,x),par)) if hasattr(par, '__iter__')  else f(self,par)
    return new_ufunc

    
class WithMixin(object):
    """Mixin to allow simple restore of an object's state
        supports the 'with' construction, guarantees that restore is called to restore the state of the model
        example:
        -------
        with ClassName(...) as something:
            # use something ...
    """
    def __enter__(self):
        return self
        
    def __exit__(self, type, value, traceback):
        self.restore()

    
      
#class LogLikelihood(object):
#    """ manage a 1-dimensional likelihood function """
#
#    def __init__(self, loglikelihood,  guess=1.0):
#        """ 
#        loglikelihood: a function of one parameter: expect to have a maximum
#        guess: a guess for setting the scale
#        """
#        self.loglike = loglikelihood
#        
#        self.wpeak=0.
#        self.maxl = self.maximum(guess)
#        self.wpeak = self.logL(self.maxl)
#        self.maxflux = self.find_maxflux(guess)
#        self.tot=None
#
#    def __str__(self):
#        h = tuple('max -sig +sig 95% TS'.split())
#        t =  (self.maxl,) + self.errors() +(self.upper_limit(),  self.TS())
#        n = len(t)
#        return '\n'.join([(n*'  %-9s') % h, ((n-1)*' %10.2e'+'%10.1f') % t])
#        
#    def logL(self, norm):
#        """ evaluate the log likelihood function
#        norm: a value, or array of values
#        """
#        if hasattr(norm, '__iter__'):
#            return np.array( map(self.loglike,norm) )
#        return self.loglike(norm)
#        
#    def __call__(self, x):
#        """ evalate the likelihood function, normalized to 1 at peak"""
#        return np.exp(self.logL(x)-self.wpeak)
#
#    def maximum(self, guess=1.0, disp=0):
#        """ find the position of the maximum likelihood
#        val : starting value for fmin
#        disp: pass to fmin
#        
#        """
#        ret = optimize.fmin( lambda x: -self.logL(x), guess, disp=disp)[0]
#        return ret if ret>0 else  0
#        
#    def find_maxflux(self, guess=1.0, tol=10):
#        """ find an appropriate maximum flux for evaluating the likelihood
#        """
#        wmax = self.logL(self.maxl)
#        v = self.maxl*2
#        if v>1e-15:
#            for i in range(10):
#                if self.logL(v) < wmax-tol: break
#                v*=np.e
#            return v
#        # case with just a limit: return flux where log likelihood is -5
#            
#        z = -self.logL(0)+5
#        g = lambda x: self.logL(x)+z
#        a = 0.1*guess
#        for i in range(10):
#            if g(a)>0 : break
#            a /= 10
#        b =10*guess
#        for i in range(10):
#            if g(b)<0: break
#            b *= 10
#        v =  optimize.brentq(g, a, b , xtol=1e-3*guess )
#        return v
#            
#    def errors(self, tol=None):
#        """ tuple with lower, upper 1-sigma uncertainty flux values (profile method)"""
#        delta = self.wpeak -0.5
#        g = lambda x: self.logL(x)-delta # function that should be zero at the 1 sigma points
#        xtol = 1e-3*self.maxflux if tol is None else tol
#        yl=yu=0
#        if self.maxl==0: return (yl,yu)
#        try:
#            yl= optimize.brentq( g, 0,         self.maxl,    xtol=xtol)
#        except: pass
#        assert g(self.maxl)*g(self.maxflux)<0, \
#            'bad: peak, max %.2e, %.2e, values: %.2e,%.2e' % \
#            (self.maxl,self.maxflux, self.logL(self.maxl),self.logL(self.maxflux))
#        yu= optimize.brentq( g, self.maxl, self.maxflux, xtol=xtol)
#        return (yl,yu) 
#            
#    def integral(self, x):
#        """ the integral of the normalized function"""
#        if self.tot is None:
#            self.tot=1.0
#            self.tot=self.integral(self.maxflux)
#        f =  lambda x : integrate.quad(self, 0, x, epsabs=1e-4,epsrel=1e-4)[0]/self.tot
#        if hasattr(x, '__iter__'):
#            return np.array( map(f,x) )
#        return f(x)
#        
#    
#    def upper_limit(self, cl = 0.95):
#        """ the flux value at the confidence level"""
#        if cl>1: cl /= 100.
#        t =0
#        try:
#            t = optimize.brentq( lambda x: self.integral(x) - cl, 0, self.maxflux, xtol=1e-3*self.maxflux)
#        except: pass
#        return t
#     
#    def TS(self):
#        return 2*( self.logL(self.maxl) - self.logL(0) )
#        
#    def plot(self,fignum=10, axes = None, **kwargs):
#        """ a simple plot of the likelihood """
#        import pylab as plt
#        if axes is None:
#            plt.close(fignum)
#            plt.figure(fignum)
#            axes = plt.gca()
#        dom = np.linspace(0, self.maxflux, 101)
#        if 'lw' not in kwargs or 'linewidth' not in kwargs: kwargs.update(lw=2)
#        axes.plot(dom, self(dom),'-', **kwargs)
#        color = kwargs.get('color', 'b')
#        a,b = self.errors()
#        if b==0: b=self.upper_limit()
#        axes.axvspan(a,b, color=color, alpha=0.25)
#        axes.grid(True)
#        axes.set_ylim((-0.1, 1.1))
#    @staticmethod
#    def test(mean=1.0, sigma = 0.1, guess=1.0):
#        def fun(x):
#            return -0.5*((x-mean)/sigma)**2
#        t =  LogLikelihood(fun, guess=guess)
#        #print t
#        return t
  