"""
analyze time variation
$Header$

"""

import os,sys, glob,pickle
import numpy as np
import pylab as plt
import scipy
from skymaps import BinnedPhotonData, SkyDir,LivetimeCube
from uw.like import pointspec_helpers,pixeldata
from uw.thb_roi import config

from uw import factory

class MonthlyData(object):
    def __init__(self,  datapath = 'monthly', **kwargs):
        self.bpdfiles = np.sort(glob.glob(config.data_join(datapath,'bpd','month_*_4bpd.fits')))
        self.ltfiles  = np.sort(glob.glob(config.data_join(datapath,'lt','month_*.fits')))
        self.n = len(self.bpdfiles)
        assert self.n>0, 'Did not find monthly data files'
        print 'Found %d monthly data files' % self.n
        self.allbpd = np.sort(glob.glob(config.data_join(datapath,'*months_4bpd.fits')))
        self.alllt = np.sort(glob.glob(config.data_join(datapath, '*months_lt.fits')))
        
    def __call__(self, n):
        if n==-1:
            return dict(number=-1, dataset=(self.allbpd[0], self.alllt[0]))
        return dict(number=n, dataset =(self.bpdfiles[n] , self.ltfiles[n])) 


   
class LogLikelihood(object):
    """ manage the 1-dimensional likelihood function for a pointlike source"""

    def __init__(self, roi, maxflux=None):
        """ setup a function to evaluate the likelihood function for a source, 
            varying only the differential flux
        roi: an ROIAnalysis object     
        """
        self.roi = roi
        self.model = roi.psm.models[0]
        self.loglike = roi.logLikelihood
        
        if maxflux is None:
            # find the maximum by examining the function
            self.maxl = self.maximum()
            self.wpeak = self.eval(self.maxl)
            self.maxflux = self.find_maxflux()
        else:
            # user specified: will make sure is ok, increase if not
            self.maxflux = maxflux
        self.tot=1.0
        self.tot = self.integral( self.maxflux)

    def __str__(self):
        h = tuple('max -sig +sig 95% TS'.split())
        t =  (self.maxl,) + self.errors() +(self.upper_limit(),  self.TS())
        n = len(t)
        return '\n'.join([(n*'  %-9s') % h, (n*' %10.2e') % t])
        
    def eval(self, norm):
        """ evaluate the log likelihood function"""
        def f(x):
            self.model.p[0] = np.log10(x) if x>1e-20 else -20
            return self.loglike(self.roi.parameters())
        
        if hasattr(norm, '__iter__'):
            return np.array( map(f,norm) )
        return f(norm)
        
    def __call__(self, x):
        """ evalate the likelihood function, normalized to 1 at peak"""
        return np.exp(-self.eval(x)+self.wpeak)

    def maximum(self, val=1e-12, disp=0):
        """ find the position of the maximum likelihood (for setup: ignore spline)
        val : starting value for fmin
        disp: pass to fmin
        
        note do it in log(flux) space
        """
        return np.exp(scipy.optimize.fmin( lambda x: self.eval(np.exp(x)), np.log(val), disp=disp)[0])
        
    def find_maxflux(self, tol=10):
        """ find the maximum value at which the log likelihood difference is small
        """
        mval = self.maximum()
        if mval<1e-14: mval=1e-14
        wmax = self.eval(mval)
        v = mval*10
        for i in range(10):
            if self.eval(v)> wmax+tol: break
            v*=np.e
        return v
            
         
    def errors(self):
        """ tuple with lower, upper 1-sigma uncertainty flux values"""
        delta = self.wpeak + np.exp(-0.5)
        g = lambda x: self.eval(x)-delta # function that should be zero at the 1 sigma points
        xtol = 1e-3*self.maxflux
        yl=yu=0
        if self.maxl==0: return (yl,yu)
        try:
            yl= scipy.optimize.brentq( g, 0,         self.maxl,    xtol=xtol)
        except: pass
        try:
            yu= scipy.optimize.brentq( g, self.maxl, self.maxflux, xtol=xtol)
        except: 
            raise
        return (yl,yu) 
            
    def integral(self, x):
        """ the integral of the normalized function"""
        f =  lambda x : scipy.integrate.quad(self, 0, x)[0]/self.tot
        if hasattr(x, '__iter__'):
            return np.array( map(f,x) )
        return f(x)
        
    
    def upper_limit(self, cl = 0.95):
        """ the flux value at the confidence level"""
        if cl>1: cl /= 100.
        t =0
        try:
            t = scipy.optimize.brentq( lambda x: self.integral(x) - cl, 0, self.maxflux, xtol=1e-3*self.maxflux)
        except: pass
        return t
     
    def TS(self):
        return 2*(self.eval(0)- self.eval(self.maxl))
        
    def plot(self,fignum=10, axes = None, **kwargs):
        """ a simple plot of the likelihood """
        if axes is None:
            plt.figure(fignum); plt.clf()
            axes = plt.gca()
        dom = np.linspace(0, self.maxflux, 101)
        axes.plot(dom, self(dom),'-', **kwargs)
        axes.grid(True)
        axes.set_ylim((-0.1, 1.1))

    

