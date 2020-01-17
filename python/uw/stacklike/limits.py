#!/usr/bin/env python
# $Header: /nfs/slac/g/glast/ground/cvs/users/mdwood/python/DMLimits.py,v 1.6 2011/03/31 17:02:25 kadrlica Exp $

"""
@author Matthew Wood <mdwood@slac.stanford.edu>
"""

__author__   = "Matthew Wood"
__date__     = "12/20/2010"

import sys
import scipy as sp
import numpy as np
import math
import scipy.optimize as opt
import scipy.interpolate as intp
import scipy.special as spf
from scipy.integrate import quad
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.stats as stats
import copy

def gauss(x,mu,sigma=1.0):
    s2 = sigma*sigma
    return 1./np.sqrt(2*s2*np.pi)*np.exp(-(x-mu)*(x-mu)/(2*s2))

def lngauss(x,mu,sigma=1.0):
    s2 = sigma*sigma

    return -0.5*np.log(2*s2*np.pi) - np.power(x-mu,2)/(2*s2)

def lgauss(x,mu,sigma=1.0):
    
    x = np.asarray(x)

    lmu = np.log10(mu)
    s2 = sigma*sigma

    # Handle scalar case
    if x.ndim == 0:
        if x <= 0:
            return -1000
        else:
            lx = np.log10(x)
            return (1./np.sqrt(2*s2*np.pi)*
                    np.exp(-np.power(lx-lmu,2)/(2*s2))/(x*np.log(10.)))
 
    lx = np.zeros(x.shape)
    v = np.zeros(x.shape)

    lx[x>0] = np.log10(x[x>0])

    v = (1./np.sqrt(2*s2*np.pi)*
         np.exp(-(lx-lmu)*(lx-lmu)/(2*s2))/(x*np.log(10.)))
    v[x<=0] = -1000

    return v

def lnlgauss(x,mu,sigma=1.0):

    x = np.asarray(x)

    lmu = np.log10(mu)
    s2 = sigma*sigma

    # Handle scalar case
    if x.ndim == 0:
        if x <= 0:
            return -1000
        else:
            lx = np.log10(x)
            return -0.5*np.log(2*s2*np.pi) - \
                np.power(lx-lmu,2)/(2*s2) - 2.302585*lx - np.log(np.log(10.))


    lx = np.zeros(x.shape)
    v = np.zeros(x.shape)

    lx[x>0] = np.log10(x[x>0])

    v = -0.5*np.log(2*s2*np.pi) - \
        np.power(lx-lmu,2)/(2*s2) - 2.302585*lx - np.log(np.log(10.))
    v[x<=0] = -1000

    return v


class LnLFn(object):
    """
    Helper class for interpolating a 1-D log-likelihood function from a
    set of tabulated values.  
    """
    def __init__(self,x,y):
        self._xmin = x[0]
        self._xmax = x[-1]
        self._ymin = y[0]
        self._ymax = y[-1]
        self._dydx_lo = (y[1]-y[0])/(x[1]-x[0])
        self._dydx_hi = (y[-1]-y[-2])/(x[-1]-x[-2])

        self._fn = intp.UnivariateSpline(copy.copy(x),copy.copy(y),s=0)

    def __call__(self,x):

        xnew = np.asarray(x)
        nx = xnew.ndim

        below_bounds = xnew < self._xmin
        above_bounds = xnew > self._xmax

        # Handle case where x is a scalar
        if xnew.ndim == 0:

            if xnew > self._xmax:
                return self._ymax + (x-self._xmax)*self._dydx_hi
            elif xnew < self._xmin:
                return self._ymin + (x-self._xmax)*self._dydx_lo
            else:
                return self._fn(xnew)
        else:

            dxhi = np.asarray(xnew-self._xmax)
            dxlo = np.asarray(xnew-self._xmin)

            # UnivariateSpline will only accept 1-D arrays so this
            # passes a flattened version of the array.
            y = self._fn(xnew.ravel())
            y.resize(xnew.shape)
            
            # If y is a rank-0 array convert it to rank-1
            if y.ndim == 0:
                y = y[np.newaxis,...]

            y[above_bounds] = (self._ymax + dxhi[above_bounds]*self._dydx_hi)
            y[below_bounds] = (self._ymin + dxlo[below_bounds]*self._dydx_lo)
            return y




class MarginalLnL(object):
    """
    MarginalLnL(lnlx,lnly,sigma=0.1,fntype='lgauss')

    Generate the marginal likelihood for the parameter x given a
    likelihood L(z) where z = x*y and y is a nuisance parameter.
    The marginalization is calculated by evaluating the 1-D integral:

    L(x) = \int L(x*y)P(y)dy

    where P(y) is the p.d.f. for y which is assumed to have the
    physical domain [0,inf].  The 1-D log-likelihood for z evaluated
    at y=y_{0} is passed as a set of tabulated values using the arrays
    `lnlx` and `lnly`.  

    This class returns a function that can be used to evaluate the
    marginal likelihood.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.

    fntype : string 
  
       Function with which to model the nuisance parameter.  Options
       are normal ('gauss') and log-normal ('lgauss').  In the case of
       the normal distribution the p.d.f. is truncated at zero.

    sigma : float
        Width parameter of the nuisance parameter distribution.

    Examples
    --------
    >>> import DMLimits as dmlim

    >>> lnlx = np.linspace(0,10,100)
    >>> lnly = lnlx*lnlx
    >>> mlnl = dmlim.MarginalLnL(lnlx,lnly,sigma=0.2,fntype='lgauss')
    >>> mlnly = mlnl(lnlx)


    """

    def __init__(self,lnlx,lnly,sigma=0.1,fntype='lgauss'):
        self._sigma=sigma
        self._pdfnorm=1

        if len(lnlx.shape) != 1 or len(lnly.shape) != 1:
            raise ValueError("")
        elif len(lnlx) < 2 or len(lnly) < 2:
            raise ValueError("")

        if fntype == 'lgauss':
            self._ypdf = lgauss
            self._ylnpdf = lnlgauss
        elif fntype == 'gauss':
            self._ypdf = gauss
            self._ylnpdf = lngauss
            self._pdfnorm = 1-0.5*(1+spf.erf(-1.0/(math.sqrt(2.0)*self._sigma)))
        else:
            raise ValueError("fntype = %s is not supported." %(fntype))

        self._fn_lnl = LnLFn(lnlx,lnly)

    def nuisancePDF(self,x):

        return self._ypdf(x,1.0,self._sigma)/self._pdfnorm

    def like(self,x,y):
        """Evaluate the 2-D likelihood in the x/y parameter space.

        Parameters
        ----------
        x : array_like
        Array of coordinates in the `x` parameter.
        
        y : array_like       
        Array of coordinates in the `y` nuisance parameter.
        """

        z = self._fn_lnl(x*y)
        return np.exp(z)*self._ypdf(y,1.0,self._sigma)/self._pdfnorm

    def lnlike(self,x,y):

        return self._fn_lnl(x*y) + \
            self._ylnpdf(y,1.0,self._sigma) - np.log(self._pdfnorm)

    def __call__(self,x):
        """Evaluate the marginal log-likelihood."""

        x = np.asarray(x)

        if x.ndim == 0:
            return np.log(quad(lambda t: self.like(x,t),0.0,np.inf)[0])
        else:            
            z = []

            for xtmp in x:
                z.append(quad(lambda t: self.like(xtmp,t),0.0,np.inf)[0])

            return np.log(z)
        
class ProfileLnL(MarginalLnL):

    def get_xmax(self,x):

        x = np.asarray(x)

        if x.ndim == 0:
            
            fn = lambda t: -self.lnlike(x,t)
            ytmp = opt.fmin(fn,1.0,disp=False)[0]
            return ytmp

    def __call__(self,x):
        """Evaluate the profiled log-likelihood."""

        x = np.asarray(x)

        if x.ndim == 0:

            fn = lambda t: -self.lnlike(x,t)
            ytmp = opt.fmin(fn,1.0,disp=False)[0]
            return self.lnlike(x,ytmp)
        else:            
            z = []
            y = []

            for xtmp in x:

                fn = lambda t: -self.lnlike(xtmp,t)
                ytmp = opt.fmin(fn,1.0,disp=False)[0]

                ztmp = self.lnlike(xtmp,ytmp)

                z.append(ztmp)
                y.append(ytmp)

            return z


class BayesianLimit(object):
    """
    BayesianLimit(lnlx,lnly,npdf=1000)

    Evaluate the upper limit on a parameter using Bayesian
    methodology.  A flat prior on the interval [0,inf] is assumed.
    A 1-D log-likelihood curve is defined using the input arrays
    `lnlx` and `lnly`.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.

    Examples
    --------
    >>> import DMLimits as dmlim

    >>> lnlx = np.linspace(0,10,100)
    >>> lnly = lnlx*lnlx
    >>> blim = dmlim.BayesianLimit(lnlx,lnly)
    >>> print (blim.getLimit(0.05))


    """

    def __init__(self,lnlx,lnly,npdf=1000):
        self._fn_lpdf = LnLFn(lnlx,lnly)
        self._fn_pdf = lambda t: math.exp(self._fn_lpdf(t))

        newlx = np.linspace(min(lnlx),max(lnlx),30)
        dx = (newlx[1]-newlx[0])
        lnlx=newlx
        self.cpdf = [0]

        s = 0
        for i, xval in enumerate(newlx):

            if i==0: 
                continue

            s += quad(self._fn_pdf,xval-dx,xval)[0]
            self.cpdf.append(s)

        self.cpdf = np.array(self.cpdf)
        self.cpdf /= s

        vp = np.linspace(0.0,1.0,npdf)

        self.fn_cpdf = intp.UnivariateSpline(lnlx,self.cpdf,s=0)

        self.xvp = [0]

        xmin = lnlx[0]
        xmax = lnlx[1]

        # Invert the cumulative probability function
        for i, p in enumerate(vp):

            if i == 0:
                continue
            elif i+1 == len(vp):
                self.xvp.append(lnlx[-1])
            else:

                rf = lambda x: self.fn_cpdf(x)-p

                while rf(xmax) < 0:
                    xmax *= 1.1

                xp = opt.brentq(rf,xmin,xmax)
                #print (i, p, xmin, xmax, rf(xmin), rf(xmax), lnlx[-1], xp)
                self.xvp.append(xp)
                xmin = xp

        self._icdf = intp.UnivariateSpline(vp,self.xvp,s=0)

    def getLimit(self,alpha=0.05):
        """Evaluate the upper limit corresponding to a C.L. of (1-alpha)%.

        Parameters
        ----------
        alpha : Upper limit confidence level.
        """

        return self._icdf(1-alpha)

class ProfileLimit(object):
    """
    ProfileLimit(lnlx,lnly)

    Evaluate the upper limit on a parameter using the MINOS
    methodology.  A 1-D log-likelihood curve for the parameter is
    defined using the input arrays `lnlx` and `lnly`.

    Parameters
    ----------
    lnlx : array_like
       Array of points at which the log-likelihood is evaluated.
        
    lnly : array_like       
       Array of log-likelihood values evalated at the points in
       `lnlx`.

    Examples
    --------
    >>> import DMLimits as dmlim

    >>> lnlx = np.linspace(0,10,100)
    >>> lnly = lnlx*lnlx
    >>> plim = dmlim.ProfileLimit(lnlx,lnly)
    >>> print (plim.getLimit(0.05)  )
    """

    def __init__(self,lnlx,lnly):

        self._xmin = lnlx[0]
        self._xmax = lnlx[-1]

        self._fn = LnLFn(lnlx,lnly)
         
        self._xmin = opt.fminbound(lambda x: -self._fn(x),lnlx[0],lnlx[-1])
        self._fmax = self._fn(self._xmin)
        #print (self._xmin, self._fmax)

    def getLimit(self,alpha):
        """Evaluate the upper limit corresponding to a C.L. of (1-alpha)%.

        Parameters
        ----------
        alpha : Upper limit confidence level.
        """

        dlnl = pow(math.sqrt(2.)*spf.erfinv(1-2*alpha),2.)/2.  
        rf = lambda x: self._fn(x)+dlnl-self._fmax
        return opt.brentq(rf,self._xmin,self._xmax)

class StackedLnL(object):
    """
    StackLnL()
    
    Stack multiple profile log-likelihoods in order to calculate
    a combined log-likelihood curve.
    """
    def __init__(self):
        self.lnl_fn = []

    def addlnl(self,lnlx,lnly):
        fn = LnLFn(lnlx,lnly)
        self.lnl_fn.append(fn)

    def eval(self,x):
        s = 0
        for fn in self.lnl_fn:
            s += fn(x)
        return s

    def __call__(self,x):
        return self.eval(x)

if __name__ == '__main__':

    x0 = 1

    fntype = 'lgauss'
    jsigma = 0.32

    lnlx = np.linspace(0.0,20.0,200)
    lnly = -np.power(lnlx-x0,2)

    mlnl = MarginalLnL(lnlx,lnly,sigma=jsigma,fntype=fntype)
    plnl = ProfileLnL(lnlx,lnly,sigma=jsigma,fntype=fntype)

    mlnly = np.array(mlnl(lnlx))
    plnly = np.array(plnl(lnlx))

    mlnl_fn = LnLFn(lnlx,mlnly)
    plnl_fn = LnLFn(lnlx,plnly)

    mlnl_xmax = opt.fminbound(lambda x: -mlnl_fn(x),lnlx[0],lnlx[-1])
    plnl_xmax = opt.fminbound(lambda x: -plnl_fn(x),lnlx[0],lnlx[-1])

    mlnl_ymax = mlnl_fn(mlnl_xmax)
    plnl_ymax = plnl_fn(plnl_xmax)

    print ('mlnl xmax ', mlnl_xmax)
    print ('plnl xmax ', plnl_xmax)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylabel('-lnL')
    ax.grid(True)

    ax.set_ylim(0,10)
    ax.set_xlim(0,10)


    plt.plot(lnlx, -lnly, 'r-', linewidth=2)
    plt.plot(lnlx, -mlnly+mlnl_ymax, 'b--', linewidth=2,label='Marginal LnL')
    plt.plot(lnlx, -plnly+plnl_ymax, 'g-.', linewidth=2,label='Profile LnL')

    ax.legend()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(True)

    xlo = 0.1
    xhi = 15.
    ylo = 0.01
    yhi = 4

    t1,t2 = sp.mgrid[xlo:xhi:.01, ylo:yhi:.01]

    z = mlnl.like(t1,t2)
    lnz = mlnl.lnlike(t1,t2)

    x = np.linspace(xlo,xhi,100)
    y = [plnl.get_xmax(t) for t in x]

    plt.plot(x,y,linewidth=2,color='r',linestyle='--')

    cs = plt.contour(-lnz.transpose(), 
                     extent=[xlo,xhi,ylo,yhi],levels=[1.0,2.3,4.6,6.],
                     cmap=mpl.cm.jet, origin='lower',vmax=3)
    plt.clabel(cs, fontsize=9, inline=1)

    fig = plt.figure()

    y = np.linspace(0.1,3.0)


    z1 = [plnl.lnlike(2.0,t) for t in y]
    z2 = [plnl.lnlike(1.0,t) for t in y]
    z3 = [plnl.lnlike(0.5,t) for t in y]
    z4 = [plnl.lnlike(0.1,t) for t in y]

    plt.plot(y,z1)
    plt.plot(y,z2)
    plt.plot(y,z3)
    plt.plot(y,z4)


    print ('95% C.L. Profile Limit:                 ', )
    print (ProfileLimit(lnlx,lnly).getLimit(0.05))
    print ('95% C.L. Profile Limit (Marginalized):  ', )
    print (ProfileLimit(lnlx,mlnly).getLimit(0.05))
    print ('95% C.L. Bayesian Limit:                ', )
    print (BayesianLimit(lnlx,lnly).getLimit(0.05))
    print ('95% C.L. Bayesian Limit (Marginalized): ', )
    print (BayesianLimit(lnlx,mlnly).getLimit(0.05))

    plt.show()
