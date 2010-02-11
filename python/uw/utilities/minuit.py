"""
Provides a  convenience class to call Minuit, mimicking the interface to scipy.optimizers.fmin.

author: Eric Wallace <wallacee@uw.edu>
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/minuit.py,v 1.7 2010/02/01 23:53:41 burnett Exp $

"""
import sys, os
# normal CMT setup does not put ROOT.py in the python path
if sys.platform == 'win32':
    import win32api
    console_title = win32api.GetConsoleTitle()
try:
    import ROOT
except:
    sys.path.append(os.path.join(os.environ['ROOTSYS'], 'bin'))
    import ROOT

if sys.platform == 'win32':
    win32api.SetConsoleTitle(console_title) #restore title bar! 
    import pyreadline # fix for tab completion
    pyreadline.parse_and_bind('set show-all-if-ambiguous on')

from ROOT import TMinuit,gMinuit,Long,Double,Minuit2

import numpy as np

def rosenbrock(x):
    """Rosenbrock function, to use for testing."""
    return (1-x[0])**2 + 100*(x[1]-x[0]**2)**2

def rosengrad(x):
    """Gradient of Rosenbrock function, for testing."""
    drdx = -2*((1-x[0])+200*x[0]*(x[1]-x[0]**2))
    drdy = 200*(x[1]-x[0]**2)
    return np.asarray([drdx,drdy])

class FCN(object):
    """Wrap a python callable as an FCN object passable to Minuit.

    __init__() params:
        fcn : A Python callable
        pars : A sequence of starting parameters for fcn
        args : An optional sequence of extra arguments to fcn
        gradient : An optional function to compute the function gradient.
                    Should take a list of parameters and return a list
                    of first derivatives. Broken!
    """
    def __init__(self,fcn,pars,args=(),gradient = None):
        self.fcn = fcn
        self.p = pars
        self.args = args
        self.npars = len(self.p)
        self.iflag = 0
        self.fval = self.fcn(self.p,*self.args)
        self.grad_fcn = gradient

    def __call__(self,nargs,grads,fval,pars,iflag):
        self.p = np.asarray([pars[i] for i in xrange(self.npars)])
        self.iflag = iflag
        self.fval = fval[0] = self.fcn(self.p,*self.args)
        if self.grad_fcn:
            grad = self.grad_fcn(self.p,*self.args)
            self.grads = grads = grad

class Minuit(object):
    """A wrapper class to initialize a minuit object with a numpy array.

    For now, assumes that all parameters are free. 

    Positional args:
        myFCN : A python callable
        params : An array (or other python sequence) of free parameters

    Keyword args:

        limits [None] : a nested sequence of (lower_limit,upper_limit) for each parameter.
        steps [[.1]*npars] : Estimated errors for the parameters, used as an initial step size.
        tolerance [.001] : Tolerance to test for convergence.  Minuit considers convergence to be
                            when the estimated distance to the minimum (edm) is <= .001*up*tolerance,
                            or 5e-7 by default.
        up [.5]  : Change in the objective function that determines 1-sigma error.  .5 if 
                          the objective function is -log(Likelihood), 1 for chi-squared.
        max_calls [10000] : Maximum number of calls to the objective function.
        param_names ['p0','p1',...] : a list of names for the parameters
        args [()] : a tuple of extra arguments to pass to myFCN and gradient.
        gradient [None] : a function that takes a list of parameters and returns a list of 
                          first derivatives of myFCN.  Assumed to take the same args as myFCN.
                          Broken!
        force_gradient [0] : Set to 1 to force Minuit to use the user-provided gradient function.
    """


    def __init__(self,myFCN,params,**kwargs):

        self.limits = np.zeros((len(params),2))
        self.steps = .1*np.ones(len(params))
        self.tolerance = .001
        self.maxcalls = 10000
        self.printMode = 0
        self.up = 0.5
        self.param_names = ['p%i'%i for i in xrange(len(params))]
        self.erflag = Long()
        self.npars = len(params)
        self.args = ()
        self.gradient = None
        self.force_gradient = 0
        self.__dict__.update(kwargs)

        self.params = np.asarray(params)
        self.fcn = FCN(myFCN,self.params,args=self.args,gradient=self.gradient)
        self.fval = self.fcn.fval
        gMinuit = TMinuit(self.npars)
        gMinuit.SetFCN(self.fcn)
        if self.gradient:
            gMinuit.mncomd('SET GRA %i'%(self.force_gradient),self.erflag)
        gMinuit.SetPrintLevel(self.printMode)


        for i in xrange(self.npars):
            gMinuit.DefineParameter(i,self.param_names[i],self.params[i],self.steps[i],self.limits[i][0],self.limits[i][1])

        gMinuit.SetErrorDef(self.up)
        self.minuit = gMinuit

    def minimize(self,method='MIGRAD'):

        self.minuit.mncomd('%s %i %f'%(method, self.maxcalls,self.tolerance),self.erflag)
        for i in xrange(self.npars):
            self.minuit.GetParameter(i,self.params[i],ROOT.Double(0))
        self.fval = self.fcn.fval
        return (self.params,self.fval)

    def errors(self,two_sided = False):
        mat = np.zeros(self.npars**2)
        self.minuit.mnhess()
        self.minuit.mnemat(mat,self.npars)
        return mat.reshape((self.npars,self.npars))
