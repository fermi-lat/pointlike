"""Provides a  convenience class to call Minuit, mimicking the interface to scipy.optimizers.fmin."""
import sys
from ROOT import TMinuit,gMinuit,Long,Double,Minuit2
import numpy as np

class FCN(object):
    """Wrap a python callable as an FCN object passable to Minuit.
    
    __init__() params:
        fcn : A Python callable
        pars : A sequence of starting parameters for fcn
        args : An optional sequence of extra arguments to fcn
    """
    def __init__(self,fcn,pars,args=()):
        self.fcn = fcn
        self.p = pars
        self.args = args
        self.npars = len(self.p)
        self.iflag = 0
        self.fval = self.fcn(self.p,*self.args)

    def __call__(self,nargs,grads,fval,pars,iflag):
        self.p = np.asarray([pars[i] for i in xrange(self.npars)])
        self.iflag = iflag
        self.fval = fval[0] = self.fcn(self.p,*self.args)

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
        self.__dict__.update(kwargs)

        #fcn = fcn_factory(myFCN,len(params))
        self.fcn = FCN(myFCN,params)
        self.fval = self.fcn.fval
        gMinuit = TMinuit(self.npars)
        gMinuit.SetFCN(self.fcn)

        for i in xrange(self.npars):
            gMinuit.DefineParameter(i,self.param_names[i],params[i],self.steps[i],self.limits[i][0],self.limits[i][1])

        gMinuit.SetErrorDef(self.up)
        self.minuit = gMinuit

    def minimize(self,method='MIGRAD'):

        self.minuit.mncomd('%s %i %f'%(method, self.maxcalls,self.tolerance),self.erflag)

    def errors(self,two_sided = False):
        mat = np.zeros(self.npars**2)
        self.minuit.mnhess()
        self.minuit.mnemat(mat,self.npars)
        return mat.reshape((self.npars,self.npars))
