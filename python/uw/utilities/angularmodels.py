"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/angularmodels.py,v 1.0 2010/07/29 13:53:17 mar0 Exp $
author: M.Roth <mar0@u.washington.edu>
"""

from uw.utilities.minuit import Minuit
from uw.utilities.CLHEP import HepRotation
import numpy as N
import scipy.integrate as si
import scipy.special as sp
import skymaps as s
import time as t

rd = 180./N.pi

################################################### MODEL CLASS         ###########################################

##   Model class
#
#    Abstraction of angular distribution models
class Model(object):

    def __init__(self,model_par,free):
        self.par = model_par
        self.free = free

################################################### END MODEL CLASS     ###########################################

################################################### PSF CLASS           ###########################################

## PSF class
#
#  Models the angular distribution of a single King function for fitting width and tails

class PSF(Model):
    ## Constructor
    #
    #  @param lims [min,max], minimum and maximum angular deviations in radians
    #  @param model_par [sig,gam], model parameters (sigma, gamma), with sigma in radians
    #  @param free [sig,gam], frees parameters to be fit

    def __init__(self,lims,model_par,free=[True,True]):
        super(PSF,self).__init__(model_par,free)
        self.lims=lims
        self.model_par=model_par
        self.free=free
        self.steps=[0.001/rd,0.01]                 #size of step taken by fitter
        self.limits=[[0.001/rd,10/rd],[1.01,5.]]   #limits of fitter (gamma>1)
        self.name='psf'
        self.header='sigma\tgamma\t'
    
    ## returns the value of the unnormalized psf for a photon
    #  @param photon CLHEP photon
    #  @params pars psf parameters, fed in by fitter
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'Wrong number of PSF parameters!'
            raise
        sig = pars[0]
        g = pars[1]
        diff = photon.srcdiff()
        f = self.psf(diff,sig,g)
        return f
    
    ## returns the value of the unnormalized psf for a photon
    #  @param diff angular deviation in radians
    #  @param sig sigma parameter in radians
    #  @param g gamma parameter
    def psf(self,diff,sig,g):
        u = diff*diff/(2*sig*sig)
        return (1-1/g)*(1+u/g)**-g

    ## return the integral for a photon based on self.lims
    #  @param photon CLHEP photon
    #  @param pars psf parameters, fed in by fitter
    def integrate(self,photon,pars):
        sig = pars[0]
        g = pars[1]
        um = self.lims[0]*self.lims[0]/(2.*sig*sig)
        ua = self.lims[1]*self.lims[1]/(2.*sig*sig)
        f0 = (1+um/g)**(-g+1)
        f1 = (1+ua/g)**(-g+1)
        return sig*sig*(f0-f1)

    ## return the integral of the psf
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        sig = self.model_par[0]
        g = self.model_par[1]
        um = delmin*delmin/(2.*sig*sig)
        ua = delmax*delmax/(2.*sig*sig)
        f0 = (1+um/g)**(-g+1)
        f1 = (1+ua/g)**(-g+1)
        return sig*sig*(f0-f1)

    ## updates King function parameters, through fitter
    #  @param pars [sigma,gamma], sigma in radians
    def update(self,pars):
        self.model_par=pars

################################################### END PSF CLASS       ###########################################

################################################### PSFALIGN CLASS      ###########################################

## PSFAlign class
#
#  Models the angular distribution of a single King function for fitting boresight alignment
class PSFAlign(PSF):

    ## Constructor
    #
    #  @param lims [min,max], minimum and maximum angular deviations in radians
    #  @param model_par [sig,gam], model parameters (sigma, gamma), with sigma in radians
    #  @param free [sig,gam], frees parameters to be fit
    def __init__(self,lims,model_par=[0,0,0],free=[False,False,False],rot=[0,0,0]):
        super(Model,self).__init__(model_par,free)
        self.rot = HepRotation(rot,False)
        self.model_par=model_par
        self.free=free
        self.lims=lims
        self.steps=[20./3600./rd,20./3600./rd,30./3600./rd]
        self.limits=[[-N.pi,N.pi],[-N.pi,N.pi],[-N.pi,N.pi]]
        self.name='psf'
        self.header='Rx\tRy\tRz\t'
    
    ## returns the value of the unnormalized psf for a photon
    #  @param photon CLHEP photon
    #  @params pars psf parameters, fed in by fitter
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'Wrong number of PSFAlign parameters!'
            raise
        sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
        g = s.IParams.gamma(float(photon.energy),int(photon.event_class))
        rx = pars[0]
        ry = pars[1]
        rz = pars[2]
        rot = HepRotation([rx,ry,rz],False).m(self.rot)
        diff = photon.diff(rot)
        f = self.psf(diff,sig,g)
        return f
    
    ## returns the value of the unnormalized psf for a photon
    #  @param diff angular deviation in radians
    #  @param sig sigma parameter in radians
    #  @param g gamma parameter
    def psf(self,diff,sig,g):
        u = diff*diff/(2*sig*sig)
        return (1-1/g)*(1+u/g)**-g

    ## return the integral for a photon based on self.lims
    #  @param photon CLHEP photon
    #  @param pars psf parameters, fed in by fitter
    def integrate(self,photon,pars):
        sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
        g = s.IParams.gamma(float(photon.energy),int(photon.event_class))
        um = self.lims[0]*self.lims[0]/(2.*sig*sig)
        ua = self.lims[1]*self.lims[1]/(2.*sig*sig)
        f0 = (1+um/g)**(-g+1)
        f1 = (1+ua/g)**(-g+1)
        return sig*sig*(f0-f1)

    ## return the integral of the psf
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
        g = s.IParams.gamma(float(photon.energy),int(photon.event_class))
        um = delmin*delmin/(2.*sig*sig)
        ua = delmax*delmax/(2.*sig*sig)
        f0 = (1+um/g)**(-g+1)
        f1 = (1+ua/g)**(-g+1)
        return sig*sig*(f0-f1)

    ## updates King function parameters, through fitter
    #  @param pars [sigma,gamma], sigma in radians
    def update(self,pars):
        self.model_par=pars

################################################### END PSF CLASS       ###########################################

################################################### BACKG CLASS         ###########################################

## Backg class
#
# Manages uniform background model
class Backg(Model):

    ## Constructor
    #
    #  @param lims [min,max], minimum and maximum angular deviations in radians
    def __init__(self,lims):
        super(Model,self).__init__([],[])
        self.model_par=[]
        self.free=[]
        self.lims=lims
        self.steps=[]
        self.limits=[]
        self.name='back'
        self.header=''
    
    ## returns 1 since background is constant in (d)**2
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'There should be no parameters for a uniform background!'
            raise
        return 1.
    
    ## returns integral of background for a photon between self.lims
    def integrate(self,photon,pars):
        um = self.lims[0]*self.lims[0]/2.
        ua = self.lims[1]*self.lims[1]/2.
        return ua-um

    ## returns integral of background between delmin and delmax
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        um = delmin*delmin/2.
        ua = delmax*delmax/2.
        return ua-um

    ## updates model parameters (none) 
    def update(self,pars):
        self.model_par=pars
################################################### END BACKG CLASS     ###########################################

################################################### HALO CLASS          ###########################################

## Halo class
#
# Manages model of gaussian halo (in d**2)
class Halo(Model):

    ## Constructor
    #  @param lims [min,max], minimum and maximum angular deviations in radians
    #  @param model_par [theta], parameter to describe gaussian width in radians
    #  @param free [True], will fit the parameter theta
    def __init__(self,lims,model_par,free=[True]):
        super(Model,self).__init__(model_par,free)
        self.model_par=model_par
        self.free=free
        self.lims=lims
        self.steps=[0.01/rd]
        self.limits=[[0.001/rd,10./rd]]
        self.name='halo'
        self.header='theta\t'
    
    ## returns value of halo distribution for a photon
    #  @param photon CLHEP photon
    #  @param pars parameter [theta], passed on by fitter
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'Wrong number of HALO parameters!'
            raise
        theta = pars[0]
        diff = photon.srcdiff()
        f0 = self.psf(diff,theta)
        return f0

    ## returns psf value of halo distributino for a specific value
    #  @param diff angular difference in radians
    #  @param theta width of distribution in radians
    def psf(self,diff,theta):
        return N.exp(-(diff)**4/(theta**4))
    
    ## returns integral of halo distribution for a photon with given parameters
    #  @param photon CLHEP photon
    #  @param pars parameter [theta], passed on by fitter
    def integrate(self,photon,pars):
        theta = pars[0]
        um = self.lims[0]*self.lims[0]/(theta*theta)
        ua = self.lims[1]*self.lims[1]/(theta*theta)
        fint = theta*theta*N.sqrt(N.pi)/4.*(sp.erf(ua)-sp.erf(um))
        return fint

    ## returns integral of halo distribution between limits
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        theta = self.model_par[0]
        um = delmin*delmin/(theta*theta)
        ua = delmax*delmax/(theta*theta)
        fint = theta*theta*N.sqrt(N.pi)/4.*(sp.erf(ua)-sp.erf(um))
        return fint

    ## updates halo parameters
    #  @params pars [theta] in radians
    def update(self,pars):
        self.model_par=pars

################################################### END HALO CLASS      ###########################################

################################################### COMPOSITE MODEL CLASS #########################################

## CompositeModel class
#
# Manages a number of different models and performs optimizes parameters
class CompositeModel(object):
    
    ## Constructor
    #  just set up all of the members
    def __init__(self):
        self.models = []  #list of Model objects
        self.par = []     #all parameters for models
        self.free = []    #corresponding free parameters for models
        self.steps=[]     #corresponding step size for fitter
        self.limits=[]    #corresponding limits for parameter values
        self.nest = []    #Estimator for number of photons from each model
        self.calls=0      #number of calls of the extended likelihood function

    ## adds a Model to be fit
    #  @param model Model object
    def addModel(self,model):
        self.models.append(model)
        for free in model.free:
            self.free.append(free)
        for par in model.model_par:
            self.par.append(par)
        for steps in model.steps:
            self.steps.append(steps)
        for limits in model.limits:
            self.limits.append(limits)

    ## maximizes the extended likelihood of the model parameters for the list of photons
    #  @param photons list of CLHEP photons
    #  @param free frees all number estimators, assumed to be true
    #  @param exp if free is set to [..False..], exp = [..Ni..], where Ni is the estimator for the number of photons in model i
    def fit(self,photons,free=[],exp=[]):
        n = len(photons)
        pars = []
        frees = []
        lims=[]
        steps=[]
        header = 'Calls\tNtotal\t'

        #set up printing headers and free parameters
        for x in range(len(self.models)):
            header = header+('N%s\t'%self.models[x].name)
            if len(free)==0:
                pars.append(n/(1.*len(self.models)))
                frees.append(True)
            else:
                pars.append(exp[x])
                frees.append(free[x])
            lims.append([0,2*n])
            steps.append(N.sqrt(n))
        
        for md in self.models:
            header = header+md.header
                
        for i,par in enumerate(self.par):
            pars.append(par)
            frees.append(self.free[i])
            lims.append(self.limits[i])
            steps.append(self.steps[i])
        header = header+'Likelihood\tsec/call'
        print header
        self.clock=t.time()
        
        #minimize likelihood for parameters with respect to the parameters
        self.minuit = Minuit(lambda x: self.extlikelihood(x,photons),pars,free=frees,limits=lims,tolerance=1e-3,strategy=1,printMode=-1)
        self.minuit.minimize()

        self.nest = self.minuit.params[:len(self.models)] # first parameters are the number estimators
        self.par=self.minuit.params[len(self.models):]    # last parameters are from the models
        lastp = 0

        #update model parameters
        for i,model in enumerate(self.models):
            prs = len(model.model_par)
            self.models[i].update(self.par[lastp:lastp+prs])
            lastp=lastp+prs

    ## returns the extended likelihood of models and number estimators
    #  @param pars [N1...Nn,p11...,pnm] 'n' number estimators (N) with 'm' parameters (p)
    #  @param photons list of CLHEP photons to fit to
    def extlikelihood(self,pars,photons):
        nest = pars[:len(self.models)]
        mpars = pars[len(self.models):]
        acc = 0

        # extended likelihood calculation
        #
        # logL = sum(Ni,1,n)-log(sum(sum(Ni*model(j|pi),1,n),1,#of photons))
        # pi are the parameters of the ith model
        # j is the jth photon
        # Ni is the number of photons associated with the ith model
        # model(j|pi) is the value of the ith model given the photon j and parameters pi
        for photon in photons:
            tacc = 0
            lastp = 0
            for i,model in enumerate(self.models):
                prs = len(model.model_par)
                f0 = model.value(photon,mpars[lastp:lastp+prs])
                fint = model.integrate(photon,mpars[lastp:lastp+prs])
                tacc = tacc + nest[i]*f0/fint
                lastp=lastp+prs
            acc = acc - N.log(tacc)
        acc = acc + sum(nest)
        self.calls=self.calls+1

        #print paramter set ever 10 function evals
        if self.calls%10==0:
            st = '%d\t'%self.calls
            st = st+'%5.0f\t'%(int(sum(nest)))
            for ne in nest:
                st = st +'%5.0f\t'%(int(ne))
            for spar in mpars:
                st = st+'%1.4f\t'%(spar)
            st = st +'%5.1f\t'%acc
            ctime = t.time()
            st = st+'%1.1f'%((ctime-self.clock)/self.calls)
            print st
        return acc

    ## integrates the total model from delmin to delmax, normalized to omin,omax (for differentials)
    ## this is useful for plotting purposes
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    #  @param omin minimum overall angle in radians
    #  @param omax maximum overall angle in radians
    def integral(self,delmin,delmax,omin,omax):
        acc = 0
        for it,model in enumerate(self.models):
            acc = acc + self.nest[it]*model.integral(delmin,delmax)/model.integral(omin,omax)
        return acc

################################################### END COMPOSITE MODEL CLASS #########################################