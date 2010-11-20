"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/stacklike/angularmodels.py,v 1.6 2010/10/21 22:35:35 mar0 Exp $
author: M.Roth <mar0@u.washington.edu>
"""

from uw.utilities.minuit import Minuit
from uw.stacklike.CLHEP import HepRotation
from uw.like import pypsf
import numpy as np
import pylab as py
import uw.utilities.image as im
import scipy.integrate as si
import scipy.special as sp
import scipy.misc as sm
import skymaps as s
import time as t
import copy as cp

rd = 180./np.pi

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
        self.mark='-'
        self.lims=lims
        self.model_par=model_par
        self.free=free
        self.steps=[0.001/rd,0.01]                 #size of step taken by fitter
        self.limits=[[0.0001/rd,360./rd],[1.000000001,1000.]]   #limits of fitter (gamma>1)
        self.name='psf'
        self.header='sigma\tgamma\t'
    
    def __call__(self,diff):
        return self.psf(diff,self.model_par[0],self.model_par[1])

    ## returns the value of the unnormalized psf for a photon
    #  @param photon CLHEP photon
    #  @params pars psf parameters, fed in by fitter
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'Wrong number of PSF parameters, got %d instead of %d'%(len(pars),len(self.model_par))
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

    ## returns the radius in radians of a particular confidence level
    #  @param cl confidence level on interval (0,1)
    def rcl(self,cl):
        if cl<0 or cl>=1:
            return -1
        sig = self.model_par[0]
        g = self.model_par[1]
        return self.recl(cl,sig,g)
     
    ## returns the radius in radians of a particular confidence level
    #  @param cl confidence level on interval (0,1)
    #  @param sig sigma parameter in radians
    #  @param g gamma parameter
    def recl(self,cl,sig,g):
        return sig*np.sqrt(2.*g*((1.-cl)**(1./(-g+1.))-1.))

    ## returns the radius in radians of a particular confidence level
    #  @param cl confidence level on interval (0,1)
    #  @param errs error in the sigma parameter in radians
    #  @param g error in the gamma parameter
    #  @param cov covariance between sigma and gamma
    def clerr(self,cl,errs,errg,cov):
        sig = self.model_par[0]
        g = self.model_par[1]
        dfds = sm.derivative(lambda x: self.recl(cl,x,g),sig,dx=sig/10.)
        dfdg = sm.derivative(lambda x: self.recl(cl,sig,x),g,dx=g/10.)
        acc = dfds*dfds*errs*errs+dfdg*dfdg*errg*errg+2*dfds*dfdg*cov
        return np.sqrt(acc)

    ## setup parameters from containment radii
    #  @param r0 first containment in radians
    #  @param r1 second containment in radians
    #  @param c0 first containment fraction [0-1]
    #  @param c1 second containment fraction [0-1]
    def fromcontain(self,r0,r1,c0,c1):
        pars = [r0/c0*0.68,2.25]
        minuit = Minuit(lambda x: self.cfunc(x[0],x[1],r0,r1,c0,c1),pars,limits=self.limits,tolerance=1e-10,strategy=2,printMode=-1)
        minuit.minimize()
        self.model_par=minuit.params
        #print r0,self.recl(c0,self.model_par[0],self.model_par[1])
        #print r1,self.recl(c1,self.model_par[0],self.model_par[1])

    ## fitting function for containment
    def cfunc(self,sigma,gamma,r0,r1,c0,c1):
        acc = ((r0-self.recl(c0,sigma,gamma))/r0)**2
        acc = acc + ((r1-self.recl(c1,sigma,gamma))/r1)**2
        return acc
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
    def __init__(self,lims,model_par=[0,0,0],free=[False,False,False],rot=[0,0,0],ebar=1e4,ec=0,irf='P6_v3'):
        super(Model,self).__init__(model_par,free)
        self.mark='-'
        self.rot = HepRotation(rot,False)
        cdb = pypsf.CALDBPsf(irf=irf+'_diff')
        if ec!=-1:
            r68=cdb.inverse_integral(ebar,ec,68.)/rd
            r95=cdb.inverse_integral(ebar,ec,95.)/rd
        else:
            r68=(cdb.inverse_integral(ebar,0,68.)+cdb.inverse_integral(ebar,1,68.))/rd/2.
            r95=(cdb.inverse_integral(ebar,0,95.)+cdb.inverse_integral(ebar,1,95.))/rd/2.
        self.PSF = PSF(lims=lims,model_par=[0.001,2.25],free=[False,False])
        self.PSF.fromcontain(r68,r95,0.68,0.95)
        self.model_par=model_par
        self.free=free
        self.lims=lims
        self.ebar=ebar
        self.steps=[20./3600./rd,20./3600./rd,30./3600./rd]
        self.limits=[[-np.pi,np.pi],[-np.pi,np.pi],[-np.pi,np.pi]]
        self.name='psf'
        self.header='Rx\tRy\tRz\t'

    ## returns the value of the unnormalized psf for a photon
    #  @param photon CLHEP photon
    #  @params pars psf parameters, fed in by fitter
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'Wrong number of PSFAlign parameters, got %d instead of %d'%(len(pars),len(self.model_par))
            raise
        rx = pars[0]
        ry = pars[1]
        rz = pars[2]
        rot = HepRotation([rx,ry,rz],False).m(self.rot)
        diff = photon.diff(rot)
        f = self.PSF(diff)
        return f
    
    ## returns the value of the unnormalized psf for a photon
    #  @param diff angular deviation in radians
    #  @param sig sigma parameter in radians
    #  @param g gamma parameter
    #def psf(self,diff,sig,g):
    #    u = diff*diff/(2*sig*sig)
    #    return (1-1/g)*(1+u/g)**-g

    ## return the integral for a photon based on self.lims
    #  @param photon CLHEP photon
    #  @param pars psf parameters, fed in by fitter
    def integrate(self,photon,pars):
        return self.PSF.integrate(photon,self.PSF.model_par)

    ## return the integral of the psf
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        return self.PSF.integral(delmin,delmax)

    ## updates King function parameters, through fitter
    #  @param pars [sigma,gamma], sigma in radians
    def update(self,pars):
        self.model_par=pars

################################################### END PSFALIGN CLASS  ###########################################

################################################### Isotropic CLASS         ###########################################

## Isotropic class
#
# Manages uniform Isotropicround model
class Isotropic(Model):

    ## Constructor
    #
    #  @param lims [min,max], minimum and maximum angular deviations in radians
    def __init__(self,lims):
        super(Model,self).__init__([],[])
        self.mark='-'
        self.model_par=[]
        self.free=[]
        self.lims=lims
        self.steps=[]
        self.limits=[]
        self.name='iso'
        self.header=''
    
    ## returns 1 since Isotropic is constant in (d)**2
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'There should be no parameters for Isotropic, got %d'%len(pars)
            raise
        return 1.
    
    ## returns integral of Isotropic for a photon between self.lims
    def integrate(self,photon,pars):
        um = self.lims[0]*self.lims[0]/2.
        ua = self.lims[1]*self.lims[1]/2.
        return ua-um

    ## returns integral of Isotropic between delmin and delmax
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        um = delmin*delmin/2.
        ua = delmax*delmax/2.
        return ua-um

    ## updates model parameters (none) 
    def update(self,pars):
        self.model_par=pars

################################################### END Isotropic CLASS     ###########################################

################################################### DIFFUSE CLASS           ###########################################

class Diffuse(Model):

    def __init__(self,lims,diff,ebar):
        super(Model,self).__init__([],[])
        self.lims=lims
        self.mark='o'
        self.name='diff'
        self.header=''
        self.model_par=[]
        self.free=[]
        self.steps=[]
        self.limits=[]
        nside = 8192#int(2.*5270./(self.lims[1]*rd))
        print nside
        self.band = s.Band(int(nside))
        self.area = self.band.pixelArea()
        self.diffuse = s.DiffuseFunction(diff,ebar)
        self.ebar = ebar
        self.func = [[],[]]

    def addsrc(self,src):
        start = t.time()
        wsdl = s.WeightedSkyDirList(self.band,src,self.lims[1],True)
        bins = np.arange(0.,100.,1.)
        bins = bins*self.lims[1]/100.
        delt = bins[1]-bins[0]
        seps = []
        mags = []
        for wsd in wsdl:
            diff = wsd.difference(src)
            mag = self.diffuse.value(wsd,self.ebar)
            seps.append(diff)
            mags.append(mag)
        seps = np.array(seps)
        mags = np.array(mags)
        first = self.func[0]==[]
        if first:
            self.func[0]=[0]
        for it in range(len(bins)-1):
            cut = mags[(seps>bins[it])&(seps<bins[it+1])]
            #integral = N(pixels)*Area(Pixels)*(sum(f(pixel)))
            fint = (self.area)*(sum(cut))
            #print len(cut),sum(cut),bins[it]+delt/2,fint
            if first:
                self.func[0].append(bins[it+1])
                self.func[1].append(fint)
            else:
                self.func[1][it] = self.func[1][it] + fint
        self.custom = Custom(self.lims,cp.copy(self.func))
        
        """grid = im.AIT(skyfun=None,size=5.)
        mask = (seps>bins[len(bins)/2])&(seps<bins[len(bins)/2+1])
        ras = np.array([wsd.ra() for wsd in wsdl])
        decs = np.array([wsd.dec() for wsd in wsdl])
        ras = ras[mask]
        decs = decs[mask]
        seldirs = [s.SkyDir(ras[it],decs[it]) for it in range(len(ras))]
        p1 = grid.plot(seldirs,marker='o',c=mags[mask],cmap=py.cm.hot)"""
        del wsdl
        del seps
        del mags
        stop = t.time()
        print 'Took %d seconds to setup the diffuse'%(int(stop-start))

    def value(self,photon,pars):
        return self.custom.value(photon,pars)
    
    def pdf(self,diff):
        return self.custom.pdf(diff)

    def integrate(self,photon,pars):
        return self.custom.integrate(photon,pars)

    def integral(self,delmin,delmax):
        return self.custom.integral(delmin,delmax)

    def update(self,pars):
        self.model_par=pars
################################################### END DIFFUSE CLASS       ###########################################

################################################### CUSTOM CLASS        ###########################################

## Custom class
#  
#  Define a custom angular distribution in counts/degree
class Custom(Model):

    ## Constructor
    #
    #  @param lims [min,max], minimum and maximum angular deviations in radians
    #  @param func [[ang],[f(ang)]] array of angular separations in radians and the pdf of that separation
    def __init__(self,lims,func):
        super(Model,self).__init__([],[])
        self.mark='o'
        self.model_par=[]
        self.free=[]
        self.lims=lims
        self.steps=[]
        self.limits=[]
        self.name='custom'
        self.header=''
        #assume bins centers
        self.seps = []
        self.bins = len(self.seps)
        acc = 0
        self.integ = []
        self.fun =[]
        self.right = []
        for it,fn in enumerate(func[1]):
            acc = acc+fn
            xbar = np.sqrt((func[0][it+1]**2+func[0][it]**2)/2.)
            delt = (func[0][it+1]**2-func[0][0]**2)/2.
            self.fun.append(fn/xbar)
            self.seps.append(xbar)
            self.right.append(func[0][it+1])
            self.integ.append(acc)
        self.fun = np.array(self.fun)/self.integ[len(self.integ)-1]
        self.integ = np.array(self.integ)/self.integ[len(self.integ)-1]
        for it,fn in enumerate(self.integ):
            delt = func[0][it+1]-func[0][it]
            self.integ[it]=fn*delt
        self.seps = np.array(self.seps)
        self.right = np.array(self.right)

    ## linear interpolater
    #  @param val x-value to be interpolated
    #  @param xv x-values to check against
    #  @param yv y-values corresponding to x-values
    def litp(self,val,xv,yv):
        idx = xv.searchsorted(val)
        mi=idx-1
        ma=idx+1
        if idx>=(len(xv)-1):
            mi = len(xv)-3
            idx = mi+1
            ma = mi+2
        if idx==0:
            mi = 0
            idx = 1
            ma= 2
        while yv[mi]==0. and mi>1:
            mi=mi-1
        while yv[ma]==0. and ma<(len(xv)-1):
            ma=ma+1

        ata = np.matrix([[xv[mi]**2+xv[idx]**2+xv[ma]**2,xv[mi]+xv[idx]+xv[ma]],[xv[mi]+xv[idx]+xv[ma],3]])
        atb = np.matrix([[xv[mi]*yv[mi]+xv[idx]*yv[idx]+xv[ma]*yv[ma]],[yv[mi]+yv[idx]+yv[ma]]])
        sol = np.linalg.inv(ata) * atb
        m = sol[0][0].item()
        b = sol[1][0].item()

        rval = m*val+b
        if rval<0:
            return 0
        else:
            return rval

    ## pdf function
    #  @param val x-value to be interpolated
    #  @param xv x-values to check against
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'There should be no parameters for a custom model, got %d'%len(pars)
            raise
        diff = photon.srcdiff()
        return self.litp(diff,self.seps,self.fun)

    ## pdf probability distribution function
    #  @param diff angular separation in radians
    def pdf(self,diff):
        return self.litp(diff,self.seps,self.fun)

    ## cdf function
    #  @param photon CLHEP photon
    #  @param pars function parameters (None)
    def integrate(self,photon,pars):
        ib = self.litp(self.lims[1],self.right,self.integ)
        ia = self.litp(self.lims[0],self.right,self.integ)
        return (ib-ia)
    
    ## returns integral of custom between delmin and delmax
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        ib = self.litp(delmax,self.right,self.integ)
        ia = self.litp(delmin,self.right,self.integ)
        return (ib-ia)

    ## updates model parameters (none) 
    def update(self,pars):
        self.model_par=pars
################################################### END CUSTOM CLASS    ###########################################

################################################### GAUSSIAN            ###########################################

## Gaussian class
#
# Manages model of gaussian in 'd' (narrower than PSF)
class Gaussian(Model):

    ## Constructor
    #  @param lims [min,max], minimum and maximum angular deviations in radians
    #  @param model_par [theta], parameter to describe gaussian width in radians
    #  @param free [True], will fit the parameter theta
    def __init__(self,lims,model_par,free=[True]):
        super(Model,self).__init__(model_par,free)
        self.mark='-'
        self.model_par=model_par
        self.free=free
        self.lims=lims
        self.steps=[0.01/rd]
        self.limits=[[0.001/rd,10./rd]]
        self.name='gaussian'
        self.header='theta\t'
    
    ## returns value of halo distribution for a photon
    #  @param photon CLHEP photon
    #  @param pars parameter [theta], passed on by fitter
    def value(self,photon,pars):
        if len(pars)!=len(self.model_par):
            print 'Wrong number of GAUSSIAN parameters!'
            raise
        theta = pars[0]
        diff = photon.srcdiff()
        f0 = self.psf(diff,theta)
        return f0

    ## returns psf value of halo distributino for a specific value
    #  @param diff angular difference in radians
    #  @param theta width of distribution in radians
    def psf(self,diff,theta):
        return np.exp(-(diff/theta)**2/2)
    
    ## returns integral of halo distribution for a photon with given parameters
    #  @param photon CLHEP photon
    #  @param pars parameter [theta], passed on by fitter
    def integrate(self,photon,pars):
        theta = pars[0]
        um = self.lims[0]
        ua = self.lims[1]
        fint = theta*theta*(np.exp(-(um/theta)**2/2)-np.exp(-(ua/theta)**2/2))
        return fint

    ## returns integral of halo distribution between limits
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        theta = self.model_par[0]
        um = delmin
        ua = delmax
        fint = theta*theta*(np.exp(-(um/theta)**2/2)-np.exp(-(ua/theta)**2/2))
        return fint

    ## updates halo parameters
    #  @params pars [theta] in radians
    def update(self,pars):
        self.model_par=pars

################################################### END GUASSIAN CLASS  ###########################################

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
        self.mark='-'
        self.model_par=model_par
        self.free=free
        self.lims=lims
        self.steps=[0.01/rd]
        self.limits=[[0.001/rd,2./rd]]
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
        return np.exp(-(diff)**4/(theta**4))
    
    ## returns integral of halo distribution for a photon with given parameters
    #  @param photon CLHEP photon
    #  @param pars parameter [theta], passed on by fitter
    def integrate(self,photon,pars):
        theta = pars[0]
        um = self.lims[0]*self.lims[0]/(theta*theta)
        ua = self.lims[1]*self.lims[1]/(theta*theta)
        fint = theta*theta*np.sqrt(np.pi)/4.*(sp.erf(ua)-sp.erf(um))
        return fint

    ## returns integral of halo distribution between limits
    #  @param delmin minimum angle in radians
    #  @param delmax maximum angle in radians
    def integral(self,delmin,delmax):
        theta = self.model_par[0]
        um = delmin*delmin/(theta*theta)
        ua = delmax*delmax/(theta*theta)
        fint = theta*theta*np.sqrt(np.pi)/4.*(sp.erf(ua)-sp.erf(um))
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
        self.mark='-'
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
            self.free.append(not free)
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
    #  @param mode minuit fit strategy: 0-quick 1-normal 2-careful
    def fit(self,photons,free=[],exp=[],mode=1,quiet=True):
        n = sum([x.weight for x in photons])
        pars = []
        frees = []
        lims=[]
        steps=[]
        header = 'Calls\tNtotal\t'

        #set up printing headers and free parameters
        if not quiet:
            print 'Fitting %d models'%len(self.models)
            print '-----------------'
        for x in range(len(self.models)):
            if not quiet:
                print self.models[x].name
            header = header+('N%s\t'%self.models[x].name)
            if len(free)==0:
                pars.append(n/(1.*len(self.models)))
                frees.append(False)
            else:
                pars.append(exp[x])
                frees.append(not free[x])
            lims.append([0,2*n])
            steps.append(np.sqrt(n))
        if not quiet:
            print '-----------------'
        for md in self.models:
            header = header+md.header
                
        for i,par in enumerate(self.par):
            pars.append(par)
            frees.append(self.free[i])
            lims.append(self.limits[i])
            steps.append(self.steps[i])
        header = header+'Likelihood\tsec/call'
        if not quiet:
            print header
        self.clock=t.time()

        self.initial = self.extlikelihood(pars,photons,quiet)
        #minimize likelihood for parameters with respect to the parameters
        #gradient=(lambda x:self.gradient(x,photons,quiet)),
        self.minuit = Minuit(lambda x: self.extlikelihood(x,photons,quiet),pars,fixed=frees,limits=lims,tolerance=1e-3,strategy=mode,printMode=-1)
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
    def extlikelihood(self,pars,photons,quiet=True):
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
            acc = acc - np.log(tacc)*photon.weight
        acc = acc + sum(nest)
        self.calls=self.calls+1

        #print parameter set every 10 function evals
        if self.calls%10==0:
            st = '%d\t'%self.calls
            st = st+'%5.0f\t'%(int(sum(nest)))
            for ne in nest:
                st = st +'%5.0f\t'%(int(ne))
            for spar in mpars:
                st = st+'%1.4f\t'%(spar)
            st = st +'%5.1f\t'%(2*(self.initial-acc))
            ctime = t.time()
            st = st+'%1.1f'%((ctime-self.clock)/self.calls)
            if not quiet:
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