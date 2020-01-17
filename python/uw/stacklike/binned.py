import skymaps as s
import pylab as py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from uw.stacklike.stacklike import *
from uw.stacklike.angularmodels import *
import uw.stacklike.stcaldb as uss
from uw.utilities.minuit import Minuit
#import uw.utilities.assigntasks as ua
from uw.like.pycaldb import CALDBManager
from uw.like.pypsf import CALDBPsf
from uw.like.quadform import QuadForm,Ellipse
from uw.stacklike.limits import BayesianLimit
from scipy.misc import derivative
from uw.like2.pipeline import cluster,engines
import scipy.optimize as so
import scipy.misc as sm
import scipy.stats as sst
import scipy.special as spec
import scipy.integrate as si
import os as os
import copy as cp
import glob as glob
import time as t
from ROOT import Double
import math
import sys
import string
import uw.pulsar.py_exposure as pe
from uw.like.Models import PowerLaw,ExpCutoff,LogParabola
from uw.stacklike import dataset
import cPickle

def format_error(v,err):
    if v>0 and err>0 and not np.isinf(err):
        logz = math.floor(math.log10(err))-1
        #print (logz)
        z = 10**logz
        err = round(err/z)*z
        v = round(v/z)*z
        if v>0:
            return '     %10s    (1 +/- %1.3f)'%(v,err/v)
        else:
            return '     %10s    (1 +/- %1.3f)'%(0,0)
    else:
        return '     %10s    (1 +/- %1.3f)'%(0,0)


####################### PSR file naming conventions ######################
#
#    pulsar source list text file:
#
#    (psr).txt
# name     ra     dec
# PSR0835-4510    128.8357   -45.1863
#
#   ft1 files:
#   (psr)on-ft1.fits
#   (psr)off-ft1.fits
#

##############################################  Namespace Defaults #############################################
cachedir = dataset.basedir+'cache/'
figdir = dataset.basedir+'figures/'
pulsdir = dataset.basedir+'data/pulsar/'
pulsars = [['velaclean',0.17/0.37,[0.0,0.11,0.63,0.69],[0.19,0.56]],['gemclean',0.2/0.3,[0.08,0.18,0.6,0.70],[0.25,0.55]]]#,('crab',0.45/0.55,[0.0,0.09,0.25,0.48,0.87,1.0],[0.09,0.25,0.48,1.0])]                                       #pulsar sourcelist name 
agnlist = ['psfstudy']
fdays=10*365
rd = 180./np.pi
INFINITY = np.Infinity#1e80#
ZERO = 1e-40 
ergs = 1.602e-6
np.seterr(all='ignore')

###############################################  CombinedLike Class ############################################

## Implements the likelihood analysis outlined by Matthew Wood
## available here: http://fermi-ts.phys.washington.edu/pointlike/common/psf_lnl.pdf
class CombinedLike(object):


    ## constructor for Combinedlikelihood class
    #  @param pulsars (sourcelist,alpha) names of pulsar list and ratio of on/off pulse
    #  @param agnlist list of agn sources
    #  @param irf intial response function
    def __init__(self,**kwargs):
        self.pulsars = pulsars
        self.agnlist = agnlist
        self.cachedir = cachedir
        self.irf = 'P7SOURCE_V6'
        self.cache = True
        self.verbose = False
        self.veryverbose = False
        self.pulse_ons = []         #pulsar on pulse angular distributions
        self.pulse_offs = []        #pulsar off pulse angular distributions
        self.agns = []              #agn angular distributions
        self.angbins=[]             #angular bins
        self.midpts = []            #midpoints of angular bins
        self.ponhists=[]            #pulsar on pulse histograms
        self.pofhists=[]            #pulsar off pulse histograms
        self.agnhists=[]            #agn histograms
        self.nbins = 0              #number of angular bins
        self.params=[]
        self.halomodel=''
        self.haloparams=[]
        self.prof=[]
        self.ctmin=0.3
        self.ctmax=1.0
        self.mode=-1
        self.__dict__.update(kwargs)
        self.TS = 0
        self.tol = 1e-3
        self.qandd=False                 #quick and dirty errors, use this if you don't care about parameter errors
        self.usegrad=False
        self.strategy=2
        self.data2mcb=False
        self.lcalls=0

    def __str__(self):
        verbose = '\n'
        verbose = verbose + '--------- PSF    -------------\n'
        verbose = verbose + 'l-edge(deg)  r-edge(deg)    fraction   fraction error\n'
        for it, x in enumerate(self.psf):
            tstr = '%8.3g   %8.3g'%(self.angbins[it],self.angbins[it+1]), '     %1.3f   (1 - %1.3f + %1.3f)'%(self.psf[it],(self.psf[it]-self.psfl[it][0])/self.psf[it],(self.psfl[it][1]-self.psf[it])/self.psf[it])#format_error(self.psf[it],self.psfe[it])
            verbose = verbose + string.join(tstr) + '\n'
        verbose = verbose +'\n'
        verbose = verbose + '--------- Pulsars-------------\n'
        for it,nj in enumerate(self.Npj):
            verbose =verbose + 'N(%s) = %1.1f - %1.1f + %1.1f : Model chisq(%1.1f)\n'%(self.pulsars[it][0],nj,nj-self.Npjl[it][0],self.Npjl[it][1]-nj,sum((self.ponhists[it]-self.ponmodels[it])**2/self.ponmodels[it]))
        verbose = verbose + '\n'
        verbose = verbose + '--------- AGN    -------------\n'
        for it,nj in enumerate(self.Naj):
            verbose = verbose + 'Npsf(%s) = %1.1f - %1.1f + %1.1f\n'%(self.agnlist[it],nj,nj-self.Najl[it][0],self.Najl[it][1]-nj)
            verbose = verbose + 'Niso(%s) = %1.1f - %1.1f + %1.1f\n'%(self.agnlist[it],self.Ni[it],self.Ni[it]-self.Nil[it][0],self.Nil[it][1]-self.Ni[it])
            if self.halomodel!='' and it==(len(self.Naj)-1):
                verbose = verbose + 'Nhalo = %1.1f - %1.1f + %1.1f\n'%(self.Nh[-1],self.Nh[-1]-self.Nhl[-1][0],self.Nhl[-1][1]-self.Nh[-1])
                verbose = verbose + 'Size(deg) = %1.3f : Model chisq(%1.1f)\n'%(self.haloparams[0]*rd,sum((self.agnhists[it]-self.agnmodels[it])**2/self.agnmodels[it]))
                verbose = verbose + 'TS = -2(LogL(Nhalo) - LogL(0)) = %1.2f\n'%(self.TS)
                tspval = self.tspval()
                verbose = verbose + 'Significance = %1.1f    p-val = %1.4f\n'%(tspval[1],tspval[0])
            verbose = verbose + '\n'
        return verbose
    
    ######################################################################
    #          Loads data from FT1 files for stacked sources             #
    ######################################################################
    ## loads photons from ft1 and bin the data
    #  @param minroi minimum angular separation (degrees)
    #  @param maxroi maximum angular separation (degrees)
    #  @param emin minimum energy (MeV)
    #  @param emax maximum energy (MeV)
    #  @param tmin minimum time range (MET)
    #  @param tmax maximum time range (MET)
    #  @param ctype conversion type (0:front, 1:back, -1:all)
    def loadphotons(self,minroi,maxroi,emin,emax,tmin,tmax,ctype):
        self.minroi = minroi
        self.maxroi = maxroi
        self.ctype = ctype
        self.emin=emin
        self.emax=emax
        self.ebar=np.sqrt(self.emin*self.emax)

        tag = '%1.2f%1.2f%1.2f%1.2f%1.0f%1.0f%1.0f%1.2f%1.2f%s'%(minroi,maxroi,emin,emax,tmin,tmax,ctype,self.ctmin,self.ctmax,self.irf)

        #load pulsar data
        thetabar = 0.
        photons = 0
        for psr in self.pulsars:

            if os.path.exists(self.cachedir+'%son%s.npy'%(psr[0],tag)):
                if self.verbose:
                    print ('Loaded %s on pulse from cache from: %son%s'%(psr[0],psr[0],tag))
                hist = np.load(self.cachedir+'%son%s.npy'%(psr[0],tag))
                self.pulse_ons.append(hist)
            else:
                sl = StackLoader(name=psr[0]+'on',ctmin=self.ctmin,ctmax=self.ctmax,quiet=not self.veryverbose)
                sl.loadphotons(minroi,maxroi,emin,emax,tmin,tmax,ctype)
                sl.getds()
                photons = photons+len(sl.photons)
                thetabar = thetabar+sum([p.ct for p in sl.photons])
                self.pulse_ons.append(np.array(cp.copy(sl.ds)))
                if self.cache:
                    sl.spickle(self.cachedir+'%son%s.pickle'%(psr[0],tag))
                    np.save(self.cachedir+'%son%s.npy'%(psr[0],tag),np.array(cp.copy(sl.ds)))
                del sl

            if os.path.exists(self.cachedir+'%soff%s.npy'%(psr[0],tag)):
                if self.verbose:
                    print ('Loaded %s off pulse from cache from: %soff%s'%(psr[0],psr[0],tag))
                hist = np.load(self.cachedir+'%soff%s.npy'%(psr[0],tag))
                self.pulse_offs.append(hist)
            else:
                sl = StackLoader(name=psr[0]+'off',ctmin=self.ctmin,ctmax=self.ctmax,quiet=not self.veryverbose)
                sl.loadphotons(minroi,maxroi,emin,emax,tmin,tmax,ctype)
                photons = photons+len(sl.photons)
                thetabar = thetabar+sum([p.ct for p in sl.photons])
                sl.getds()
                self.pulse_offs.append(np.array(cp.copy(sl.ds)))
                if self.cache:
                    sl.spickle(self.cachedir+'%son%s.pickle'%(psr[0],tag))
                    np.save(self.cachedir+'%soff%s.npy'%(psr[0],tag),np.array(cp.copy(sl.ds)))
                del sl

        #load agn data
        self.backs = []
        self.sagns =[]
        for it,lists in enumerate(self.agnlist):
            if os.path.exists(self.cachedir+'%s%s.npy'%(lists,tag)):
                if self.verbose:
                    print ('Loaded %s from cache: %s%s'%(lists,lists,tag))
                hist = np.load(self.cachedir+'%s%s.npy'%(lists,tag))
                sl = cPickle.load(open(self.cachedir+'%s%s.pickle'%(lists,tag)))#StackLoader(name=lists,ctmin=self.ctmin,ctmax=self.ctmax,quiet=not self.veryverbose)
                #sl.bindata()
                sl.solveback()
                self.sagns.append(sl.Npsf)
                self.backs.append(sl.Nback)
                print (len(sl.ds),sl.Npsf+sl.Nback)
                sl.makeplot('agnplot%d'%(it),bin=40)
                del sl
                self.agns.append(hist)
            else:
                sl = StackLoader(name=lists,ctmin=self.ctmin,ctmax=self.ctmax,quiet=not self.veryverbose)
                sl.loadphotons(minroi,maxroi,emin,emax,tmin,tmax,ctype)
                photons = photons+len(sl.photons)
                thetabar = thetabar+sum([p.ct for p in sl.photons])
                sl.getds()
                self.agns.append(np.array(cp.copy(sl.ds)))
                #sl.bindata()
                sl.solveback()
                self.sagns.append(sl.Npsf)
                self.backs.append(sl.Nback)
                if self.cache:
                    sl.spickle(self.cachedir+'%s%s.pickle'%(lists,tag))
                    np.save(self.cachedir+'%s%s.npy'%(lists,tag),np.array(cp.copy(sl.ds)))
                del sl
        self.thetabar = 0.5*(self.ctmin+self.ctmax) if photons==0 else thetabar/photons
        self.photons=photons
            

    ######################################################################
    #   Bins the FT1 data into angular bins that are adaptive or fixed   #
    ######################################################################
    ## adaptively bins angular distributions in angle or sqrt
    #  @param bins number of angular bins, -1 for sqrt(N) bins
    def bindata(self,bins=8,halomod='',par=0):
        alldata = []
        chist=[]
        #adaptive binning

        if halomod=='':
            if bins>0:
                
                if sum([len(puls) for puls in self.pulse_ons])>2*len(self.pulse_ons)*bins:
                    print ('Using Pulsars for angular binning: %d photons'%(sum([len(puls) for puls in self.pulse_ons])))

                    #determine needed bins by combining all data (except background)

                    for puls in self.pulse_ons:
                        for sep in puls:
                            alldata.append(sep)
                            chist.append(1.)
                    for it1,puls in enumerate(self.pulse_offs):
                        for sep in puls:
                            alldata.append(sep)
                            chist.append(-self.pulsars[it1][1])
                    #for sep in self.agns[0]:
                    #    alldata.append(sep)
                    alldata = np.array(alldata)
                    key = np.argsort(alldata)
                    chist = np.array(chist)[key]
                    alldata = alldata[key]
                    chist = np.array([sum(chist[:x+1]) for x in range(len(chist))])

                # sqrt(N) binning
                else:
                    print ('Using AGN for angular binning')
                    maxph = 0
                    bestph = 0
                    for it,tagn in enumerate(self.agns):
                        if len(tagn)>maxph:
                            maxph,bestph=len(tagn),it
                    for sep in self.agns[bestph]:
                        alldata.append(sep)
                        chist.append(1.)
                    """if len(alldata)>20:
                        bins = int(np.sqrt(len(alldata))/2.)
                    else:
                        bins = min(32,int(np.sqrt(len(self.agns[0]))/2.))
                    xbins = (np.arange(0,bins,1)/(1.*bins))**2
                    if len(alldata)>0:
                        minimum = min(min(alldata),min(self.agns[0]))*rd
                    else:
                        minimum = min(self.agns[0])*rd
                    xbins =  minimum + xbins*(self.maxroi-minimum)
                    self.angbins = xbins/rd"""

                    alldata = np.array(alldata)
                    key = np.argsort(alldata)
                    chist = np.array(chist)[key]
                    alldata = alldata[key]
                    chist = np.array([sum(chist[:x+1])-self.backs[bestph]*(alldata[x]*rd/self.maxroi)**2 for x in range(len(chist))])
                chist = chist/max(chist)
                cumm = np.array([(1.*x+1.)/len(alldata) for x in range(len(alldata))])      #cumulative dist function
                ct = (1.*np.arange(0,bins+1,1))/bins                                          #cumulative fractions, [0,1/bins,2/bins...(bins-1)/bins]
                mask = np.array([max(0,len(chist[chist<x])-1) for x in ct])
                xbins = alldata[mask]#np.array([alldata[max(0,len(chist[chist<x])-1)] for x in ct])         #bin edges corresponding to fractions
                self.angbins = xbins
                cbins = []
                for it in range(len(self.angbins)):
                    flag = True
                    for par in cbins:
                        if par == self.angbins[it]:
                            flag=False
                    if flag:
                        cbins.append(self.angbins[it])
                self.angbins = np.array(cbins)
            else:
                self.angbins = np.histogram(self.agns[0],bins=np.sqrt(len(self.agns[0])))[1]
        else:
            halomodel = eval(halomod)
            mod = halomodel(lims=[self.minroi/rd,self.maxroi/rd],model_par=par)
            self.angbins = np.array([mod.rcontain(x) for x in np.linspace(1e-2/rd,1.,bins+1)])
        #did you already do it?
        if self.ponhists==[]:
            for it1,puls in enumerate(self.pulse_ons):
                if len(puls)>0:
                    self.ponhists.append(np.histogram(puls,self.angbins)[0])
                    self.pofhists.append(np.histogram(self.pulse_offs[it1],self.angbins)[0])
            for it1,agns in enumerate(self.agns):
                if len(agns)>0:
                    self.agnhists.append(np.histogram(agns,self.angbins)[0])

        self.angbins = self.angbins*rd           #convert to degrees
        self.nbins = len(self.angbins)-1         #lop of bin edge for bin center
        self.iso = np.array([(self.angbins[it+1]**2-self.angbins[it]**2)/(max(self.angbins)**2-min(self.angbins)**2) for it in range(self.nbins)])  #isotropic model
        self.midpts = np.array([(self.angbins[it+1]+self.angbins[it])/2. for it in range(self.nbins)])           #bin midpoint
        self.widths = np.array([(self.angbins[it+1]-self.angbins[it])/2. for it in range(self.nbins)])           #bin widths
        self.areas = np.array([self.angbins[it+1]**2-self.angbins[it]**2 for it in range(self.nbins)])           #jacobian
        self.hmd = np.zeros(self.nbins)

    ######################################################################
    #    Generate the profile likelihood in one parameter                #
    ######################################################################
    ## calculates the profile likelihood for one parameter
    #  @param ip parameter index [0:len(params)-1]
    #  @param x value of parameter to be profiled
    def profile(self,ip,x,cverb=False,**kwargs):
        self.__dict__.update(kwargs)
        tol = abs(0.0001/self.likelihood(np.log10(self.params)))

        if np.isscalar(x):
            if x>self.limits[ip][1] or x<self.limits[ip][0]:
                return INFINITY
            params = cp.copy(self.params)
            fixed = cp.copy(self.fixed)
            free = cp.copy(self.free)
            self.params[ip]=x
            self.fixed[ip]=True
            self.free[ip]=False
            direc = self.loggrad2(np.log10(self.params))
            direc = direc/np.sqrt(sum(direc*direc))
            direc=direc[self.free]
            best = so.fmin_powell(lambda z: self.likelihood(self.setargs(z)),np.log10(self.params[self.free]),ftol=tol,disp=0,full_output=1)
            #self.params = 10**self.setargs(best[0])
            #minuit = Minuit(self.likelihood,np.log10(params),#gradient=self.gradient,force_gradient=1,
            #                 fixed=fixed,limits=np.log10(self.limits),strategy=2,tolerance=tol,printMode=(self.mode if not self.verbose else 0))
            #minuit.minimize()
            #fval = minuit.fval
            fval = best[1]
            if cverb:
                print (x,self.lmax-fval)
            self.tpars=10**(self.setargs(best[0]))
            self.params=params
            self.fixed=fixed
            self.free=free
            #del minuit
            return fval
        else:
            params = cp.copy(self.params)
            fixed = cp.copy(self.fixed)
            vals = []
            for xval in x:
                params = cp.copy(self.params)
                fixed = cp.copy(self.fixed)
                free = cp.copy(self.free)
                self.params[ip]=cp.copy(xval)
                self.fixed[ip]=True
                self.free[ip]=False
                direc = self.loggrad2(np.log10(self.params))
                direc = direc/np.sqrt(sum(direc*direc))
                direc=direc[self.free]
                best = so.fmin_powell(lambda z: self.likelihood(self.setargs(z)),cp.copy(np.log10(self.params[self.free])),ftol=tol,disp=0,full_output=1)
                #self.params = 10**self.setargs(best[0])
                #minuit = Minuit(self.likelihood,np.log10(params),#gradient=self.gradient,force_gradient=1,
                #                 fixed=fixed,limits=np.log10(self.limits),strategy=2,tolerance=tol,printMode=(self.mode if not self.verbose else 0))
                #minuit.minimize()
                #fval = minuit.fval
                fval = best[1]
                vals.append(fval)
                if cverb:
                    print (xval,self.lmax-fval)
                self.tpars=10**(self.setargs(best[0]))
                self.params=params
                self.fixed=fixed
                self.free=free
            return np.array(vals)

    #######################################################################
    #  Find upper and lower changes to profile likelihood for a parameter #
    #######################################################################
    ## Finds the values of the parameter corresponding to profile logl changes of 'delt', 0.5 correspond to 1-sigma
    #  @param ip parameter index [0:len(params)-1]
    #  @param delt change in likelihood (positive value), 0.5 corresponds to 1 sigma, 2 to 2 sigma, etc
    #  @param startpar optional guess for parameter value corresponding to 'delt'
    def findlims(self,ip,delt,startpar=-1):
        #print ('Checking p%d of %d'%(ip+1,len(self.params)))
        par = self.params[ip]
        self.lmax = self.profile(ip,par)

        #check to see if limits exist already, otherwise step the default size, unless the parameter is particularly small
        if self.prof==[]:
            step = 1 if par>1e-4 else (2./par-1)*10
            #print ('Using default step')
        else:
            if self.prof[ip][1]>0:
                step = self.prof[ip][1]/par-1.
                #print ('Using upper limit for step')
            else:
                step = 1 if par>1e-4 else (2./par-1)*10
                #print ('Using default step')
        if par*(1+step)>self.limits[ip][1]:
            step=0.9*self.limits[ip][1]/par-1.

        #UPPER LIMIT
        #check if a parameter step has been specified
        if startpar<0:
            xr = [np.log(par),np.log(par)+np.log(1+step/2.),np.log(par*(1+step))]
            yr = [0,self.lmax-self.profile(ip,np.exp(xr[1])),self.lmax-self.profile(ip,np.exp(xr[2]))]
        else:
            xr = [np.log(startpar),np.log(startpar+1),np.log(startpar+2)]
            yr = [self.lmax-self.profile(ip,np.exp(xr[0])),self.lmax-self.profile(ip,np.exp(xr[1])),self.lmax-self.profile(ip,np.exp(xr[2]))]

        xi,yi = cp.copy(xr),cp.copy(yr)          #copy initial values

        #try the log-profiler first
        bestup,bestval,failed = self.effprofiler(xr,yr,delt,ip,True,True)
        if failed:
            #print ('Log fitter failed, trying linear')
            xmin = self.limits[ip][0]*1.1
            xmax = self.limits[ip][1]*0.5
            xr = [xmin,0.5*(xmin+xmax),xmax]#np.exp(np.array(xi))
            yr = [self.lmax-self.profile(ip,xr[0]),self.lmax-self.profile(ip,xr[1]),self.lmax-self.profile(ip,xr[2])]#yi
            bestup,bestval,failed = self.effprofiler(xr,yr,delt,ip,True,False)

        #LOWER LIMIT
        #use upper limit as a guess for lower
        step = -bestup/par+1 if par>1e-4 else -0.5
        if par*(1+step)<self.limits[ip][0]:
            step=1.1*self.limits[ip][0]/par-1.
        if startpar<0:
            minpt = max(ZERO,par*(1+step))
            xr = [np.log(minpt),np.log(par*(1+step/2.)),np.log(par)]
            yr = [self.lmax-self.profile(ip,np.exp(xr[0])),self.lmax-self.profile(ip,np.exp(xr[1])),0]
        else:
            xr = [np.log(startpar),np.log(startpar+1),np.log(startpar+2)]
            yr = [self.lmax-self.profile(ip,np.exp(xr[0])),self.lmax-self.profile(ip,np.exp(xr[1])),self.lmax-self.profile(ip,np.exp(xr[2]))]

        xi,yi = cp.copy(xr),cp.copy(yr)          #copy initial values

        #try the log-profiler first
        bestlow,bestval,failed = self.effprofiler(xr,yr,delt,ip,False,True)
        if failed:
            #print ('Log fitter failed, trying linear')
            xmin = self.limits[ip][0]*1.1
            xmax = self.limits[ip][1]*0.5
            xr = [xmin,0.5*(xmin+xmax),xmax]#np.exp(np.array(xi))
            yr = [self.lmax-self.profile(ip,xr[0]),self.lmax-self.profile(ip,xr[1]),self.lmax-self.profile(ip,xr[2])]#yi
            bestlow,bestval,failed = self.effprofiler(xr,yr,delt,ip,False,False)
        if failed:
            bestlow = 2*par-bestup
        return[bestlow,bestup]
    
    #######################################################################
    #  Limit finding helper class                                         #
    #######################################################################
    ## efficient likelihood profiler - based on the goldensearch method of fitting parabolas to the likelihood profile
    #  xr initial x positions
    #  yr initial likelihood profile values corresponding to xr
    #  delt value of likelihood change sought
    #  ip parameter index
    #  upper upper limit search, otherwise lower limit
    #  log xr is log(param)
    def effprofiler(self,xr,yr,delt,ip,upper=True,log=True):
        xi,yi=xr,yr
        cdelt = delt*10
        best = INFINITY
        bestpar = np.log(self.params[ip]) if log else self.params[ip]
        steps=0
        chisq=1.
        xvals = [xr[0],xr[1],xr[2]]
        yvals = [yr[0],yr[1],yr[2]]
        sols=[]

        #don't take too many steps, should converge in less than 20 if likelihood profile is sensible
        #chisq is just the relative variance of the x-positions, should exit if bracketed properly
        #exit when the change in likelihood is within 2% of the correct change in likelihood
        while steps<20 and abs(cdelt+delt)>delt/50. and chisq>0.001:

            #calculate the variance in the x-position
            xbar = np.exp(sum(xr)/len(xr)) if log else sum(xr)/len(xr)
            chisq = np.sqrt(sum((np.exp(np.array(xr))-xbar)**2))/xbar if log else np.sqrt(sum((np.array(xr)-xbar)**2))/xbar

            #solve the parabola
            try:
                A,B,C = self.quadsolve(xr,yr)
                sols.append([A,B,C])
            except:
                return bestpar,best,True
            
            #make sure profile will pass through the requested y-value
            try:
                disc = B**2-4*A*(C+delt)
                if disc>0:
                    
                    #make sure we're picking the right position
                    if (A<=0 and upper) or (A>0 and not upper):
                        xr0 = -B/(2*A)+abs(np.sqrt(disc)/(2*A))
                    elif (A>=0 and upper) or (A<0 and not upper):
                        xr0 = -B/(2*A)-abs(np.sqrt(disc)/(2*A))
                    else:
                        return bestpar,best,True
                else:
                    return bestpar,best,True
            except:
                return bestpar,best,True

            #calculate likelihood at the best new position
            yr0 = self.lmax-self.profile(ip,np.exp(xr0)) if log else self.lmax-self.profile(ip,xr0)
            estpar = np.exp(-B/(2*A)) if log else -B/(2*A)     #check the predicted max
            xvals.append(xr0)
            yvals.append(yr0)
            cdelt = yr0
            format = '%10.3f '*13
            print (format%(xr[0],xr[1],xr[2],xr0,yr[0],yr[1],yr[2],cdelt,best,bestpar,chisq,estpar,self.params[ip]))
            steps+=1
            if abs(cdelt+delt)<best:
                best=abs(cdelt+delt)
                bestpar=xr0
            
            #exchange the furthest point in analysis with new x-value
            if xr0<xr[0]:
                xr[2],yr[2]=xr[1],yr[1]
                xr[1],yr[1]=xr[0],yr[0]
                xr[0]=xr0
                yr[0]=yr0
                continue
            if xr0>xr[0] and xr0<xr[1]:
                xr[2],yr[2]=xr[1],yr[1]
                xr[1]=xr0
                yr[1]=yr0
                continue
            if xr0>xr[1] and xr0<xr[2]:
                xr[0],yr[0]=xr[1],yr[1]
                xr[1]=xr0
                yr[1]=yr0
                continue
            if xr0>xr[2]:
                xr[0],yr[0]=xr[1],yr[1]
                xr[1],yr[1]=xr[2],yr[2]
                xr[2]=xr0
                yr[2]=yr0

        #some plotting
        if True:
            py.figure(4,figsize=(8,8))
            py.clf()
            parfunc = lambda x: x[1][0]*(x[0]**2) + x[1][1]*x[0] + x[1][2]
            xspc = np.linspace(min(xvals),max(xvals),100)
            for it,xv in enumerate(xvals):
                py.plot(xv,yvals[it],'o')
            for sol in sols:
                py.plot(xspc,parfunc([xspc,sol]))
            py.legend(range(len(yvals)+len(sols)),loc=3 if upper else 4)
            flag = 'upper' if upper else 'lower'
            py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/tmp/%d%s.png'%(ip,flag))
        failed = (abs(cdelt+delt)>delt/50.) and chisq>0.1
        bestpar = np.exp(bestpar) if log else bestpar
        
        #return the best fit parameter, the difference from the delt, and if the fit failed
        return bestpar,best,failed
    
    def testst(self):
        npars = len(self.params)
        if not self.fixed[npars-1]:
            self.TS = 2*(self.profile(npars-1,ZERO)-self.profile(npars-1,self.params[-1]))
        else:
            self.TS = 0
        return self.TS

    ## parabolic solver
    #  @param xr x-values
    #  @param yr y-values
    def quadsolve(self,xr,yr):
        Am = np.matrix([[xr1**2,xr1,1.] for xr1 in xr])
        bm = np.matrix([[yr1] for yr1 in yr])
        sol = np.linalg.inv(Am)*bm
        A,B,C = sol[0][0].item(),sol[1][0].item(),sol[2][0].item()
        return [A,B,C]

    ######################################################################
    #    Generate the profile likelihood in one parameter                #
    ######################################################################
    def getalllims(self):
        lims = []
        for it,par in enumerate(self.params):
            lims.append(self.findlims(it,0.5))
        return np.array(lims)

    ######################################################################
    #    Generate the profile likelihood in one parameter                #
    ######################################################################
    def find_logl_change(self,it,x,ref,upper=True,dlogl=0.5):
        x10 = 10**x
        ref10 = 10**ref
        if x10>self.limits[it][1] or x10<self.limits[it][0]:
            return INFINITY
        #print (x,ref,x10,ref10)
        fact = x10>ref10 if upper else x10<ref10
        prf = self.profile(it,x10)
        dl = abs(self.lmax-prf+dlogl) if fact else INFINITY
        #print (x10,ref10,self.lmax,prf,dl)
        t.sleep(0.05)
        return dl

    def freeze(self,ip,val):
        self.fixed[ip]=True
        self.params[ip]=val

    def thaw(self,ip,val):
        self.fixed[ip]=False
        self.params[ip]=val


    ######################################################################
    #    Sets up results of fitting                                      #
    ######################################################################
    def setmembers(self):
        psrs = len(self.ponhists)
        agns = len(self.agnhists)
        #scale = sum(self.psf)                                                                                          #normalize psf
        nmu = self.nbins - 1        
        self.psf = cp.copy(self.params[:nmu])
        self.psfe = np.array([self.errs[i] for i in range(nmu)])                             #PSF errors
        # Compute the value and error of the nth PSF weight
        self.psf = np.append(self.psf,[1-np.sum(self.psf)])
        self.psfl = [prof for prof in self.prof[:nmu]]
        dxdmu = np.zeros(len(self.params))
        dxdmu[:nmu] = -1
        mun_err = np.sqrt(sum([psfe**2 for psfe in self.psfe]))#np.dot(np.dot(self.cov,dxdmu),dxdmu))        
        self.psfe = np.append(self.psfe,[mun_err])
        self.psfl.append([0,0])
        self.psfl = np.array(self.psfl)
        self.Npj  = cp.copy(self.params[nmu:nmu+psrs])               #PSR number est
        self.Npje = self.errs[nmu:nmu+psrs]                 #PSR number est errors
        self.Npjl = self.prof[nmu:nmu+psrs]
        self.Naj  = cp.copy(self.params[nmu+psrs:nmu+psrs+agns])     #AGN number est
        self.Naje = self.errs[nmu+psrs:nmu+psrs+agns]           #AGN number errors
        self.Najl = self.prof[nmu+psrs:nmu+psrs+agns]
        self.Ni   = cp.copy(self.params[nmu+psrs+agns:nmu+psrs+agns+agns]) #isotropic number est
        self.Nie  = self.errs[nmu+psrs+agns:nmu+psrs+agns+agns]         #isotropic number est errors
        self.Nil  = self.prof[nmu+psrs+agns:nmu+psrs+agns+agns]

        if self.halomodel!='':
            self.Nh  = cp.copy(self.params[nmu+psrs+agns+agns:])  #halo est
            self.Nhe = self.errs[nmu+psrs+agns+agns:]           #halo err                               #halo err
            self.Nhl = self.prof[nmu+psrs+agns+agns:]
            self.testst()
        
        ########  calculate background estimators  ###########
        self.vij = []
        self.vije = []

        #loop over pulsars
        for it1,row in enumerate(self.ponhists):

            tvij = []
            tvije = []
            N = self.Npj[it1]               #pulsar number estimator
            a = self.pulsars[it1][1]        #ratio of on-off

            #loop over angular bins
            for it2,n in enumerate(row):
                #n on pulse in current bin
                m = self.psf[it2]                                        #PSF in current bin
                b = self.pofhists[it1][it2]                              #off pulse in current bin
                v = self.backest(a,n,b,m,N)                              #background number estimator
                tvij.append(v)

                #estimate errors by propagation
                dvdm = sm.derivative(lambda x: self.backest(a,n,b,x,N),m,m/10.)   #background/psf derivative
                dvdN = sm.derivative(lambda x: self.backest(a,n,b,m,x),N,N/10.)   #background/psr number estimator derivative
                #cov = self.cov[it2][nmu+it1]                              #covariance of psf and number estimator
                Ne = self.Npje[it1]
                me = self.psfe[it2]
                ve = np.sqrt((dvdm*me)**2+(dvdN*Ne)**2)      #naive error propagation
                tvije.append(ve)

            self.vij.append(tvij)
            self.vije.append(tvije)

        self.vij = np.array(self.vij)                                        #PSR background number estimator
        self.vije = np.array(self.vije)                                      #PSR background number estimator errors

        self.ponmodels = np.array([self.Npj[it1]*self.psf+self.pulsars[it1][1]*self.vij[it1] for it1 in range(len(self.ponhists))])
        self.pofmodels = np.array([self.vij[it1] for it1 in range(len(self.pofhists))])
        self.agnmodels = np.array([self.Naj[it1]*self.psf + self.Ni[it1]*self.iso + self.Nh[it1]*self.hmd for it1 in range(len(self.agnhists))])


    def makehalo(self,pars):
        halomodel = eval(self.halomodel)
        mod = halomodel(lims=[min(self.angbins)/rd,max(self.angbins)/rd],model_par=pars)
        mint = mod.integral(min(self.angbins)/rd,max(self.angbins)/rd)
        self.hmd = np.array([mod.integral(self.angbins[it]/rd,self.angbins[it+1]/rd)/mint for it in range(self.nbins)])
        self.hmd = self.hmd/sum(self.hmd)

    ######################################################################
    #    Sets up parameters                                              #
    ######################################################################
    def setuppars(self,**kwargs):
        self.__dict__.update(kwargs)
        psrs = ''
        for psl in self.pulsars:
            psrs = psrs + '%s\t'%psl[0]

        agns = ''
        for agnl in self.agnlist:
            agns = agns + '%s\t'%agnl

        pars = ''
        for par in self.haloparams:
            pars = pars + '%1.6f\t'%par

        if self.verbose:
            print ('')
            print ('**********************************************************')
            print ('*                                                        *')
            print ('*             Combined Likelihood Analysis               *')
            print ('*                                                        *')
            print ('**********************************************************')
            print ('Pulsars: %s'%psrs)
            print ('AGNs: %s'%agns)
            print ('Using halo model: %s'%self.halomodel)
            print ('Using parameters: %s'%pars)
            print ('**********************************************************')

        #sf = CALDBPsf(CALDBManager(irf=self.irf))
        #
        #fint = psf.integral(self.ebar,self.ctype,max(self.angbins)/rd,min(self.angbins)/rd)
        il = uss.IrfLoader(self.irf)
        rpars = il.params(self.ebar,self.ctype)#average_psf(self.emin,self.emax,self.ctmin,self.ctmax,self.ctype)
        rpars[5] = rpars[5]*rpars[2]
        rpars[0]*=rd
        rpars[3]*=rd
        self.psfm = self.makepsf(rpars)#np.array([psf.integral(self.ebar,self.ctype,self.angbins[it+1]/rd,self.angbins[it]/rd)/fint for it in range(self.nbins)])
        self.psfm = np.append(self.psfm,[1-np.sum(self.psfm)])
        self.psfm = np.array([max(ZERO,ps) for ps in self.psfm])
        psrs = len(self.ponhists)
        agns = len(self.agnhists)
        self.Nh=[ZERO for agn in self.agnlist]
        self.Nhe=[ZERO for agn in self.agnlist]
        if self.params==[]:
            #set up initial psf parameters (uniform)
            self.params = [self.psfm[x] for x in range(self.nbins-1)]
            self.limits = [[ZERO,1] for x in range(self.nbins-1)]
            self.fixed = [False for x in range(self.nbins-1)]

            alims = [max(sum(agnhist),1) for agnhist in self.agnhists]

            #pulsar estimators
            for hist in self.ponhists:
                self.params.append(max(ZERO,sum(hist)))
                self.limits.append([ZERO,sum(hist)*10])
                self.fixed.append(False)

            #agn estimators
            for it,hist in enumerate(self.agnhists):
                #est = sum(hist)-hist[-1]/self.areas[-1]*(self.maxroi**2-self.minroi**2)
                self.params.append(self.sagns[it])#est)
                self.limits.append([ZERO,sum(hist)*10])
                self.fixed.append(False)

            #iso estimator
            for it,hist in enumerate(self.agnhists):
                #est = hist[-1]/self.areas[-1]*(self.maxroi**2-self.minroi**2)
                self.params.append(self.backs[it])#est)
                self.limits.append([ZERO,sum(hist)*10])
                self.fixed.append(len(self.agnhists)==0)

        else:
            cnum = len(self.params)
            self.params=[par for par in self.params[:cnum-len(self.agnlist)]]
            self.fixed=[par for par in self.fixed[:cnum-len(self.agnlist)]]
            self.limits=[par for par in self.limits[:cnum-len(self.agnlist)]]
            #print (len(self.params),len(self.fixed),len(self.limits))

        if self.halomodel=='':
            self.hmd=np.zeros(self.nbins)
        else:
            halomodel = eval(self.halomodel)
            if self.haloparams[0]<0:
                return INFINITY
            self.makehalo(self.haloparams)

        for it,agn in enumerate(self.agnlist):
            self.params.append(self.Nh[it] if (it+1)==agns else ZERO)
            self.limits.append([ZERO,sum(self.agnhists[it])*10])
            self.fixed.append(not(self.halomodel!='' and it==(len(self.agnlist)-1)))
            #print (len(self.params),len(self.fixed),len(self.limits))

        if self.verbose:
            print ('Setting up Minuit and maximizing')
            header = ''+string.join(['Np  \t' for nj in self.pulsars])+string.join(['Na  \t' for nj in self.agnlist])+string.join(['Ni  \t' for ni in self.agnlist])+string.join(['Nh  \t' for nh in self.agnlist])+'like'
            print (header)

        """if self.halomodel=='':
            for it in range(self.nbins-1):
                self.fixed[it]=False
            for it in range(psrs):
                self.fixed[self.nbins-1+it]=False
            for it in range(agns):
                self.fixed[self.nbins-1+psrs+it]=(it==agns-1)
                self.fixed[self.nbins-1+psrs+it+agns]=(it==agns-1)
                self.fixed[self.nbins-1+psrs+it+agns]=True
        else:
            for it in range(self.nbins-1):
                self.fixed[it]=True
            for it in range(psrs):
                self.fixed[self.nbins-1+it]=True
            for it in range(agns):
                self.fixed[self.nbins-1+psrs+it]= not (it==agns-1)
                self.fixed[self.nbins-1+psrs+it+agns]= not (it==agns-1)
                self.fixed[self.nbins-1+psrs+it+agns]= not (it==agns-1)"""

        self.params=np.array(self.params)
        self.fixed = np.array(self.fixed)
        self.free = np.array([not fix for fix in self.fixed])
        self.limits = np.array(self.limits)


    ######################################################################
    #    Sets up Minuit and determines maximum likelihood parameters     #
    ######################################################################
    ## maximizes likelihood
    # possible keyword arguments
    # @param halomodel string corresponding to any angular model in angularmodels.py
    # @param haloparams array of model parameters for halomodel
    def fit(self,custom=False,**kwargs):
        self.__dict__.update(kwargs)

        if not custom:
            self.setuppars()

        self.free = np.array([not fix for fix in self.fixed])
        ival = self.likelihood(np.log10(self.params))
        tol = abs(0.0001/ival)
        if self.verbose:
            print ('Using ftol: %1.3e with l0=%1.1f'%(tol,ival))
        ############  setup Minuit and optimize  ###############
        #print (self.gradient(np.log10(self.params),False))
        #print (self.params)
        self.best = so.fmin_powell(lambda z: self.likelihood(self.setargs(z)),np.log10(self.params[self.free]),ftol=tol,disp=0,full_output=1)
        #best = so.fmin_ncg(lambda z: self.likelihood(self.setargs(z)),np.log10(self.params[self.free]),fprime=lambda z: self.loggrad2(self.setargs(z)),disp=0,full_output=1)
        self.params = 10**self.setargs(self.best[0])
        #print (self.params)
        #if self.usegrad:
        #    self.minuit = Minuit(self.likelihood,np.log10(self.params),gradient=self.gradient,force_gradient=1,
        #                     fixed=self.fixed,limits=np.log10(self.limits),strategy=2,tolerance=tol,printMode=self.mode)
        #else:
        self.minuit = Minuit(self.likelihood,np.log10(cp.copy(self.params)),#gradient=self.gradient,force_gradient=1,
                             fixed=self.fixed,limits=np.log10(self.limits),strategy=2,tolerance=tol,printMode=(self.mode if not self.verbose else 0))
        #self.minuit = POWELL(lambda x: self.likelihood(x),np.log10(self.params),np.log10(self.limits),self.fixed)#so.fmin_powell(lambda x: self.likelihood(x),np.log10(self.params),disp=0,full_output=1)
        self.minuit.minimize()
        self.lmax = self.best[1]#self.minuit.fval#self.minuit[1]#
        if self.verbose:
            print ('Likelihood value: %1.1f'%self.minuit.fval#[1]#)
            print ('**********************************************************')
        #print (self.gradient(self.minuit.params,False))
        ###########  get parameters and errors #################
        #self.params = 10**(self.minuit.params)#[0])#
        self.cov = self.minuit.errors()
        #print (self.cov)
        if self.qandd:
            self.errs = 10**(np.sqrt(np.diag(self.cov)))
            self.prof = np.array([[self.params[it]/self.errs[it],self.params[it]*self.errs[it]] for it in range(len(self.errs))])
            self.errs-=1
            self.errs*=self.params
        else:
            self.prof = np.array([self.findlims(it,0.5) if not self.fixed[it] else [self.params[it],self.params[it]] for it in range(len(self.params))])
            self.errs = np.array([np.sqrt(abs(self.prof[it][1]-self.params[it])*abs(self.prof[it][0]-self.params[it])) for it in range(len(self.params))])

        self.setmembers()

        return self.lmax

    
    def setargs(self,xarr):
        it = 0
        retarr = []
        #print (xarr)
        for it2,par in enumerate(self.params):
            if self.fixed[it2]:
                retarr.append(np.log10(par))
            else:
                retarr.append(xarr[it])
                it+=1
        return np.array(retarr)
    
    #######################################################################
    #                  Add Halo into data                                 #
    #######################################################################
    def addhalo(self,halomodel,haloparams,nhalo):
        halomodel = eval(halomodel)
        if haloparams[0]<0:
            return INFINITY
        mod = halomodel(lims=[min(self.angbins)/rd,max(self.angbins)/rd],model_par=haloparams)
        mint = mod.integral(min(self.angbins)/rd,max(self.angbins)/rd)
        hmd = np.array([mod.integral(self.angbins[it]/rd,self.angbins[it+1]/rd)/mint for it in range(self.nbins)])
        hmd = hmd/sum(hmd)*nhalo
        newdata = sst.poisson.rvs(hmd)
        print (newdata)
        self.agnhists[-1]+= newdata
        

    ######################################################
    #       Generate Poisson distributed MC ang bins     #
    ######################################################
    def data2mc(self):
        if not self.data2mcb:
            self.backups = [cp.copy(self.ponhists),cp.copy(self.pofhists),cp.copy(self.agnhists)]
            self.data2mcb=True
        self.ponhists = np.array([[0 if bin==0 else sst.poisson.rvs(bin) for bin in ponhis] for ponhis in self.ponmodels])
        self.pofhists = np.array([[0 if bin==0 else sst.poisson.rvs(bin) for bin in pofhis] for pofhis in self.pofmodels])
        self.agnhists = np.array([[0 if bin==0 else sst.poisson.rvs(bin) for bin in agnhis] for agnhis in self.agnmodels])
        """for it in range(len(self.ponhists)):
            print (self.backups[0][it])
            print (self.ponhists[it])
            print (self.backups[1][it])
            print (self.pofhists[it])
        for it in range(len(self.agnhists)):
            print (self.backups[2][it])
            print (self.agnhists[it])"""

    #######################################################
    #        Calculate sample P-values                    #
    #######################################################
    def realize(self,num):
        self.params[-1]=ZERO
        self.fixed[-1]=True
        TS0 = self.TS
        l0 = self.fit(custom=True,qandd=True)
        print (self)
        self.fixed[-1]=False
        self.data2mcb=False
        print (TS0,self.Nh[-1])
        models = [cp.copy(self.ponmodels),cp.copy(self.pofmodels),cp.copy(self.agnmodels)]
        itr = 0
        TSs = []
        for it in range(num):
            self.data2mc()
            self.fit(custom=True,qandd=True)
            TS1 = self.TS#self.profile(len(self.params)-1,ZERO)-self.profile(len(self.params)-1,self.Nh[-1])
            TSs.append(TS1)
            tTSs = np.array(TSs)
            if itr>0:
                pm = self.wilson(itr+1,it+1)
                pv = (itr+1.)/(it+1.)
                print ('%1.3f %1.1f '%(TS1,self.Nh[-1]),itr,it,'%1.3f %1.3f %1.3f'%(pm[0],pv,pm[1]),' (>1):%d (>2):%d (>3):%d'%(len(tTSs[tTSs>1]),len(tTSs[tTSs>2]),len(tTSs[tTSs>3])))
            else:
                print ('%1.3f %1.1f '%(TS1,self.Nh[-1]),itr,it,' (>1):%d (>2):%d (>3):%d'%(len(tTSs[tTSs>1]),len(tTSs[tTSs>2]),len(tTSs[tTSs>3])))
            if TS1>(TS0):
                itr+=1
            self.ponmodels,self.pofmodels,self.agnmodels=models
        self.mc2data()
        self.params[-1]=self.Naj[-1]
        self.fixed[-1]=False
        self.fit(custom=True,qandd=True)
        return itr

    def wilson(self,sig,trial):
        p = 1.*sig/trial
        denom = (1.+1./trial)
        cent = (p+1./(2*trial))/denom
        pm = np.sqrt(p*(1-p)/trial+1./(4.*(trial)**2))/denom
        return [cent-pm,cent+pm]

    ########################################################
    #      Go from MC to data                              #
    ########################################################
    def mc2data(self):
        self.ponhists = self.backups[0]
        self.poffhists = self.backups[1]
        self.agnhists =self.backups[2]
        self.data2mcb=False

    #######################################################################
    #                  Summary of Likelihood Analysis                     #
    #######################################################################
    ## printResults() - summary
    def printResults(self):
        print (str(self))

    def gradient2(self,params):
        grad=[]
        for it in range(len(params)):
            derv = derivative(lambda x: self.likeone(params,it,x),params[it],dx=1e-6)
            grad.append(derv)
        return np.array(grad)

    def loggrad2(self,lparams):
        return (10**lparams)*np.log(10.)*self.gradient2(10**lparams)

    def likeone(self,params,which,val):
        pars = cp.copy(params)
        pars[which]=val
        return self.likelihood(np.log10(pars))

    ######################################################################
    #      Likelihood function from (3) in paper                         #
    ######################################################################
    ## likelihood function
    #  @param params likelihood function parameters
    def likelihood(self,params,verb=False):#npij,bpij,naij,mi,Npj,Naj,Ni,Nh=0):

        params = 10**params

        if verb:
            print (params)
            t.sleep(0.05)

        psrs = len(self.ponhists)
        agns = len(self.agnhists)

        npij = self.ponhists
        bpij = self.pofhists
        naij = self.agnhists

        nmu = self.nbins - 1

        mi = params[:nmu]
        Npj = params[nmu:nmu+psrs]
        Naj = params[nmu+psrs:nmu+psrs+agns]
        Ni = params[nmu+psrs+agns:nmu+psrs+agns+agns]
        Nh = params[nmu+psrs+agns+agns:]
        acc = 0
        #set to true for slow,verbose output
        if (1-np.sum(mi))<0:
            if self.veryverbose:
                print ('Bad PSF: %1.3f'%(np.sum(mi)))
                print (mi)
            return 0
        mi = np.append(mi,[1-np.sum(mi)])
        if self.veryverbose:
            print ('**************************************************************')
            print ('--------------------------------------------------------------')
            print ('                        Pulsars                               ')
            print ('--------------------------------------------------------------')
            print ('mu\tn\tmod\tb\tvi\tNi\tcont1\tcont2\tacc')

        ########################################
        #          first sum in (3)            #
        ########################################
        #loop over pulsars
        for it1,row in enumerate(npij):
            if self.veryverbose:
                print ('--------------------------------------------------------------')
                print ('                        %s                               '%self.pulsars[it1][0])
                print ('--------------------------------------------------------------')
            N = Npj[it1]                    #pulsar number estimator
            a = self.pulsars[it1][1]        #ratio of on-off

            #loop over angular bins
            for it2,n in enumerate(row):
                #n on pulse in current bin
                m = mi[it2]                                     #PSF in current bin
                b = bpij[it1][it2]                              #off pulse in current bin

                v = self.backest(a,n,b,m,N)                     #get background estimator

                #catch negative log terms
                lterm = N*m+a*v
                #if lterm >0. or v<0.:
                #if self.verbose:
                #    print (lterm,v,'Pulsar')
                #   lterm = 1 if lterm<0 else lterm
                #v = 1 if v<0 else v
                #return INFINITY
                cont1 =  N*m + a*v
                if lterm>0:
                    cont1 += -n*np.log(lterm)
                acc = acc + cont1
                
                cont2 = v
                if v>0:
                    cont2 += -b*np.log(v)
                acc = acc + cont2

                if self.veryverbose:
                    print ('%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f'%(m,n,lterm,b,v,N*m,cont1,cont2,acc))
                    t.sleep(0.05)
        
        if self.veryverbose:
            print ('--------------------------------------------------------------')
            print ('                        AGN                                   ')
            print ('--------------------------------------------------------------')
            print ('mu\tn\tlterm\tNi\tIi\tHi\tcont1\tacc')
        ########################################
        #         second sum in (3)            #
        ########################################
        #loop over agn
        for it1,row in enumerate(naij):

            if self.veryverbose:
                print ('--------------------------------------------------------------')
                print ('                        %s                               '%self.agnlist[it1])
                print ('--------------------------------------------------------------')
            #loop over angular bins
            for it2,bin in enumerate(row):

                #make sure log term is proper
                lterm = Naj[it1]*mi[it2]+Ni[it1]*self.iso[it2] + Nh[it1]*self.hmd[it2]
                #lterm = 1 if lterm<0 else 
                if lterm<0.:
                    if self.verbose:
                        print (it2,Naj[it1],mi[it2],Ni[it1],self.iso[it2],Nh[it1],self.hmd[it2],lterm,'AGN')
                    #    lterm=1
                    #return INFINITY

                cont1 = Naj[it1]*mi[it2] + Ni[it1]*self.iso[it2] + Nh[it1]*self.hmd[it2]
                if lterm>0:
                    cont1 += -bin*np.log(lterm)
                acc = acc + cont1

                if self.veryverbose:
                    print ('%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f'%(mi[it2],bin,lterm,Naj[it1]*mi[it2],Ni[it1]*self.iso[it2],Nh[it1]*self.hmd[it2],cont1,acc))
                    t.sleep(0.05)
        if self.verbose:
            if (self.lcalls%20)==0:
                header = ''+string.join(['Np  \t' for nj in self.pulsars])+string.join(['Na  \t' for nj in self.agnlist])+string.join(['Ni  \t' for ni in self.agnlist])+string.join(['Nh  \t' for nh in self.agnlist])+'like'
                print (header)
            vals = ''+string.join(['%1.1f\t'%nj for nj in Npj])+string.join(['%1.1f\t'%nj for nj in Naj])+string.join(['%1.1f\t'%nj for nj in Ni])+string.join(['%1.1f\t'%nj for nj in Nh])+'%1.4f'%(acc)
            print (vals)
            t.sleep(0.05)
        self.lcalls+=1
        return acc

    ######################################################################
    #      Find single or double PSF fits to fractions                   #
    ######################################################################
    ## finds sigma gamma and fraction
    # @param double fits a double PSF
    def fitpsf(self,double=False,pverbose=False):
        psf = uss.IrfLoader(self.irf)#CALDBPsf(CALDBManager(irf=self.irf))
        de = 0.45
        stest = psf.rcontain(psf.params(self.ebar,self.ctype)[:-1],0.39)*rd#psf.inverse_integral(self.ebar,self.ctype,39.)
        if double:
            if True:
                #nc,nt,gc,gt,sc,st,w = psf.get_p(self.ebar,self.ctype)
                sc,gc,nc,st,gt,nt,ea=psf.params(self.ebar,self.ctype)
                pars = [sc*rd,gc,0.5,st*rd,gt]
                lims = [[0,100],[1,100],[0,100],[1,100],[0.5-de,0.5+de]]
            else:
                #gc,si,w = psf.get_p(self.ebar,self.ctype)
                sc,gc,nc,st,gt,nt,ea=psf.params(self.ebar,self.ctype)
                pars = [sc*rd,gc*1.2,0.5,sc*rd,gc[0]*0.8]
                lims = [[0,100],[1,100],[0,100],[1,100],[0.5-de,0.5+de]]
        else:
            if True:
                #nc,nt,gc,gt,sc,st,w = psf.get_p(self.ebar,self.ctype)
                sc,gc,nc,st,gt,nt,ea=psf.params(self.ebar,self.ctype)
                pars = [stest,gc]
                lims = [[0,10],[1,10]]
            else:
                #gc,si,w = psf.get_p(self.ebar,self.ctype)
                sc,gc,nc,st,gt,nt,ea=psf.params(self.ebar,self.ctype)
                pars = [stest,gc]
                lims = [[0,10],[1,10]]
        print (pars)
        tolerance = abs(0.0025/self.psflikelihood(pars))
        print (tolerance)
        tmin = so.fmin_powell(lambda x: self.psflikelihood(x,pverbose) if x[0]>0 else 1e40,pars,ftol=tolerance,disp=0,full_output=1)
        #tmin = so.fmin_powell(lambda x: self.psfchisq(x) if x[0]>0 else 1e40,pars,disp=0,full_output=1)
        print (tmin[0])

        if not double:
            psf1 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=tmin[0])
            print (psf1.rcl(0.68),psf1.rcl(0.95))
            print (psf.rcontain(psf.params(self.ebar,self.ctype)[:-1],0.68)*rd,psf.rcontain(psf.params(self.ebar,self.ctype)[:-1],0.95)*rd)
        else:
            tp68=psf.rcontain(psf.params(self.ebar,self.ctype)[:-1],0.68)*rd
            tp95=psf.rcontain(psf.params(self.ebar,self.ctype)[:-1],0.95)*rd
            psf1 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=tmin[0][0:2])
            psf2 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=tmin[0][3:5])
            alph = tmin[0][2]
            tpars = [par for par in tmin[0]]
            tpars.append(1.-tpars[2])
            ftot = alph*psf1.integral(0,40) + (1-alph)*psf2.integral(0,40)
            #print (psf1.rcl(0.68),psf1.rcl(0.95))
            #print (getfrac(0.68),getfrac(0.95))   #alph*psf1.rcl(0.68)+(1-alph)*psf2.rcl(0.68),alph*psf1.rcl(0.95)+(1-alph)*psf2.rcl(0.95)
            print (tp68,tp95)
            tp68=psf.rcontain(tpars,0.68)
            tp95=psf.rcontain(tpars,0.95)
            print (tp68,tp95)
            #print (psf.inverse_integral(self.ebar,self.ctype,68.),psf.inverse_integral(self.ebar,self.ctype,95.))
        return np.insert(tmin[0],0,self.ebar)

    ######################################################################
    #      Find single or double PSF fits to fractions                   #
    ######################################################################
    ## finds sigma gamma and fraction
    # @param double fits a double PSF
    def psfchisq(self,pars):
        if len(pars)>2:
            pars2 = [par for par in pars]
            pars2.append(1-pars[2])
        tints = self.makepsf(pars2)
        #print (tints-self.psf[:-1])
        chisqa = [((tints[it]-self.psf[it])/tints[it])**2 if self.psf[it]>0 else np.Infinity for it in range(len(tints))]
        chisq = sum(chisqa)
        #print (string.join(['%1.6f '%par for par in pars]),chisq)
        #t.sleep(0.25)
        return chisq

    def makepsf(self,pars):
        if len(pars)==2:
            psf1 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=[pars[0],pars[1]])
            fint = np.array([psf1.integral(self.angbins[it],self.angbins[it+1]) for it in range(self.nbins-1)])/psf1.integral(psf1.lims[0],psf1.lims[1])
        else:
            psf1 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=[pars[0],pars[1]])
            psf2 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=[pars[3],pars[4]])
            ftot = pars[2]*psf1.integral(psf1.lims[0],psf1.lims[1])+pars[5]*psf2.integral(psf1.lims[0],psf1.lims[1])
            fint = np.array([pars[2]*psf1.integral(self.angbins[it],self.angbins[it+1])+pars[5]*psf2.integral(self.angbins[it],self.angbins[it+1]) for it in range(self.nbins-1)])/ftot
        return fint

    ######################################################################
    #   Calculates the likelihood of a particular single King function   #
    ######################################################################
    ## 
    # @param pars either 2 params [sigma,gamma] for single king, or 6 for double [score,gcore,ncore,stail,gtail,ntail]
    def psflikelihood(self,pars,pverbose=False):

        if len(pars)==2:
            sig,gam=pars
            if sig<0 or gam<=1.:
                return 0
        else:
            sig,gam,nc,sig2,gam2=pars
            if sig<0 or sig2<0 or nc<0 or nc>1 or gam<=1 or gam2<=1 or gam>5 or gam2>5:
                if self.verbose:
                    print ('Bad PSF')
                return 0
            pars = [par for par in pars]
            pars.append(1.-nc)
        fint = self.makepsf(pars)

        params = cp.copy(self.params)
        fixed = cp.copy(self.fixed)
        free = cp.copy(self.free)
        for it,fin in enumerate(fint):
            self.fixed[it]=True
            self.params[it]=fin
            self.free[it]=False
            #self.limits[it]=[fin,fin]

        ival = self.likelihood(np.log10(self.params))
        tol = abs(0.01/ival)
        tpars,fval = so.fmin_powell(lambda z: self.likelihood(self.setargs(z)),np.log10(self.params[self.free]),ftol=tol,disp=0,full_output=1)[0:2]
        if pverbose:
            tstr = ''.join(['%1.3f\t'%(par) for par in pars])
            tstr += ''.join(['%1.1f\t'%(10**par) for par in tpars])
            print tstr,'%1.1f'%(self.lmax-fval)
        
        self.params = cp.copy(params)
        self.fixed = cp.copy(fixed)
        self.free = cp.copy(free)
        del params
        del fixed
        del free
        del tpars
        #memory issues?
        #del self.minuit
        del fint
        #del psf1
        #if len(pars)!=2:
        #    del psf2
        #print (pars,fval)
        #t.sleep(0.05)
        return fval

    ######################################################################
    #      Makes plots of PSR, AGN fits, PSF and background residuals    #
    ######################################################################
    ## outputs a plot of the distributions
    # @param name output PNG filename
    def makeplot(self,name):
        py.ioff()
        py.figure(35,figsize=(32,16))
        py.clf()
    
        scale = 0.25

        #########  Pulsars  ############
        ax = py.subplot(2,4,1)
        ax.set_yscale("log", nonposy='clip')
        ax.set_xscale("log", nonposx='clip')
        py.title('Pulsars')
        names = []
        pts = []
        mi = 1e40
        ma = 0
        amask = self.areas>0
        for it,hist in enumerate(self.ponhists):
            p1 = py.errorbar(self.midpts,(hist)/self.areas,xerr=self.widths,marker='o',ls='None')
            p3 = py.errorbar(self.midpts,self.ponmodels[it]/self.areas,xerr=self.widths,yerr=np.sqrt(self.ponmodels[it])/self.areas,marker='o',ls='None')
            p2 = py.errorbar(self.midpts,(self.pofhists[it])/self.areas,xerr=self.widths,marker='o',ls='None')
            names.append(self.pulsars[it][0]+' ON')
            pts.append(p1[0])
            names.append(self.pulsars[it][0]+' model')
            pts.append(p3[0])
            names.append(self.pulsars[it][0]+' OFF')
            pts.append(p2[0])
            mi = min(mi,max((self.pofhists[it][amask])/self.areas[amask]))
            ma = max(ma,max((self.Npj[it]*self.psf[amask]+self.pulsars[it][1]*self.vij[it][amask])/self.areas[amask]))
        py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$dN/d\theta^{2}$')
        py.grid()
        mi = max(mi,1./max(self.areas[amask]))
        py.xlim(min(self.midpts-self.widths)*(1-scale),max(self.midpts+self.widths)*(1+scale))
        py.ylim(0.25*mi,2*ma)
        py.legend(pts,names,loc=3)

        #########  AGN plots  #############
        ax = py.subplot(2,4,2)
        ax.set_yscale("log", nonposy='clip')
        ax.set_xscale("log", nonposx='clip')
        py.title('AGN')
        names = []
        pts = []
        mi = 1e40
        ma = 0
        for it,hist in enumerate(self.agnhists):
            model = self.Naj[it]*self.psf+self.Ni[it]*self.iso + self.Nh[it]*self.hmd
            #print (self.Naj[it],self.psf,self.Ni[it],self.iso,self.Nh[0],self.hmd)
            modelerrs = np.sqrt((self.Naj[it]*self.psf/sum(self.psf)*np.sqrt((self.psfe/self.psf)**2+(self.Naje[it]/self.Naj[it])**2))**2+(self.Nie[it]*self.iso)**2+(self.Nhe[it]*self.hmd)**2)
            back = self.Ni[it]*self.iso
            backerrs = self.Nie[it]*self.iso
            p1 = py.errorbar(self.midpts,hist/self.areas,xerr=self.widths,marker='o',ls='None')
            p2 = py.errorbar(self.midpts,self.agnmodels[it]/self.areas,xerr=self.widths,yerr=np.sqrt(self.agnmodels[it])/self.areas,marker='o',ls='None')

            names.append(self.agnlist[it]+' Data')
            pts.append(p1[0])
            names.append(self.agnlist[it]+' Model')
            pts.append(p2[0])

            mi = min(mi,min((hist[amask])/self.areas[amask])*0.25)
            ma = max(ma,max((model[amask])/self.areas[amask]))

            if self.halomodel!='':
                #p4 = py.errorbar(self.midpts,self.Nh*self.hmd,xerr=self.widths,yerr=(self.Nhe*self.hmd)/self.areas,marker='o',ls='None')
                #names.append('Halo')
                #pts.append(p4[0])
                p5 = py.errorbar(self.midpts,(self.Nh[it]*self.hmd+self.Ni[it]*self.iso)/self.areas,xerr=self.widths,marker='o',ls='None')
                names.append('Halo+Iso')
                pts.append(p5[0])
            else:
                p3 = py.errorbar(self.midpts,back/self.areas,xerr=self.widths,yerr=np.sqrt(back)/self.areas,marker='o',ls='None')
                names.append(self.agnlist[it]+' Iso')
                pts.append(p3[0])
        py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$dN/d\theta^{2}$')
        py.xlim(min(self.midpts-self.widths)*(1-scale),max(self.midpts+self.widths)*(1+scale))
        if mi==1e40:
            mi = min(hist[amask]/self.areas[amask])
        py.ylim(max(0.25/max(self.areas),mi),2*ma)
        py.grid()
        py.legend(pts,names,loc=3)

        ############  PSF residuals plot  ###############
        ax = py.subplot(2,4,3)
        names = []
        pts = []
        ax.set_xscale("log", nonposx='clip')
        err = np.ones(self.nbins)
        p1 = py.errorbar(self.midpts,(self.psf-self.psfm)/self.psfm,xerr=self.widths,yerr=self.psfe/self.psfm,ls='None',marker='o')[0]
        cmask = self.psfm>0
        chisq = sum(((self.psf[cmask]-self.psfm[cmask])**2/(self.psfe[cmask])**2))
        
        py.grid()
        ma = max(abs((self.psf-self.psfm+self.psfe)/self.psfm))
        ma = max(ma,max(abs((self.psf-self.psfm-self.psfe)/self.psfm)))
        py.xlim(min(self.midpts-self.widths)*(1-scale),max(self.midpts+self.widths)*(1+scale))
        py.ylim(-1.5*ma,1.5*ma)
        py.legend([p1],['PSF estimator'])
        py.title('PSF residuals from %s'%self.irf)
        py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$(\rm{Data - Model})/Model$')
        py.figtext(0.54,0.6,'Chisq (dof) = %1.1f (%1.0f)'%(chisq,self.nbins))
        
        ##############  PSR Background estimators  ######
        ax = py.subplot(2,4,5)
        ax.set_xscale("log", nonposx='clip')
        names = []
        pts = []
        ma = 0.
        tchisq = 'Name: Chisq (%d)\n'%self.nbins
        for it,psr in enumerate(self.pofhists):
            try:
                pt1 = py.errorbar(self.midpts,(psr-self.vij[it])/self.vij[it],xerr=self.widths,yerr=1./np.sqrt(self.vij[it]),ls='None',marker='o')[0]
                names.append(self.pulsars[it][0])
                pts.append(pt1)
                mask = psr>0
                up = max((psr[mask]-self.vij[it][mask]+self.vije[it][mask])/self.vij[it][mask])
                down = min((psr[mask]-self.vij[it][mask]-self.vije[it][mask])/self.vij[it][mask])
                ma = max(ma,up)
                ma = max(ma,abs(down))
                cmask = self.vije[it]>0
                chisq = sum(((psr[cmask]-self.vij[it][cmask])**2/self.vij[it][cmask]))
                tchisq = tchisq + '%s: %1.1f (%d)\n'%(self.pulsars[it][0],chisq,len(psr[cmask]))
            except:
                print ('Bad plotting' )
        py.grid()
        py.title('PSR Background Estimator Residuals')
        py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$(\rm{Data - Model})/Model$')
        py.xlim(min(self.midpts-self.widths)*(1-scale),max(self.midpts+self.widths)*(1+scale))
        py.ylim(-1.5*ma,1.5*ma)
        py.figtext(0.15,0.15,tchisq)
        py.legend(pts,names)

        ##############  PSR model fractional differences  ######
        ax = py.subplot(2,4,6)
        ax.set_xscale("log", nonposx='clip')
        names = []
        pts = []
        ma = 0.
        tchisq='Name: Chisq (%d)\n'%self.nbins
        for it,psr in enumerate(self.ponhists):
            try:
                model = self.Npj[it]*self.psf + self.pulsars[it][1]*self.vij[it]
                modelerr = np.sqrt((self.Npj[it]*self.psfe)**2 + (self.Npje[it]*self.psf)**2 + (self.pulsars[it][1]*self.vije[it])**2)
                pt1 = py.errorbar(self.midpts,(psr-model)/model,xerr=self.widths,yerr=1./np.sqrt(model),ls='None',marker='o')[0]
                names.append(self.pulsars[it][0])
                pts.append(pt1)
                mask = psr>0
                up = max((psr[mask]-model[mask]+modelerr[mask])/model[mask])
                down = min((psr[mask]-model[mask]-modelerr[mask])/model[mask])
                ma = max(ma,up)
                ma = max(ma,abs(down))
                cmask = modelerr>0
                chisq = sum(((psr[cmask]-model[cmask])**2/model[cmask]))
                tchisq = tchisq + '%s: %1.1f (%d)\n'%(self.pulsars[it][0],chisq,len(psr[cmask]))
            except:
                print ('Bad plotting' )
        py.grid()
        py.title('PSR On-pulse Estimator Residuals')
        py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$(\rm{Data - Model})/Model$')
        py.xlim(min(self.midpts-self.widths)*(1-scale),max(self.midpts+self.widths)*(1+scale))
        py.ylim(-1.5*ma,1.5*ma)
        py.figtext(0.35,0.15,tchisq)
        py.legend(pts,names)

        ##############  AGN model fractional differences  ######
        ax = py.subplot(2,4,7)
        ax.set_xscale("log", nonposx='clip')
        names = []
        pts = []
        ma = 0.
        tchisq='Name: Chisq (%d)\n'%self.nbins
        for it,agn in enumerate(self.agnhists):
            try:
                model = self.Naj[it]*self.psf + self.Ni[it]*self.iso + self.Nh[it]*self.hmd
                modelerr = np.sqrt((self.Naj[it]*self.psfe)**2 + (self.Naje[it]*self.psf)**2 + (self.Nie[it]*self.iso)**2 + (self.Nhe[it]*self.hmd)**2)
                pt1 = py.errorbar(self.midpts,(agn-model)/model,xerr=self.widths,yerr=1./np.sqrt(model),ls='None',marker='o')[0]
                names.append(self.agnlist[it])
                pts.append(pt1)
                mask = agn>0
                up = max((agn[mask]-model[mask]+modelerr[mask])/model[mask])
                down = min((agn[mask]-model[mask]-modelerr[mask])/model[mask])
                ma = max(ma,up)
                ma = max(ma,abs(down))
                cmask = modelerr>0
                chisq = sum(((agn[cmask]-model[cmask])**2/model[cmask]))
                tchisq = tchisq + '%s: %1.1f (%d)\n'%(self.agnlist[it],chisq,len(agn[cmask]))
            except:
                print ('Bad plotting' )
        py.grid()
        py.title('AGN Model Estimator Residuals')
        py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$(\rm{Data - Model})/Model$')
        py.xlim(min(self.midpts-self.widths)*(1-scale),max(self.midpts+self.widths)*(1+scale))
        py.ylim(-1.5*ma,1.5*ma)
        py.legend(pts,names)
        py.figtext(0.54,0.15,tchisq)
        py.figtext(0.72,0.35,str(self),fontsize=18)
        py.savefig(name+'.png')

############################################   Helper functions  ##############################################################


    ######################################################################
    #    Determination of maximum likelihood estimator of vij from (3)   #
    ######################################################################
    ## the maximum likelhood estimator of the pulsar background
    # @param a ratio of on/off phase window             (observation)
    # @param n number of photons in on window bin       (observation)
    # @param b bumber of photons in off window bin      (observation)
    # @param m PSF in bin                               (derived)
    # @param N number of photons associated with pulsar (derived)
    def backest(self,a,n,b,m,N):

        sterm = 4*a*(1+a)*b*m*N+(m*N-a*(b+n-m*N))**2    #discriminant
        #catch negative discriminant
        if sterm<0.:
            print ('Unphysical Solution: %1.4f'%sterm)

        #calculate background estimator analytically
        v = a*(b+n)-m*N-a*m*N+np.sqrt(sterm)
        v = v/(2.*a*(1+a))
        return v

    ######################################################################
    #    Gradient of maximum likelihood estimator of vij from (3)        #
    ######################################################################
    ## the maximum likelhood estimator of the pulsar background
    # @param a ratio of on/off phase window             (observation)
    # @param n number of photons in on window bin       (observation)
    # @param b bumber of photons in off window bin      (observation)
    # @param m PSF in bin                               (derived)
    # @param N number of photons associated with pulsar (derived)
    def gradback(self,a,n,b,m,N):
        sterm = 4*a*(1+a)*b*m*N+(m*N-a*(b+n-m*N))**2    #discriminant
        #catch negative discriminant
        if sterm<0.:
            print ('Unphysical Solution: %1.4f'%sterm)

        #calculate gradient of background estimator analytically
        grad1 = -N-a*N+(4*a*(1+a)*b*N+2*(m*N-a*(b+n-m*N))*(N-a*(-N)))/(2*np.sqrt(sterm))    #psf derivative
        grad2 = -m-a*m+(4*a*(1+a)*b*m+2*(m*N-a*(b+n-m*N))*(m-a*(-m)))/(2*np.sqrt(sterm))    #number estimator derivative
        v = np.array([grad1,grad2])
        v = v/(2.*a*(1+a))
        return v

    def loggrad(self,params):
        return self.gradient(params)*np.array(params)*np.log(10.)

    ######################################################################
    #         Gradient of maximum likelihood from (3)                    #
    ######################################################################
    ## The gradient of the maximum likelhood estimator
    # @param params likelihood parameters defined by likelihood function
    # @param verb verbose slow output of gradient calculation
    def gradient(self,params,verb=False):
        #verb=self.verbose
        params = 10**params
        #limits checks
        for it,param in enumerate(params):
            if self.fixed[it]:
                params[it]=self.params[it]
            if param>self.limits[it][1]:
                params[it]=self.limits[it][1]
            if param<self.limits[it][0]:
                params[it]=self.limits[it][0]
            #print ('Failed limit check')
            #return INFINITY
            if self.fixed[it]:
                params[it]=self.params[it]

        psrs = len(self.ponhists)
        agns = len(self.agnhists)

        npij = self.ponhists
        bpij = self.pofhists
        naij = self.agnhists

        nmu = self.nbins - 1

        mi = params[:nmu]
        Npj = params[nmu:nmu+psrs]
        Naj = params[nmu+psrs:nmu+psrs+agns]
        Ni = params[nmu+psrs+agns:nmu+psrs+agns+agns]
        Nh = params[nmu+psrs+agns+agns:]

        mun = 1-np.sum(mi)
        mi = np.append(mi,[mun])
        grad = []
        if verb:
            print ('Obs\tmod\tfact\tnum\tacc')
            print ('----------------------')

        #PSF gradient
        for it in range(nmu):
            acc = 0
            flag = False
            
            #loop over pulsars
            for it2 in range(psrs):
                alpha = self.pulsars[it2][1]                                                                                
                backg = self.backest(alpha,npij[it2][it],bpij[it2][it],mi[it],Npj[it2])                                     #vij
                denom = Npj[it2]*mi[it]+alpha*backg                                                                         #ith model
                grad0 = self.gradback(alpha,npij[it2][it],bpij[it2][it],mi[it],Npj[it2])[0]                                 #d(vij)/d(mu_i)
                estn = self.backest(alpha,npij[it2][nmu],bpij[it2][nmu],mi[nmu],Npj[it2])                                   #vnj
                gradn = self.gradback(alpha,npij[it2][nmu],bpij[it2][nmu],mi[nmu],Npj[it2])[0]                              #d(vnj)/d(mu_n)
                denomn = Npj[it2]*mi[nmu]+alpha*estn                                                                        #nth model
                fact0 = 0 if npij[it2][it]==0 else (Npj[it2]+alpha*grad0)*(npij[it2][it]/(denom) - 1.)                      #ith signal contribution
                fact1 = 0 if npij[it2][nmu]==0 else -(Npj[it2]+alpha*gradn)*(npij[it2][nmu]/denomn-1.)                      #nth signal contribution
                fact2 = 0 if bpij[it2][it]==0 else grad0*(bpij[it2][it]/backg-1.)                                           #ith background contribution
                
                #catch bad gradients
                if (denom <=0 and npij[it2][it]>0) or (backg<=0 and bpij[it2][it]>0) or (denomn<=0 and npij[it2][nmu]>0):
                    flag=True
                acc = acc + fact0 + fact1 + fact2
                if verb:
                    print (npij[it2][it],denom,npij[it2][nmu],denomn,bpij[it2][it],backg,acc)
                    t.sleep(0.25)

            #loop over AGN
            for it2 in range(agns):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh[it2]*self.hmd[it]                                             #ith model
                denomn = Naj[it2]*mi[nmu]+Ni[it2]*self.iso[nmu]+Nh[it2]*self.hmd[nmu]                                         #nth model
                fact0 = 0 if naij[it2][it]==0 else Naj[it2]*(naij[it2][it]/(denom) - 1.)                                    #ith contribution
                fact1 = 0 if naij[it2][nmu]==0 else -Naj[it2]*(naij[it2][nmu]/denomn - 1.)                                  #nth contribution

                #catch bad gradients
                if (denom <=0 and naij[it2][it]>0) or (denomn<=0 and naij[it2][nmu]>0):
                    flag=True
                acc = acc + fact0 + fact1
                if verb:
                    print (naij[it2][it],denom,naij[it2][nmu],denomn,acc)
                    t.sleep(0.25)

            if flag:
                grad.append(0)
            else:
                grad.append(-acc)
            if verb:
                print ('----------------------')

        #Pulsar Number estimator gradient
        for it2 in range(psrs):
            alpha = self.pulsars[it2][1]
            flag = False
            acc = 0
            for it in range(self.nbins):
                denom = Npj[it2]*mi[it]+alpha*self.backest(alpha,npij[it2][it],bpij[it2][it],mi[it],Npj[it2])               #ith model
                grad1 = self.gradback(alpha,npij[it2][it],bpij[it2][it],mi[it],Npj[it2])[1]                                 #d(vij)/d(Npj)
                
                fact = 0 if npij[it2][it]==0 else (mi[it]+alpha*grad1)*(npij[it2][it]/denom - 1.)

                if (denom <=0 and npij[it2][it]>0):
                    flag=True
                acc = acc + fact
                if verb:
                    print (npij[it2][it],denom,mi[it],alpha*grad1,acc)
                    t.sleep(0.25)
            if flag:
                grad.append(0)
            else:
                grad.append(-acc)
            if verb:
                print ('----------------------')

        #AGN number estimator gradient
        for it2 in range(agns):
            acc = 0
            flag = False
            for it in range(self.nbins):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh[it2]*self.hmd[it]
                fact = 0 if naij[it2][it]==0 else mi[it]*(naij[it2][it]/denom - 1.)
                if (denom <=0 and naij[it2][it]>0):
                    flag=True
                acc = acc + fact
                if verb:
                    print (naij[it2][it],denom,acc)
                    t.sleep(0.25)
            if flag:
                grad.append(0)
            else:
                grad.append(-acc)
            if verb:
                print ('----------------------')

        #Isotropic number estimator gradient for AGN
        for it2 in range(agns):
            acc = 0
            flag = False
            for it in range(self.nbins):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh[it2]*self.hmd[it]
                fact = 0 if naij[it2][it]==0 else self.iso[it]*(naij[it2][it]/denom - 1.)
                if (denom <=0 and naij[it2][it]>0):
                    flag=True
                acc = acc + fact
                if verb:
                    print (naij[it2][it],denom,acc)
                    t.sleep(0.25)
            if flag:
                grad.append(0)
            else:
                grad.append(-acc)
            if verb:
                print ('----------------------')
        

        #Halo number estimator gradient for AGN
        flag = False
        for it2 in range(agns):
            acc = 0
            for it in range(self.nbins):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh[it2]*self.hmd[it]
                fact = 0 if naij[it2][it]==0 else self.hmd[it]*(naij[it2][it]/denom - 1.)
                if (denom <=0 and naij[it2][it]>0):
                    flag=True
                acc = acc + fact
                if verb:
                    print (naij[it2][it],denom,acc)
                    t.sleep(0.25)
            if flag:
                grad.append(0)
            else:
                grad.append(-acc)
            if verb:
                print ('----------------------')
        #for it,grd in enumerate(grad):
        #    print (grd,params[it],grd*params[it])
        return np.array(grad)

    ######################################################################
    #        Simple Error calculation from likelihood (1D)               #
    ######################################################################
    ## Error calculation about the maximum for one parameter (no covariance)
    # @param num number of parameter to find error
    def errors(self,num):
        eig = np.zeros(len(self.minuit.params))
        eig[num]=1.
        if self.mode>0:
            disp = 1
        else:
            disp = 0
        
        #find points on either side of maximum where likelihood has decreased by 1/2
        err1 = so.fmin_powell(lambda x: abs(self.likelihood(np.log10(self.params+x[0]*eig))-self.lmax-0.5),
                              [self.minuit.params[num]*0.01],full_output=1,disp=disp)
        err1 = abs(err1[0])
        err2 = so.fmin_powell(lambda x: abs(self.likelihood(np.log10(self.params-x[0]*eig))-self.lmax-0.5),
                              [self.minuit.params[num]*0.01],full_output=1,disp=disp)
        err2 = abs(err2[0])

        #try to catch badly formed likelihood surfaces
        if self.likelihood(np.log10(self.params-err1*eig))==INFINITY:
            return err2*err2
        if self.likelihood(np.log10(self.minuit.params+err2*eig))==INFINITY:
            return err1*err1
        return err1*err2

    ######################################################################
    #      Covariant error estimation from the likelihood surface        #
    ######################################################################
    ## Estimates the true error by projecting the surface in multiple dimensions
    # @param num number of parameter to find error
    def finderrs(self,num):
        ma = 0.
        err = np.sqrt(self.errs2[num])
        rt = np.sqrt(2)

        #py.figure(2,figsize=(16,16))
        #py.clf()
        rows = int(np.sqrt(len(self.params)))+1
        #go through all parameters
        #find the quadratic form of the likelihood surface
        #determine the increase in the variance from covariance of parameters
        for it in range(len(self.params)):

            if it!=num:
                #py.subplot(rows,rows,it+1)
                eigx = np.zeros(len(self.params))
                eigx[num]=err
                eigy = np.zeros(len(self.params))
                eigy[it]=np.sqrt(self.errs2[it])

                #calulate likelihood along ring around maximum (1-sigma in likelihood)
                px = (self.likelihood(np.log10(self.params+eigx))-self.lmax)
                pxpy = (self.likelihood(np.log10(self.params+(eigx+eigy)/rt))-self.lmax)
                py1 = (self.likelihood(np.log10(self.params+eigy))-self.lmax)
                mxpy = (self.likelihood(np.log10(self.params+(-eigx+eigy)/rt))-self.lmax)
                mx = (self.likelihood(np.log10(self.params-eigx))-self.lmax)
                mxmy = (self.likelihood(np.log10(self.params+(-eigx-eigy)/rt))-self.lmax)
                my = (self.likelihood(np.log10(self.params-eigy))-self.lmax)
                pxmy = (self.likelihood(np.log10(self.params+(eigx-eigy)/rt))-self.lmax)
                q = [px,pxpy,py1,mxpy,mx,mxmy,my,pxmy]
                
                """gridpts = 12
                minmax = 5.
                vals = np.arange(-gridpts,gridpts+1,1)
                vals = vals*minmax/gridpts
                z = np.array([self.likelihood(self.minuit.params+eigx*x+eigy*y)-self.minuit.fval for x in vals for y in vals]).reshape(len(vals),-1)
                py.contour(vals,vals,z,[0.5,2.,4.5,8.,12.5])
                #py.colorbar()
                py.xlim(-minmax,minmax)
                py.ylim(-minmax,minmax)"""
                #find quadratic fit to likelihood surface
                try:
                    el = Ellipse(q)
                    """x,y=el.contour(1/np.sqrt(2))
                    py.plot(x,y,'-')
                    ma = max(max(x),max(y))*1.5
                    py.xlim(min(-2,-ma),max(2,ma))
                    py.ylim(min(-2,-ma),max(2,ma))"""


                    pars = el.qf.p
                    a,b,c = (pars[0]-pars[4]*pars[4]/(4.*pars[2])),pars[1],(-0.5-pars[3]*pars[3]/(4.*pars[2]))

                    #parameter we're interested in is x-axis
                    #find the x-value tangent at quadratic surface = 0.5
                    xmin = abs(-b/(2*a)+np.sqrt(b**2-4*a*c)/(2*a))
                    ymin = -(pars[3]+pars[4]*xmin)/(2*pars[2])
                    terr = abs(xmin)

                    #is it the largest?
                    ma = max(ma,terr)

                    #print (xmin,ymin)
                    #t.sleep(0.25)
                except:
                    pass
                    #print ('Caught poorly formed quad surface')
        #py.savefig('likes%d.png'%num)
        return (ma*err)**2

    def tspval(self):
        pval = 1.-sst.stats.chisqprob(max(self.TS,0),1)/2.
        bestp = so.fmin_powell(lambda x: (spec.erf(x[0]/np.sqrt(2))-pval)**2,[1.],full_output=1,disp=0)
        return [1.-pval,bestp[0]]


############################################   End of CombinedLike Class ######################################################




############ unit test ##########################
## test function
# @param bins number of adaptive bins
# @param ctype conversion type 0:front, 1:back
# @param emin minimum energy
# @param emax maximum energy
# @param days number of days of data to examine from start of P6 data
def test(bins=12,ctype=0,emin=1000,emax=1778,days=730,irf='P6_v3_diff',maxr=-1,sel='[0:2]',agnlis=['agn-psf-study-bright'],double=False,ctlim=[0.4,1.0],tol=1e-1,verbose=False):
    psf = CALDBPsf(CALDBManager(irf=irf))
    ebar = np.sqrt(emin*emax)
    psrs = ''
    if maxr<0:
        maxr = psf.inverse_integral(ebar,ctype,99.5)*1.5             #use 1.5 times the 99.5% containment as the maximum distance
    cl = CombinedLike(irf=irf,mode=-1,pulsars = eval('pulsars'+sel+''),agnlist=agnlis,verbose=verbose,ctmin=ctlim[0],ctmax=ctlim[1])
    cl.loadphotons(0,maxr,emin,emax,239557417,239517417+days*86400,ctype)
    cl.bindata(bins)
    cl.tol=tol
    f0 = cl.fit()
    params = cp.copy(cl.fitpsf(double))
    print (str(cl))
    for psr in cl.pulsars:
        psrs = psrs + '%s_'%psr[0]
    #cl.makeplot('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%s%s'%(emin,emax,ctype,maxr,bins,psrs,('').join(cl.agnlist)))
    return params,cl,f0

################# halo analysis ###################
## halo runs the comparison between pulsars and AGN
#  @param bins number of angular bins
#  @param ctype conversion type 0-front,1-back
#  @param emin minimum energy
#  @param emax maximum energy
#  @param days number of days since MET start to examine
#  @param irf reference response function to start with
#  @param maxr maxROI in degrees, -1 will choose based on IRF
#  @param sel string to select subset of pulsars, default is all [vela,geminga]
#  @param agnlis list of AGN to compare, if more than one specified, the last one will be examined for halos and the others will be calibration sources
#  @param model angular model for halo to check ['CDisk','CHalo'], see uw.stacklike.angularmodels for more information
def halo(bins=12,ctype=0,emin=1000,emax=1778,days=fdays,irf='P6_v3_diff',maxr=-1,sel='[0:2]',agnlis=['agn-psf-study-bright'],model='CDisk',verbose=False,veryverbose=False,ret=False):
    import cPickle
    #setup the binned likelihood object
    psf = CALDBPsf(CALDBManager(irf=irf))
    ebar = np.sqrt(emin*emax)
    if maxr<0:
        maxr = psf.inverse_integral(ebar,ctype,99.5)*1.5                       #use 1.5 times the 99.5% containment as the maximum distance
    #test different sizes of halos for the current likelihood
    #pts = np.arange(0,nps+1,1)
    pts = [0.1,0.5,1.0]
    for pt in pts:
        cl = CombinedLike(irf=irf,mode=-1,pulsars = eval('pulsars'+sel+''),agnlist=agnlis,verbose=verbose,veryverbose=veryverbose,qandd=True)
        cl.loadphotons(0,maxr,emin,emax,239557417,239517417+days*86400,ctype)
        cl.bindata(bins)#,halomod=model,par=[pt/rd,ebar,ctype])
        cl.fit(qandd=True)
        npsrs = len(cl.ponhists)
        nagns = len(cl.agnhists)
        flag1 = agnlis[-1]=='1es0229p200' or agnlis[-1]=='1es0347-121'
        if flag1:
            cl.params[cl.nbins+npsrs+nagns-2]=ZERO  #allpars[it+npsrs+len(psfm)]
            cl.fixed[cl.nbins+npsrs+nagns-2]=True
        print (cl.fixed)
        if ret:
            return cl

        #null likelihood - no halos!
        f0 = cl.fit(custom=flag1,qandd=True)
        print (cl)

        #find best single king fit parameters
        #fitpars = cl.fitpsf()
        psrs = ''

        nmu = cl.nbins-1

        ip = nmu+npsrs+nagns

        #find the reference number of AGN photons for no halo
        cl.psfphotons = cl.Naj[-1]
        psfphotons = cl.Naj[-1]
        cl.plflux = cl.prof[ip-1][0]
        cl.puflux = cl.prof[ip-1][1]

        #find the PSF parameters 
        #bestpsf =so.fmin_powell(lambda x: cl.profile(ip-1,x[0]),[cl.Naj[-1]],disp=0,full_output=1)
        #print (cl)
        print (cl.plflux,cl.psfphotons,cl.puflux)

        for psr in cl.pulsars:
            psrs = psrs + '%s_'%psr[0]
        agns  = ''
        for agn in cl.agnlist:
            agns = agns + '%s_'%agn

        r05 = psf.inverse_integral(ebar,ctype,5.)
        r34 = psf.inverse_integral(ebar,ctype,34.)
        r68 = psf.inverse_integral(ebar,ctype,68.)
        r95 = psf.inverse_integral(ebar,ctype,95.)
        r99 = psf.inverse_integral(ebar,ctype,99.)
        testpsf = PSF(lims=[0,1],model_par=[0.05,2.25])
        print (r68,r95)
        testpsf.fromcontain([r68,r95],[0.68,0.95])
        fitpars = cl.fitpsf() if ebar < 10000 else [cl.ebar,testpsf.model_par[0],testpsf.model_par[1]]
        print (fitpars)
        nps = 10

        npars = len(cl.params)
        agnback = cp.copy(cl.agnhists[-1])
        uplims2 = []
        uplims1 = []
        detect = []

        of = open('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/uplimemi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%s%s_%s.txt'%(emin,emax,ctype,maxr,bins,psrs,agns,model),'w')
        print ('%1.3f'%f0, file=of)
        print ('%1.1f'%psfphotons, file=of)
        fracs = np.linspace(0.05,0.95,19)
        #print (string.join(['%1.2f'%(frac) for frac in fracs],'\t'), file=of)
        print ('#size(rad)\tmaxL\tNha \tcl68\tcl95\tTS\tNpsf', file=of)

        if flag1:
            cl.setuppars(halomodel=model,haloparams=[pt/rd,ebar,ctype,fitpars[1]/rd,fitpars[2]])
            cl.params[cl.nbins+npsrs+nagns-2]=ZERO#allpars[it+npsrs+len(psfm)]
            cl.fixed[cl.nbins+npsrs+nagns-2]=True
            cl.params[-1]=(cl.Naj[-1]+cl.Ni[-1])*0.5
            cl.fixed[-1]=False
            print (cl.fixed)
        else:
            cl.setuppars(halomodel=model,haloparams=[pt/rd,ebar,ctype,fitpars[1]/rd,fitpars[2]])
            cl.params[cl.nbins+npsrs-1]=cl.Naj[-1]/2.#allpars[it+npsrs+len(psfm)]
            cl.fixed[cl.nbins+npsrs-1]=False
            cl.params[-1]=cl.Naj[-1]/2.
            cl.fixed[-1]=False
        #find best parameters for the fit to the halo model 
        cl.fit(custom=True,qandd=True)
        print (cl)
        print (cl.fixed)
        lmax = cl.lmax

        #best = so.fmin_powell(lambda x:cl.profile(npars-1,x),cl.Nh[-1],disp=0,full_output=1,ftol=abs(0.001/cl.lmax))
        #cl.lmax=best[1].item()
        #cl.Nh[-1]=best[0].item()
        #cl.profile(npars-1,cl.Nh[-1])
        #cl.params=cl.tpars
        #cl.setmembers()
        #get the TS of the halo
        #cl.TS = 2*(cl.profile(npars-1,ZERO)-cl.profile(npars-1,cl.Nh[-1]))
        TSs = []
        TSs.append(cl.TS)

        #find the 5 sigma upper limit
        if ebar<10000:
            bestup = cl.findlims(npars-1,12.5)  #,step)#so.fmin_powell(lambda x: minup(x[0],cl,lmax,npars,cl.Nh[-1]),[cl.Nh[-1]+step],disp=0,full_output=1,maxiter=2)
        else:
            bestup = so.fmin_powell(lambda x: (cl.profile(npars-1,x,False)-cl.lmax-12.5)**2 if x>cl.Nh[-1] else -cl.lmax,(cl.Nh[-1]+cl.Naj[-1]),disp=0,ftol=0.1,full_output=1)
            bestup = [bestup[0].item(),bestup[0].item()]
        print (bestup[1])
        bestup = [bestup[1],cl.profile(npars-1,bestup[1])-cl.lmax]

        def printstate():
            print ('New iteration')
        
        assert cl.lmax!=np.Infinity

        usig = so.fmin_powell(lambda x: (cl.profile(npars-1,x,False)-cl.lmax-0.5)**2 if x>cl.Nh[-1] else -cl.lmax,(cl.Nh[-1]+cl.Naj[-1]),disp=0,ftol=0.1,full_output=1)
        delt = usig[0].item()-cl.Nh[-1]
        lsig =  so.fmin_powell(lambda x: abs(cl.profile(npars-1,x,False)-cl.lmax-0.5)**2 if x<cl.Nh[-1] else -cl.lmax, max(ZERO,cl.Nh[-1]-delt),disp=0,ftol=0.1,full_output=1)
        #profile the likelihood near the maximum between +/- 5 sigma, if Nhalo is negative, stop at 0
        xr = np.linspace(ZERO,bestup[0],20)
        #if bestup[0]<cl.prof[npars-1][1]:
        #    xr = np.linspace(ZERO,5*cl.prof[npars-1][1]-4*cl.params[npars-1],20)
        fx,frac,tpars = [],[],[]
        for x in xr:
            fx.append(lmax-cl.profile(npars-1,x))
            frac.append(x/(cl.tpars[ip-1]+x))
            tpars.append(cl.tpars)
        
        tpars = np.array(tpars).transpose()
        py.figure(45,figsize=(16,16))
        py.clf()
        py.subplot(2,1,1)
        for it in range(cl.nbins-1):
            py.plot(xr,tpars[it])
        py.subplot(2,2,3)
        for it in range(npsrs):
            py.plot(xr,tpars[cl.nbins-1+it])
        py.subplot(2,2,4)
        for it in range(3*nagns):
            py.plot(xr,tpars[cl.nbins-1+it+npsrs])
        

        py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s_%s_profparams.png'%(emin,emax,ctype,maxr,bins,pt,psrs,agns,model,irf))
        it=0
        fx,frac = np.array(fx),np.array(frac)

        #np.save('xr.npy',xr)
        #np.save('fx.npy',fx)
        #np.save('frac.npy',frac)
        #max bad fits
        mask = np.array([not np.isnan(fxr) for fxr in fx])
        xr,fx=xr[mask],fx[mask]

        print (xr)
        print (frac)
        print (fx)
        fx2 = cp.copy(fx)

        #calculate the 68 and 95% upperlimits
        cl1 = BayesianLimit(xr,fx).getLimit(0.32)
        cl2 = BayesianLimit(xr,fx).getLimit(0.05)

        mask = frac<0.90
        pos = len(frac[mask])-1
        xr2 = np.linspace(ZERO,xr[pos],20)
        fx,frac = [],[]
        for x in xr2:
            fx.append(lmax-cl.profile(npars-1,x))
            frac.append(x/(cl.tpars[ip-1]+x))
        fx,frac = np.array(fx),np.array(frac)

        print (frac)
        print (fx)
        clfrac = BayesianLimit(frac,fx).getLimit(0.05)

        #output Nhalo, and the two upper limits
        detect.append(cl.Nh[-1]/(psfphotons))
        print (cl2)
        uplims2.append(cl2/(psfphotons))
        uplims1.append(cl1/(psfphotons))
        #print (cl2[0],psfphotons,cl2[0]/(psfphotons))
        cl.ul68 = cl1
        cl.ul95 = cl2
        cl.Nhe[-1]= (bestup[0]-cl.Nh[-1])/5.
        cl.lflux = lsig[0].item()
        cl.uflux = usig[0].item()
        realnum=100
        if cl.lflux<0:
            cl.pval=0.5
        else:
            suc = 0#cl.realize(realnum)
            cl.pval=0.5#suc/(realnum*1.)
            #cl.pvall,cl.pvalu=cl.wilson(suc,realnum)
        if (cl.Nh[-1]+cl.Naj[-1])<1.:
            xrg = np.arange(ZERO,cl.Ni[-1],25)
        else:
            xrg = np.arange(ZERO,cl.Nh[-1]+cl.Naj[-1],25)
        #ypts = np.array([cl.lmax-cl.profile(npars-1,xr) for xr in xrg])
        py.figure(25,figsize=(8,8))
        py.clf()
        py.plot(xr,fx2)
        py.plot(cl.Nh[-1],0,'rd')
        py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s_%s_prof.png'%(emin,emax,ctype,maxr,bins,pt,psrs,agns,model,irf))
        py.clf()
        py.plot(frac,fx,'b+')
        py.plot(max(cl.Nh[-1]/(cl.Nh[-1]+cl.Naj[-1]),0),0,'rd')
        py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s_%s_proffrac.png'%(emin,emax,ctype,maxr,bins,pt,psrs,agns,model,irf))
        print ('Lflux:',cl.lflux)
        print ('Nhalo:',cl.Nh[-1])
        print ('Uflux:',cl.uflux)
        print ('U95frac:',clfrac)
        print ('U95flux:',cl.ul95)
        print ('Not Sure:',cl.uflux/(cl.Naj[-1]+cl.uflux))
        print ('P-val halo:',cl.pval)
        del cl.minuit
        cl.frac=clfrac
        cl.makeplot('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s_%s'%(emin,emax,ctype,maxr,bins,pt,psrs,agns,model,irf))
        pfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s_%s.pickle'%(emin,emax,ctype,maxr,bins,pt,psrs,agns,model,irf),'w')
        cPickle.dump(cl,pfile)
        pfile.close()
        print ('%1.5f\t%1.3f\t%1.4f\t%1.4f\t%1.4f\t%1.3f'%(pt/rd,lmax,cl.Nh[-1],cl1,cl2,cl.TS), file=of)

    #make some summary plots

    of.close()

def gettask(agnbins=12,ctype=0,emin=1000,emax=3162,lis='agn_redshift2_lo_bzb',sel='[0:2]',mod='CHalo',days=fdays,irf='P6_v11_diff',verbose=False,veryverbose=False,ret=False,suplist=''):
    agnlists = [[],[],[],['agn_redshift2_lo_bzb']]
    return 'ub.halo(bins=%d,ctype=%d,emin=%d,emax=%d,days=%d,irf=\'%s\',maxr=4,sel=\'%s\',agnlis=[%s\'%s\'],model=\'%s\',verbose=%s,veryverbose=%s,ret=%s)'%(agnbins,ctype,emin,emax,days,irf,sel,suplist,lis,mod,str(verbose),str(veryverbose),str(ret))

def makealltasks():
    model=['CDisk','CHalo']
    emins = [1000,3162,10000,31623]#100,178,316,562,
    emaxs = [3162,10000,31623,100000]#178,316,562,1000,
    sel = ['[0:2]','[0:2]','[0:2]','[0:0]']
    suplists = ['','','','\'agnlobzb\',']
    tasks = [gettask(16,0,emins[x],emaxs[x],lis,sel=sel[x],mod=mod,irf=irf,days=fdays,suplist=suplists[x]) 
        for x in range(len(emins)) for lis in ['agnlobzb','agnhibzb','1es0229','1es0347'] for mod in model]#,'PKS0754p100','agn_redshift2_lo','agn_redshift2_hi'
    return tasks

## runallhalos - runs the halo analysis on multiple engines
#  @param model array of models to test
def runallhalos():
    #machines = 'tev01 tev02 tev03 tev04 tev05 tev06 tev07 tev08 tev09 tev10 tev11'.split()#
    #cluster.start_master()
    #os.system('ipcluster start&')
    #import time as t
    #t.sleep(60)
    from IPython.parallel import Client
    setup_string = 'import uw.stacklike.binned as ub;reload(ub);from uw.stacklike.binned import *'
    clfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/allhalos.txt','w')
    clfile.close()
    
    tasks = makealltasks()
    #'agn_redshift2_lo_bzb','agn_redshift2_lo_bzb\',\'agn_redshift2_hi_bzb','agn_redshift2_lo_bzb\',\'1es0229p200','agn_redshift2_lo_bzb\',\'1es0347-121','agn_redshift2_lo_bzb\',\'agn_redshift2_hi_bzb'
    """engines = 8#len(tasks)/len(machines)
    ua.setup_mec(engines=engines,machines=machines,clusterfile='/phys/groups/tev/scratch1/users/Fermi/mar0/allhalos.txt',clobber=True)
    t.sleep(30)
    
    at = ua.AssignTasks(setup_string,tasks,log=logfile,timelimit=1000000,progress_bar=True,ignore_exception=True)
    at(30)
    ua.kill_mec()"""
    
    rc = Client()
    dview = rc[:]

    logfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/python/mec.log','w')
    #engine = engines.Engines(log=logfile,quiet=False,maxerrors=len(tasks))
    dview.execute('import uw.stacklike.binned as ub')
    dview.execute('reload(ub)')
    dview.execute('tasks = ub.makealltasks()')
    dview.scatter('nums',range(len(tasks)))
    dview.execute('x=[eval(tasks[y]) for y in nums]')
#    engine.execute(setup_string)
#    for task in tasks:
#        engine.submit(task,[0])



def makeallplots():
    machines = 'tev01 tev02 tev03 tev04 tev05 tev06 tev07 tev08 tev09 tev10 tev11'.split()#
    setup_string = 'import uw.stacklike.binned as ub;reload(ub);from uw.stacklike.binned import *'
    tasks = ['ub.makeplots(%s,True,\'agn_redshift2_lo_bzb\',True)'%(0,emins[x],emaxs[x],lis,mod) for x in range(len(emins)) for lis in ['\'agn_redshift2_hi_bzb\'','\'1es0229p200\'','\'1es0347-121\'','\'PKS0754p100\''] for mod in model]
    print (tasks)
    engines = 4#len(tasks)/len(machines)
    ua.setup_mec(engines=engines,machines=machines,clobber=True)
    t.sleep(30)
    logfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/python/mec.log','w')
    at = ua.AssignTasks(setup_string,tasks,log=logfile,timelimit=10000,progress_bar=True,ignore_exception=True)
    at(30)
    ua.kill_mec()


## makeplots - take output information from the halo() analysis and makes the upper limit plots
#  @param iname name of sourcelist where halos were examined, ie 'agn-psf-study-bright'
#  @param exp calculate halo components with respect to exposure, else fraction of PSF photons
#  @param ext if there was an agnlist used as an additional calibration, needs to be added here
#  @param useupper only plot upper limits
def makeplots(iname,uexp=False,ext='',useupper=False):
    import cPickle
    import matplotlib as mpl
    mpl.rcParams['font.size']=16
    mpl.rcParams['font.family']='serif'
    psf = CALDBPsf(CALDBManager(irf='P6_v11_diff'))

    #get list of sources and setup exposure calculation if desired
    if uexp:
        s.EffectiveArea.set_CALDB(os.environ['CALDB']+'/data/glast/lat')
        ea = s.EffectiveArea('P6_v11_diff_front')
        if anname=='PASS6':
            ltc = s.LivetimeCube('/phys/groups/tev/scratch1/users/Fermi/mar0/data/pulsar/psr-livetime.fits')
        else:
            ltc = s.LivetimeCube('/phys/groups/tev/scratch1/users/Fermi/mar0/data/7.3src/all-lt.fits')
        exp = s.Exposure(ltc,ea)
        slist = file('/phys/users/mar0/sourcelists/%s.txt'%iname)
        header = slist.readline()
        sdlist = [s.SkyDir(float(line.split()[1]),float(line.split()[2])) for line in slist]
    else:
        exp=1

    #set up figures and energy binning for plots
    energy = [1000,1778,3162,5623,10000,17783,31623,1e5]
    energy2 = [1000,3162,10000,31623,1e5]
    prepat = ['vela_gem_','vela_gem_','vela_gem_','agn_redshift2_lo_bzb_']
    fig = py.figure(9,figsize=(7,6))
    py.clf()
    ctype=0
    name = ext+'_'+iname+'_' if ext!='' else iname+'_'
    width = [0.1,0.5,1.0]
    model= ['Disk','Halo']
    model2 = ['Disk','Gaussian']

    def marker(cl,exp=1.,useexp=False,useupper=False,ctype=0):
        xl, xh = cl.emin,cl.emax#cl.emin*(10**(1/16.)), cl.emax/(10**(1/16.))
        lflux = cl.lflux
        uflux = cl.uflux
        cn = cl.Nh[-1]
        #cl.fit(halomodel='')
        yc = cl.psfphotons
        yl= cl.plflux
        yh = cl.puflux
        yl = lflux/exp if useexp else lflux
        yh = uflux/exp if useexp else uflux
        yc = cn/exp if useexp else cn

        bc = cl.ebar

        if ( cl.lflux>1.) and not useupper:
            if not useexp:
                nmu = cl.nbins-1
                npsrs = len(cl.ponhists)
                nagns = len(cl.agnhists)
                ip = nmu+npsrs+nagns
                hn = len(cl.params)-1
                cl.profile(hn,yh)
                yh = yh/(cl.tpars[ip-1]+yh)
                cl.profile(hn,yl)
                yl = yl/(cl.tpars[ip-1]+yl)
                cl.profile(hn,yc)
                yc = yc/(cl.tpars[ip-1]+yc)
            print (yl,yc,yh,cl.psfphotons/exp,cl.plflux/exp,cl.puflux/exp,cl.lflux,cl.uflux)
            return [[xl,xh], [yc, yc],[bc,bc],[yl,yh]]
        else:
            x,y = [bc, cl.ul95/exp] if useexp else [bc,cl.frac]#[bc,cl.ul95/(cl.Naj[-1]+cl.ul95)]
            arrow = [y, y*0.8, y*0.8, y*0.6, y*0.8, y*0.8] if useexp else [y, y-0.025, y-0.025, y-0.05, y-0.025, y-0.025]
            arrow2 = [x, x,     x*1.1, x,     x/1.1, x] if useexp else [x, x,     x*1.03, x,     x/1.03, x]
            print ('95\% flux',y)
            return [[xl,xh], [y,y],arrow2,arrow,[xl,xh], [y,y],arrow2,arrow]

    fact = 1.0
    lss = ['-','--']
    fignames = [r'$\rm{[a]}$',r'$\rm{[b]}$',r'$\rm{[c]}$']

    for iteration,wid in enumerate(width):
        pts = []
        names = []
        py.clf()
        ofile = open(figdir+'%s_%1.1f_summary.txt'%(iname,wid),'w')
        for it,mod in enumerate(model):

            pickles=[]
            for emin,emax,pren in zip(energy2[:-1],energy2[1:],prepat):
                template = '/phys/groups/tev/scratch1/users/Fermi/mar0/figures/emi*%d_*%d_ec%d*bins16*%1.1f*%s__C%s_%s.pickle'%(emin,emax,ctype,wid,pren+iname,mod,irf)
                #print (template)
                pickle = np.sort(glob.glob(template))
                for pick in pickle:
                    pickles.append(pick)
            print (pickles)
            mx=0
            mi=1e40

            for it2,pickle in enumerate(pickles):
                
                if True:
                #try:
                    print ('Reading %s'%(pickle))
                    pfile = open(pickle)
                    cl = cPickle.load(pfile)
                    if type(cl.frac)==type([]):
                        cl.frac=cl.frac[)0]
                    print ('TS: ',cl.TS)
                    print ('Pval: ',cl.pval
                    texp = sum([exp.value(sds,cl.ebar) for sds in sdlist])/(cl.ebar*ergs) if uexp else 1
                    lines = marker(cl,texp,uexp,useupper)
                    tmx = max(lines[3])
                    tmi = min(lines[3])
                    mx = tmx if tmx>mx else mx
                    mi = tmi if tmi<mi and tmi>0 else mi
                    if it2==0:
                        p1 = py.plot(fact*np.array(lines[0]),lines[1],linewidth=1,color='k',ls=lss[it],label=r'$\rm{%s} $'%(model2[it]))
                    else:
                        p1 = py.plot(fact*np.array(lines[0]),lines[1],linewidth=1,color='k',ls=lss[it])
                    py.plot(fact*np.array(lines[2]),lines[3],linewidth=1,color='k',ls=lss[it])
                    #py.plot(fact*np.array(lines[4]),lines[5],'k--',linewidth=1,)
                    #py.plot(fact*np.array(lines[6]),lines[7],'k--',linewidth=1,color='k')
                    if uexp:
                        exstr = '%1.2e}'%(cl.ul95/texp)
                        exstr.replace('e','$x$^{')
                        print ('$%d-%d$ & $%1.1f$ & $%s$'%(cl.emin,cl.emax,max(0,cl.TS),exstr), file=ofile)
                    else:
                        print ('$%d-%d$ & $%1.1f$ & $%1.3f$'%(cl.emin,cl.emax,max(cl.TS,0),cl.frac), file=ofile)
                else:
                #except:
                    print ('Plotting failed')

            if len(pickles)>0:
                pts.append(p1)
                names.append(r'$\rm{%s} $'%(model2[it]))
        fact = fact/1.0
        #py.grid()
        if uexp:
            py.loglog()
            ypos = 0.75
            specs,lowers,uppers = [],[],[]
            for sd in sdlist:
                diff,ebrs,spec,lower,upper=getspec(sd)
                #for it4 in range(len(spec)):
                #    print (ebrs[it4],lower[it4],spec[it4],upper[it4])
                if specs==[]:
                    #texp = np.array([sum([exp.value(sds,ebr) for sds in sdlist])/(ebr*ergs) for ebr in ebrs])
                    specs=spec#/texp
                    lowers=lower
                    uppers=upper
                else:
                    specs+=spec#/texp
                    lowers+=lower
                    uppers+=upper
            if False:
                texp = np.array([sum([exp.value(sds,ebr) for sds in sdlist])/(ebr*ergs) for ebr in ebrs])
                llimit = 1.0/texp
                pa = py.plot(ebrs[specs>llimit],specs[specs>llimit],'r-',linewidth=2)
                pb = py.plot(ebrs[lowers>llimit],lowers[lowers>llimit],'r-',linewidth=1)
                pc = py.plot(ebrs[uppers>llimit],uppers[uppers>llimit],'r-',linewidth=1)
                pts.append(pa)
                names.append(r'$ Fermi\/\rm{Spectrum} $')
        else:
            py.semilogx()
            ypos = 0.75

            #plots specific to the HESS sources, with observations and models
        """if iname=='1es0229p200' or iname=='1es0347-121':
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/%s_spectrum.txt'%iname)
            data = 10**data.transpose()
            p5 = py.plot(data[0]*1000,data[1],'k',linewidth=2)
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/%s_tev.txt'%iname)
            data = 10**data.transpose()
            p6 = py.plot(data[0]*1000,data[1],'kd',linewidth=2,mfc='w')"""

        def expcutoff(en,norm,eref,pindex,ecut,sup=1.):
            return norm*((en/eref)**pindex)*np.exp(-(en/ecut)**sup)
        def uncpl(en,n0,e0,g0,sn0,sg0):
            return np.sqrt((sn0/n0)**2 + (np.log(en/e0)*sg0)**2)
        if (iname=='1es0229p200' or iname=='1es0347-121') and False:
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/%s_spectrum.txt'%iname)
            data = 10**data.transpose()
            p5 = py.plot(data[0]*1000,data[1],'k',linewidth=2,label='Total Model')
            names.append('Total Model')
            pts.append(p5)
        if iname=='1es0229p200':
            #py.cla()
            #py.loglog()
            energies = np.logspace(2,3,60)
            energies2 = np.logspace(3,5,60)
            energies1 = np.logspace(2,7,60)
            plf,plg,ple = 1.338e-13,-2.72,1e3
            splf,splg = 1.0e-13,0.5
            plg2 = -1.76
            splg2 = 0.35
            spectrum = [expcutoff(energies,plf*ergs*energies**2,ple,plg,1e14,2),expcutoff(energies2,plf*ergs*energies2**2,ple,plg2,1e14,2)]
            specerr = [uncpl(energies,plf,ple,plg,splf,splg)*spectrum[0],uncpl(energies2,plf,ple,plg2,splf,splg2)*spectrum[1]]
            energies = np.logspace(2,5,120)
            # py.plot(energies,np.hstack(spectrum),'k:',label=r'$\rm{gtlike\/model}$')
            # p11 = py.plot(energies,np.hstack(spectrum)+np.hstack(specerr),'k-',linewidth=2,label=r'$\rm{gtlike\/error}$')
            #py.plot(energies,np.hstack(spectrum)-np.hstack(specerr),'k-',linewidth=2)
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/%s_HESS.txt'%iname).transpose()
            p6 = py.errorbar(data[0]*1000,data[3],xerr=[data[1]*1000,data[2]*1000],yerr=[data[4],data[5]],ls='None',marker='d',mfc='w',color='k',label=r'$ \rm{H.E.S.S.} $')
            names.append('H.E.S.S.')
            pts.append(p6[0])
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/%s_VERITAS.txt'%iname).transpose()
            p7 = py.errorbar(data[0]*1000,data[5],xerr=[data[1]*1000,data[2]*1000],yerr=[data[6],data[7]],ls='None',marker='d',mfc='k',color='k',label=r'$ \rm{VERITAS} $')
            pts.append(p7[0])
            names.append('VERITAS')
            """data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/1es0229p200_model.txt')
            data = 10**data
            data = data.transpose()
            p8 = py.plot(data[0],data[1],'k-.')
            energies = np.logspace(2,7,30)
            spectrum = expcutoff(energies,10**(-11.09),1e6,0.8,1e14)
            print (expcutoff(1.e6,10**(-11.09),1.e6,0.8,1.e100))
            py.plot(energies,spectrum,'k-.',linewidth=1,label=r'$\rm{Deabsorbed}$')"""
        if iname=='1es0347-121':
            #py.cla()
            #py.loglog()
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/1ES0347-121_tev.txt')
            data = 10**data
            npts = len(data)/3
            pts1,upp,low=data[:npts].transpose(),data[npts:npts*2].transpose(),data[2*npts:].transpose()
            names.append('H.E.S.S.')
            fct = 10**(1./12.)
            p7 = py.errorbar(pts1[0]/1e6,pts1[1],xerr=[(1-1/fct)*pts1[0]/1e6,(fct-1)*pts1[0]/1e6],yerr=[pts1[1]-low[1],upp[1]-pts1[1]],ls='None',marker='d',mfc='w',color='k',label=r'$ \rm{H.E.S.S.}$ ')
            pts.append(p7[0])
            """data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/1ES0347-121_model.txt')
            data = 10**data
            data = data.transpose()
            p8 = py.plot(data[0]/1.e6,data[1],'k-.')
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/1ES0347-121_model2.txt')
            data = 10**data
            data = data.transpose()
            p9 = py.plot(data[0]/1.e6,data[1],'k-',linewidth=2)"""
            """
            data = np.loadtxt('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/1ES0347-121_model3.txt')
            data = 10**data
            data = data.transpose()
            p10 = py.plot(data[0],data[1]*ergs*1000,'k-.',linewidth=2,label=r'$\rm{10^{-15}\/G}$')"""
            energies = np.logspace(2,3,60)
            energies2 = np.logspace(3,5,60)
            energies1 = np.logspace(2,7,60)
            #plf,plg,ple = 2.35e-13,-1.66784,1e3
            #splf,splg = 9.42165e-14,0.195687
            plf,plg,ple = 2.159e-13,-2.07,1e3
            splf,splg = 8.8e-14,0.6
            plg2 = -1.62
            splg2 = 0.20
            spectrum = expcutoff(energies1,10**(-11.0686),1e6,0.45,1e14,2)
           # py.plot(energies1,spectrum,'k-.',linewidth=1,label=r'$\rm{Deabsorbed}$')
            spectrum = [expcutoff(energies,plf*ergs*energies**2,ple,plg,1e14,2),expcutoff(energies2,plf*ergs*energies2**2,ple,plg2,1e14,2)]
            specerr = [uncpl(energies,plf,ple,plg,splf,splg)*spectrum[0],uncpl(energies2,plf,ple,plg2,splf,splg2)*spectrum[1]]
            energies = np.logspace(2,5,120)
            #py.plot(energies,np.hstack(spectrum),'k:',label=r'$\rm{gtlike\/model}$')
            #p11 = py.plot(energies,np.hstack(spectrum)+np.hstack(specerr),'k-',linewidth=2,label=r'$\rm{gtlike\/error}$')
            #py.plot(energies,np.hstack(spectrum)-np.hstack(specerr),'k-',linewidth=2)
            #pts.append(p11)
            names.append('Finke et al. 2010')

        prop = mpl.font_manager.FontProperties(size=14)
        py.legend(prop=prop,numpoints=1)#pts,names,prop=prop)
        if uexp:
            py.ylabel(r'$\nu F_{\nu}\rm{\/(erg\/s^{-1}cm^{-2})}\/$')
        else:
            py.ylabel(r'$\rm{Fraction}$')
        py.xlabel(r'$\rm{Energy\/(MeV)}$')
        if uexp:
            py.ylim(1e-14,1e-9)
            py.xlim(100,1e7)
        else:
            py.ylim(1e-3,1)
            py.xlim(1000,1e5)
        py.figtext(0.25,ypos,fignames[iteration])
        axes = py.gca()
        axes.set_aspect('equal')
        py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/%s%s_%s_summary.eps'%(name,('%1.1f'%wid).replace('.',''),irf))#,mod
        py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/%s%s_%s_summary.png'%(name,('%1.1f'%wid).replace('.',''),irf))#,mod
        ofile.close()
    py.close(9)

## find the nearest catalog source and use its spectra, WARNING do not use if it's not in the catalog
#  @param sd SkyDir corresponding to source location
def getspec(sd):
    LOGP=True
    SRCLIB = pf.open('/phys/groups/tev/scratch1/users/Fermi/catalog/P72Y_uw26.fits')
    SRCTB = SRCLIB[1].data
    diff=INFINITY
    for it,src in enumerate(SRCTB):
        ra,dec=src.field('RA'),src.field('DEC')
        sd2 = s.SkyDir(float(ra),float(dec))
        tdiff=sd.difference(sd2)
        if tdiff<diff:
            diff=tdiff
            index=it
    myflux = SRCTB.field('Flux_Density')[index]; myindex = SRCTB.field('Spectral_Index')[index]; mycutoff = SRCTB.field('Cutoff_Energy')[index]
    myuflux = SRCTB.field('Unc_Flux_Density')[index]; myuindex = SRCTB.field('Unc_Spectral_Index')[index];myucutoff = SRCTB.field('Cutoff_Energy_Unc')[index]
    mypivot = SRCTB.field('Pivot_Energy')[index]
    myra = SRCTB.field('RA')[index]; mydec = SRCTB.field('DEC')[index]
    ebins = np.logspace(1,5.5,40)
    ebars = np.sqrt(ebins[:-1]*ebins[1:])
    print (myra,mydec,diff*57.,myflux,myindex,mypivot)

    bfun = lambda r: r**-myindex * np.sqrt(myuflux**2 + (myflux*np.log(r))**2 * (myuindex**2))

    if LOGP:
        mybeta = SRCTB.field('beta')[index]
        myubeta = SRCTB.field('Unc_beta')[index]
    else:
        mybeta = 0
        myubeta = 0
    if not (np.isnan(mybeta) or np.isinf(mybeta)) and LOGP and mybeta!=0:
        ec = LogParabola(e0=mypivot,p=[myflux,myindex,mybeta,mypivot])
        print ('Used Logparabola')
        upper = (ebars**2)*ergs*(ec(ebars)+bfun(ebars/mypivot))
        lower = (ebars**2)*ergs*(ec(ebars)/(1+bfun(ebars/mypivot)/ec(ebars)))
    else:
        if not  (np.isnan(mycutoff) or np.isinf(mycutoff)):
            ec = ExpCutoff(e0=mypivot,p=[myflux,myindex,mycutoff])
            print ('Used Expcutoff')
        else:
            print ('Used PowerLaw')
            ec = PowerLaw(e0=mypivot,p=[myflux,myindex])
        upper = (ebars**2)*ergs*(ec(ebars)+bfun(ebars/mypivot))
        lower = (ebars**2)*ergs*(ec(ebars)/(1+bfun(ebars/mypivot)/ec(ebars)))
    iflux = np.array([(ebars[it]**2)*ergs*ec(ebars[it]) for it in range(len(ebars))])#ec.i_flux(ebins[it],ebins[it+1],e_weight=1,cgs=True)
    return diff*rd,ebars,iflux,lower,upper

def sourcelistinfo(lis):
    ff = pf.open('/phys/groups/tev/scratch1/users/Fermi/catalog/gll_psc_v03.fit')
    iof = open('/phys/users/mar0/sourcelists/'+lis+'.txt')
    of = open('/phys/users/mar0/sourcelists/'+lis+'_tb.txt','w')
    tb = ff[1].data
    header = iof.readline()
    tbflux = np.array([float(src.field('Flux1000_3000')) for src in tb])
    for lines in iof:
        line = lines.split()
        name = line[0]
        sd = s.SkyDir(float(line[1]),float(line[2]))
        for src in tb:
            tsd = s.SkyDir(float(src.field('RA')),float(src.field('DEC')))
            diff = tsd.difference(sd)*180/np.pi
            if diff<0.25:
                signif = float(src.field('Sqrt_TS10000_100000')) #np.sqrt(float(src.field('Sqrt_TS1000_3000'))**2+float(src.field('Sqrt_TS3000_10000'))**2+**2)
                fluxes = float(src.field('Flux1000_3000'))+float(src.field('Flux3000_10000'))+float(src.field('Flux10000_100000'))
                string = '%s & %1.4f & %1.4f & %1.4f & %1.4f & %1.2f & %1.1f & %1.0f & %s\\\\'%(src.field('Source_Name'),float(src.field('RA')),float(src.field('DEC')),float(src.field('GLON')),float(src.field('GLAT')),float(src.field('Spectral_Index')),fluxes*1e9,signif,src.field('CLASS1'))
                print (string.replace('-','$-$'), file=of)
            if diff>0.25 and diff<4.:
                fluxes = float(src.field('Flux1000_3000'))
                rank = len(tbflux[tbflux>fluxes])
                if rank<300:
                    print ('Found a source %s within 2 degrees of %s with a flux of %1.2e (%d)'%(src.field('Source_Name'),name,fluxes,rank))

################## plots for the paper ###################
############## pulsar and agn bins     ###################
def tabulatepulsars():
    maxr,days,ctype,sel,irf,bins=4,730,0,'[0:2]','P7SOURCE_V6',12
    ebins = [1000,3162,10000,31623,100000]
    sels = ['[0:2]','[0:2]','[0:2]','[0:0]']
    agnlists = [[],[],[],['agn_redshift2_lo_bzb']]
    cls = []
    tbfile = open('pulsartable.tex','w')
    ffile = open('pulsartable.txt','w')
    it = 0
    for emin,emax in zip(ebins[:-1],ebins[1:]):

        header = r'\tablehead{Bin edges (deg) & $m_i$ '
        cl = CombinedLike(irf=irf,mode=-1,pulsars = eval('pulsars'+sels[it]+''),agnlist=agnlists[it],qandd=True)
        print (agnlists[it])
        print (eval('pulsars'+sels[it]+''))
        cl.loadphotons(0,maxr,emin,emax,239557417,239517417+days*86400,ctype)
        cl.bindata(bins)
        cl.fit(qandd=True)
        psr =  len(cl.pulsars)>0
        print (r'\begin{deluxetable}{%s}'%('r'*(2+5*len(cl.pulsars)+4*len(cl.agns))), file=tbfile)
        print (r'\tabletypesize{\scriptsize}', file=tbfile)
        print (r'\tablecaption{Statistics for %s in the energy range $%d-%d$ for the analysis in Section 4.%d}'%('the Vela and Geminga pulsars' if psr else 'the low-redshift BL Lacs',emin,emax,2 if psr else 3), file=tbfile)
        print (r'\tablewidth{0pt}', file=tbfile)
        for puls in cl.pulsars:
            header += r'& $\nu_i^{off}$ & $n_i^{off}$ & $N_{psr}m_i$ & $\nu_i^{on}$ & $n_i^{on}$'
        for agn in cl.agns:
            header += r'& $N^{agn}m_i$ & $N^{iso}b_i$ & $\nu_i^{agn}$ & $n_i^{agn}$'
        header+=r'}'
        print (header, file=tbfile)
        print (r'\startdata', file=tbfile)
        for it2 in range(cl.nbins):
            line = r'$%1.3f-%1.3f$ & %1.3f '%(cl.angbins[it2],cl.angbins[it2+1],cl.psf[it2])
            for it3,puls in enumerate(cl.pulsars):
                line += r'& $%1.1f$ & %d & $%1.1f$ & $%1.1f$ & %d'%(cl.vij[it3][it2],cl.pofhists[it3][it2],cl.Npj[it3]*cl.psf[it2],cl.Npj[it3]*cl.psf[it2]+puls[1]*cl.vij[it3][it2],cl.ponhists[it3][it2])
            for it3,agn in enumerate(cl.agns):
                line += r'& $%1.1f$ & $%1.1f$ & $%1.1f$ & %d'%(cl.Naj[it3]*cl.psf[it2],cl.Ni[it3]*cl.iso[it2],cl.Naj[it3]*cl.psf[it2]+cl.Ni[it3]*cl.iso[it2],cl.agnhists[it3][it2])
            line += r'\\'
            print (line, file=tbfile)
        print (r'\enddata', file=tbfile)
        print (r'\end{deluxetable}', file=tbfile)
        print (cl)
        cls.append(cl)
        it+=1
    tbfile.close()
    ffile.close()
    return cls

############### light curves for pulsars  ##################
def makephaseplots(bins=100):
    import uw.utilities.fitstools as uuf

    vtb = pf.open(pulsdir+'vela-ft1.fits')[1].data
    vtb = vtb[vtb.field('THETA')<66.4]
    ras,decs,phases = vtb.field('RA'),vtb.field('DEC'),vtb.field('PULSE_PHASE')
    mask = uuf.rad_mask(ras,decs,s.SkyDir(128.8391,-45.1792),4.,True)
    vtb = vtb[mask]
    pmask = np.sort(phases[mask])
    steps = len(pmask)/bins
    pbins = [0 if frac ==0 else pmask[(frac)*steps] for frac in range(bins)]
    pbins.append(1.)
    pbins=np.array(pbins)
    parea = np.array([pbins[it+1]-pbins[it] for it in range(len(pbins)-1)])
    phist = np.histogram(pmask,bins=pbins)
    py.figure(figsize=(8,4))
    widths = parea/2
    mids = np.array([(pbins[it+1]+pbins[it])*0.5 for it in range(len(pbins)-1)])
    for it,energy in enumerate([100]):#,316,1000,3162]):
        #hist = np.histogram(phases[vtb.field('ENERGY')>energy],bins=bins)
        #mid = (hist[1][:-1]+hist[1][1:])*0.5
        xhists = pbins#,pbins[:-1]+1,2])
        yhists = np.hstack([phist[0]/parea,phist[0][0]/parea[0]])/1e6#,phist[0]/parea,phist[0][0]/parea[0]])/1e6
        print (phist[0],len(phist[0]))
        onmask = ((pbins>0.1) & (pbins<0.15))# | ((pbins>0.5) & (pbins<0.6))
        p1 = py.step(np.hstack([0.105,xhists[onmask],0.15]),np.hstack([0,yhists[onmask],0]),color='grey',where='post',fillstyle='full',linewidth=2)
        onmask = (pbins>0.7)
        p2 = py.step(np.hstack([0.707,xhists[onmask],1.0]),np.hstack([0,yhists[onmask],0]),color='black',where='post',fillstyle='full',linewidth=2)
        py.legend((r'$\rm{On\/pulse}$',r'$\rm{Off\/pulse}$'))
        py.step(xhists,yhists,where='post',color='black')
        onmask = ((pbins>0.1) & (pbins<0.15))# | ((pbins>0.5) & (pbins<0.6))
        py.step(np.hstack([0.105,xhists[onmask],0.15]),np.hstack([0,yhists[onmask],0]),color='grey',where='post',fillstyle='full',linewidth=2)
        onmask = (pbins>0.7)
        py.step(np.hstack([0.707,xhists[onmask],1.0]),np.hstack([0,yhists[onmask],0]),color='black',where='post',fillstyle='full',linewidth=2)
        onmask = ((pbins>0.5) & (pbins<0.6))
        py.step(np.hstack([0.5,xhists[onmask],0.597]),np.hstack([0,yhists[onmask],0]),color='grey',where='post',fillstyle='full',linewidth=2)
        #py.step(np.hstack([mid,mid+1]),np.hstack([hist[0]-min(hist[0])+1,hist[0]-min(hist[0])+1]),where='mid',color='%1.2f'%(1.-1./(it+1)))
        #py.step(,hist[0]-min(hist[0])+1,where='mid',color='%1.2f'%(1.-1./(it+1)))
    #py.semilogy()
    py.xlabel(r'$\rm{Pulsar\/Phase}$')
    py.ylabel(r'$\rm{10^{6}\/Events/Bin\/Width}$')
    py.ylim(0.01*max(yhists),1.2*max(yhists))
    py.savefig(figdir+'velaphase.png')
    py.savefig(figdir+'velaphase.eps')
    py.clf()
    vtb = pf.open(pulsdir+'gem-ft1.fits')[1].data
    vtb = vtb[vtb.field('THETA')<66.4]
    ras,decs,phases = vtb.field('RA'),vtb.field('DEC'),vtb.field('PULSE_PHASE')
    mask = uuf.rad_mask(ras,decs,s.SkyDir(98.47855,17.7739),4.,True)
    vtb = vtb[mask]
    pmask = np.sort(phases[mask])
    steps = len(pmask)/bins
    pbins = [0 if frac ==0 else pmask[(frac)*steps] for frac in range(bins)]
    pbins.append(1.)
    pbins=np.array(pbins)
    parea = np.array([pbins[it+1]-pbins[it] for it in range(len(pbins)-1)])
    phist = np.histogram(pmask,bins=pbins)
    #py.figure(figsize=(16,8))
    widths = parea/2
    mids = np.array([(pbins[it+1]+pbins[it])*0.5 for it in range(len(pbins)-1)])
    for it,energy in enumerate([100]):#,316,1000,3162]):
        #hist = np.histogram(phases[vtb.field('ENERGY')>energy],bins=bins)
        #mid = (hist[1][:-1]+hist[1][1:])*0.5
        xhists = pbins#,pbins[:-1]+1,2])
        yhists = np.hstack([phist[0]/parea,phist[0][0]/parea[0]])/1e6#,phist[0]/parea,phist[0][0]/parea[0]])/1e6
        print (phist[0],len(phist[0]))
        onmask = ((pbins>0.6) & (pbins<0.68))
        py.step(np.hstack([0.6,xhists[onmask],0.68]),np.hstack([0,yhists[onmask],0]),color='grey',where='post',fillstyle='full',linewidth=2)
        onmask = (pbins>0.25) & (pbins<0.55)
        py.step(np.hstack([0.25,xhists[onmask],0.547]),np.hstack([0,yhists[onmask],0]),color='black',where='post',fillstyle='full',linewidth=2)
        py.legend((r'$\rm{On\/pulse}$',r'$\rm{Off\/pulse}$'))
        py.step(xhists,yhists,where='post',color='black')
        onmask = ((pbins>0.1) & (pbins<0.17))# | ((pbins>0.5) & (pbins<0.6))
        py.step(np.hstack([0.105,xhists[onmask],0.17]),np.hstack([0,yhists[onmask],0]),color='grey',where='post',fillstyle='full',linewidth=2)
        onmask = ((pbins>0.6) & (pbins<0.68))
        py.step(np.hstack([0.6,xhists[onmask],0.68]),np.hstack([0,yhists[onmask],0]),color='grey',where='post',fillstyle='full',linewidth=2)
        onmask = (pbins>0.25) & (pbins<0.55)
        py.step(np.hstack([0.25,xhists[onmask],0.547]),np.hstack([0,yhists[onmask],0]),color='black',where='post',fillstyle='full',linewidth=2)
        #py.step(np.hstack([mid,mid+1]),np.hstack([hist[0]-min(hist[0])+1,hist[0]-min(hist[0])+1]),where='mid',color='%1.2f'%(1.-1./(it+1)))
        #py.step(,hist[0]-min(hist[0])+1,where='mid',color='%1.2f'%(1.-1./(it+1)))
    py.xlabel(r'$\rm{Pulsar\/Phase}$')
    py.ylabel(r'$\rm{10^{6}\/Events/Bin\/Width}$')
    py.ylim(0.01*max(yhists),1.2*max(yhists))
    py.savefig(figdir+'gemphase.png')
    py.savefig(figdir+'gemphase.eps')