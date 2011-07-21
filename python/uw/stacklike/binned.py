import skymaps as s
import pylab as py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from uw.stacklike.stacklike import *
from uw.stacklike.angularmodels import *
from uw.utilities.minuit import Minuit
import uw.utilities.assigntasks as ua
from uw.like.pycaldb import CALDBManager
from uw.like.pypsf import CALDBPsf
from uw.like.quadform import QuadForm,Ellipse
from uw.stacklike.limits import BayesianLimit
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

def format_error(v,err):
    if v>0 and err>0:
        logz = math.floor(math.log10(err))-1
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

pulsdir = '/phys/groups/tev/scratch1/users/Fermi/mar0/data/pulsar/'               #directory with pulsar ft1 data
agndir = '/phys/groups/tev/scratch1/users/Fermi/mar0/data/6.3/'                   #directory with agn ft1 data
srcdir = '/phys/users/mar0/sourcelists/'                                          #directory with lists of source locations (name ra dec)
cachedir = '/phys/users/mar0/cache/'
pulsars = [('vela',1.0),('gem',1.0)]                                                          #pulsar sourcelist name 
agnlist = ['agn-psf-study-bright']
irf = 'P6_v8_diff'
rd = 180./np.pi

###############################################  CombinedLike Class ############################################

## Implements the likelihood analysis outlined by Matthew Wood
## available here: http://fermi-ts.phys.washington.edu/pointlike/common/psf_lnl.pdf
class CombinedLike(object):


    ## constructor for Combinedlikelihood class
    #  @param pulsdir pulsar ft1 data directory
    #  @param agndir agn ft1 data directory
    #  @param srcdir sourcelist directory (see Stacklike documentation)
    #  @param pulsars (sourcelist,alpha) names of pulsar list and ratio of on/off pulse
    #  @param agnlist list of agn sources
    #  @param irf intial response function
    def __init__(self,**kwargs):
        self.pulsdir = pulsdir
        self.agndir = agndir
        self.srcdir = srcdir
        self.pulsars = pulsars
        self.agnlist = agnlist
        self.cachedir = cachedir
        self.irf = irf
        self.cache = True
        self.verbose = False
        self.pulse_ons = []         #pulsar on pulse angular distributions
        self.pulse_offs = []        #pulsar off pulse angular distributions
        self.agns = []              #agn angular distributions
        self.angbins=[]             #angular bins
        self.midpts = []            #midpoints of angular bins
        self.ponhists=[]            #pulsar on pulse histograms
        self.pofhists=[]            #pulsar off pulse histograms
        self.agnhists=[]            #agn histograms
        self.nbins = 0              #number of angular bins
        self.halomodel=''
        self.haloparams=[]
        self.ctmin=0.3
        self.ctmax=1.0
        self.mode=-1
        self.__dict__.update(kwargs)
        self.TS = 0

    def __str__(self):
        verbose = '\n'
        verbose = verbose + '--------- PSF    -------------\n'
        verbose = verbose + 'l-edge(deg)  r-edge(deg)    fraction   fraction error\n'
        for it, x in enumerate(self.psf):
            tstr = '%8.3g   %8.3g'%(self.angbins[it],self.angbins[it+1]), format_error(self.psf[it],self.psfe[it])
            verbose = verbose + string.join(tstr) + '\n'
        verbose = verbose +'\n'
        verbose = verbose + '--------- Pulsars-------------\n'
        for it,nj in enumerate(self.Npj):
            verbose =verbose + 'N(%s) = %1.0f (1 +/- %1.3f)\n'%(self.pulsars[it][0],nj,self.Npje[it]/nj)
        verbose = verbose + '\n'
        verbose = verbose + '--------- AGN    -------------\n'
        for it,nj in enumerate(self.Naj):
            verbose = verbose + 'Npsf(%s) = %1.0f (1 +/- %1.3f)\n'%(self.agnlist[it],nj,self.Naje[it]/nj)
            verbose = verbose + 'Niso(%s) = %1.0f (1 +/- %1.3f)\n'%(self.agnlist[it],self.Ni[it],self.Nie[it]/self.Ni[it])
        verbose = verbose + '\n'
        if self.halomodel!='':
            verbose = verbose + 'Nhalo = %1.0f (1 +/- %1.3f)\n'%(self.Nh,self.Nhe/self.Nh)
            verbose = verbose + 'Size(deg) = %1.3f\n'%(self.haloparams[0]*rd)
            verbose = verbose + 'TS = -2(LogL(Nhalo) - LogL(0)) = %1.2f\n'%(self.TS)
            tspval = self.tspval()
            verbose = verbose + 'Significance = %1.1f    p-val = %1.4f'%(tspval[1],tspval[0])
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

        tag = '%1.2f%1.2f%1.2f%1.2f%1.0f%1.0f%1.0f%1.2f%1.2f'%(minroi,maxroi,emin,emax,tmin,tmax,ctype,self.ctmin,self.ctmax)

        #load pulsar data
        thetabar = 0.
        photons = 0
        for psr in self.pulsars:

            if os.path.exists(self.cachedir+'%son%s.npy'%(psr[0],tag)):
                if self.verbose:
                    print 'Loaded %s on pulse from cache'%(psr[0])
                hist = np.load(self.cachedir+'%son%s.npy'%(psr[0],tag))
                self.pulse_ons.append(hist)
            else:
                sl = StackLoader(lis=psr[0],irf=self.irf,srcdir=self.srcdir,useft2s=False,ctmin=self.ctmin,ctmax=self.ctmax)
                sl.files = [self.pulsdir + psr[0]+ 'on-ft1.fits']
                sl.loadphotons(minroi,maxroi,emin,emax,tmin,tmax,ctype)
                sl.getds()
                photons = photons+len(sl.photons)
                thetabar = thetabar+sum([p.ct for p in sl.photons])
                self.pulse_ons.append(np.array(cp.copy(sl.ds)))
                if self.cache:
                    np.save(self.cachedir+'%son%s.npy'%(psr[0],tag),np.array(cp.copy(sl.ds)))
                del sl

            if os.path.exists(self.cachedir+'%soff%s.npy'%(psr[0],tag)):
                if self.verbose:
                    print 'Loaded %s off pulse from cache'%(psr[0])
                hist = np.load(self.cachedir+'%soff%s.npy'%(psr[0],tag))
                self.pulse_offs.append(hist)
            else:
                sl = StackLoader(lis=psr[0],irf=self.irf,srcdir=self.srcdir,useft2s=False,ctmin=self.ctmin,ctmax=self.ctmax)
                sl.files = [self.pulsdir + psr[0]+ 'off-ft1.fits']
                sl.loadphotons(minroi,maxroi,emin,emax,tmin,tmax,ctype)
                photons = photons+len(sl.photons)
                thetabar = thetabar+sum([p.ct for p in sl.photons])
                sl.getds()
                self.pulse_offs.append(np.array(cp.copy(sl.ds)))
                if self.cache:
                    np.save(self.cachedir+'%soff%s.npy'%(psr[0],tag),np.array(cp.copy(sl.ds)))
                del sl

        #load agn data
        for lists in self.agnlist:
            if os.path.exists(self.cachedir+'%s%s.npy'%(lists,tag)):
                if self.verbose:
                    print 'Loaded %s from cache'%(lists)
                hist = np.load(self.cachedir+'%s%s.npy'%(lists,tag))
                self.agns.append(hist)
            else:
                sl = StackLoader(lis=lists,irf=self.irf,srcdir=self.srcdir,useft2s=False,ctmin=self.ctmin,ctmax=self.ctmax)
                sl.files = glob.glob(self.agndir+'*.fits')
                sl.loadphotons(minroi,maxroi,emin,emax,tmin,tmax,ctype)
                photons = photons+len(sl.photons)
                thetabar = thetabar+sum([p.ct for p in sl.photons])
                sl.getds()
                self.agns.append(np.array(cp.copy(sl.ds)))
                if self.cache:
                    np.save(self.cachedir+'%s%s.npy'%(lists,tag),np.array(cp.copy(sl.ds)))
                del sl
        self.thetabar = 0.5*(self.ctmin+self.ctmax) if photons==0 else thetabar/photons
        self.photons=photons
            

    ######################################################################
    #   Bins the FT1 data into angular bins that are adaptive or fixed   #
    ######################################################################
    ## adaptively bins angular distributions in angle or sqrt
    #  @param bins number of angular bins, -1 for sqrt(N) bins
    def bindata(self,bins=8):
        alldata = []
        
        #determine needed bins by combining all data (except background)
        chist=[]
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


        #adaptive binning
        if bins>0 and len(alldata)>20:
            chist = np.array([sum(chist[:x+1]) for x in range(len(chist))])
            chist = chist/max(chist)
            cumm = np.array([(1.*x+1.)/len(alldata) for x in range(len(alldata))])      #cumulative dist function
            ct = (1.*np.arange(0,bins+1,1))/bins                                          #cumulative fractions, [0,1/bins,2/bins...(bins-1)/bins]
            mask = np.array([max(0,len(chist[chist<x])-1) for x in ct])
            xbins = alldata[mask]#np.array([alldata[max(0,len(chist[chist<x])-1)] for x in ct])         #bin edges corresponding to fractions
            self.angbins = xbins

        # sqrt(N) binning
        else:
            if len(alldata)>20:
                bins = int(np.sqrt(len(alldata)))
            else:
                bins = min(32,int(np.sqrt(len(self.agns[0]))/1.5))
            xbins = (np.arange(0,bins,1)/(1.*bins))**2
            if len(alldata)>0:
                minimum = min(min(alldata),min(self.agns[0]))*rd
            else:
                minimum = min(self.agns[0])*rd
            xbins =  minimum + xbins*(self.maxroi-minimum)
            self.angbins = xbins/rd

        #did you already do it?
        if self.ponhists==[]:
            for it1,puls in enumerate(self.pulse_ons):
                self.ponhists.append(np.histogram(puls,self.angbins)[0])
                self.pofhists.append(np.histogram(self.pulse_offs[it1],self.angbins)[0])
            for it1,agns in enumerate(self.agns):
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
    def profile(self,ip,x,**kwargs):
        self.__dict__.update(kwargs)
        if np.isscalar(x):
            params = cp.copy(self.params)
            params[ip] = x
            fixed = cp.copy(self.fixed)
            fixed[ip] = True

            self.minuit = Minuit(self.likelihood,params,#gradient=self.gradient,force_gradient=1,
                                 fixed=fixed,limits=self.limits,strategy=2,tolerance=0.0001,printMode=self.mode)
            self.minuit.minimize()
            return self.minuit.fval
        else:
            params = cp.copy(self.params)
            vals = []
            for xval in x:
                params[ip] = xval
                fixed = cp.copy(self.fixed)
                fixed[ip] = True

                self.minuit = Minuit(self.likelihood,params,#gradient=self.gradient,force_gradient=1,
                                     fixed=fixed,limits=self.limits,strategy=2,tolerance=0.0001,printMode=self.mode)
                self.minuit.minimize()
                vals.append(self.minuit.fval)
            return np.array(vals)

    ######################################################################
    #    Sets up Minuit and determines maximum likelihood parameters     #
    ######################################################################
    ## maximizes likelihood
    # possible keyword arguments
    # @param halomodel string corresponding to any angular model in angularmodels.py
    # @param haloparams array of model parameters for halomodel
    def fit(self,**kwargs):
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
            print ''
            print '**********************************************************'
            print '*                                                        *'
            print '*             Combined Likelihood Analysis               *'
            print '*                                                        *'
            print '**********************************************************'
            print 'Pulsars: %s'%psrs
            print 'AGNs: %s'%agns
            print 'Using halo model: %s'%self.halomodel
            print 'Using parameters: %s'%pars
            print '**********************************************************'

        psf = CALDBPsf(CALDBManager(irf=self.irf))
        fint = psf.integral(self.ebar,self.ctype,max(self.angbins)/rd,min(self.angbins)/rd)
        self.psfm = np.array([psf.integral(self.ebar,self.ctype,self.angbins[it+1]/rd,self.angbins[it]/rd)/fint for it in range(self.nbins)])

        #set up initial psf parameters (uniform)
        self.params = [self.psfm[x] for x in range(self.nbins-1)]
        self.limits = [[-1,1] for x in range(self.nbins-1)]
        self.fixed = [False for x in range(self.nbins-1)]

        alims = sum([sum(agnhist) for agnhist in self.agnhists])

        psrs = len(self.ponhists)
        agns = len(self.agnhists)

        #pulsar estimators
        for hist in self.ponhists:
            self.params.append(sum(hist)/2.)
            self.limits.append([0,sum(hist)*100])
            self.fixed.append(False)

        #agn estimators
        for hist in self.agnhists:
            self.params.append(alims/2.)
            self.limits.append([0,alims*100])
            self.fixed.append(False)

        #iso estimator
        self.params.append(1.)
        self.limits.append([0,alims*100])
        self.fixed.append(len(self.agnhists)==0)

        self.Nh=0
        self.Nhe=1e-40
        if self.halomodel=='':
            self.hmd=np.zeros(self.nbins)
        else:
            halomodel = eval(self.halomodel)
            if self.haloparams[0]<0:
                return np.Infinity
            mod = halomodel(lims=[min(self.angbins)/rd,max(self.angbins)/rd],model_par=self.haloparams)
            mint = mod.integral(min(self.angbins)/rd,max(self.angbins)/rd)
            self.hmd = np.array([mod.integral(self.angbins[it]/rd,self.angbins[it+1]/rd)/mint for it in range(self.nbins)])
            self.hmd = self.hmd/sum(self.hmd)

        self.params.append(self.Nh)
        self.limits.append([0,alims*100])
        self.fixed.append(self.halomodel=='')
        if self.verbose:
            print 'Setting up Minuit and maximizing'
        ############  setup Minuit and optimize  ###############
        self.minuit = Minuit(self.likelihood,self.params,#gradient=self.gradient,force_gradient=1,
                             fixed=self.fixed,limits=self.limits,strategy=2,tolerance=0.0001,printMode=self.mode)
        self.minuit.minimize()
        if self.verbose:
            print 'Likelihood value: %1.1f'%self.minuit.fval
            print '**********************************************************'
        #self.gradient(self.minuit.params,True)
        ###########  get parameters and errors #################

        self.cov = self.minuit.errors()
        self.errs = np.sqrt(np.diag(self.cov))
        self.params = cp.copy(self.minuit.params)

        """if any(self.errs<0):
            self.minuit.minuit.mnmnos()
            self.errs = []
            for it in range(len(self.minuit.params)):
                eplus,eminus,ecurv,gcc = Double(),Double(),Double(),Double()
                self.minuit.minuit.mnerrs(it,eplus,eminus,ecurv,gcc)
                self.errs.append(float(ecurv))
            self.errs = np.array(self.errs)"""

        #for i in range(len(self.minuit.params)):
        #    print '%12.4g +/- %12.4g'%(self.minuit.params[i], self.errs[i])

        #for i in range(len(self.cov)):
        #    for j in range(len(self.cov[i])):
        #        print '%12.4g'%self.cov[i][j],
        #    print

        #if any(self.errs<0):
        self.errs2 = [self.errors(it) for it in range(len(self.minuit.params))]#self.errs[it][it] for it in range(len(self.minuit.params))]
        #self.errs2 = self.errs
        self.errs2 = [self.finderrs(it) for it in range(len(self.minuit.params))]
        self.errs=np.sqrt(self.errs2)

        #scale = sum(self.psf)                                                                                          #normalize psf
        nmu = self.nbins - 1        
        self.psf = self.params[:nmu]
        self.psfe = np.array([self.errs[i] for i in range(nmu)])                             #PSF errors
        # Compute the value and error of the nth PSF weight
        self.psf = np.append(self.psf,[1-np.sum(self.psf)])
        dxdmu = np.zeros(len(self.params))
        dxdmu[:nmu] = -1
        mun_err = np.sqrt(np.dot(np.dot(self.cov,dxdmu),dxdmu))        
        self.psfe = np.append(self.psfe,[mun_err])

        self.Npj  = self.params[nmu:nmu+psrs] #PSR number est
        self.Npje = self.errs[nmu:nmu+psrs]         #PSR number est errors
        self.Naj  = self.params[nmu+psrs:nmu+psrs+agns]   #AGN number est
        self.Naje = self.errs[nmu+psrs:nmu+psrs+agns]           #AGN number errors
        self.Ni   = self.params[nmu+psrs+agns:nmu+psrs+agns+agns] #isotropic number est
        self.Nie  = self.errs[nmu+psrs+agns:nmu+psrs+agns+agns]         #isotropic number est errors
        
        if self.halomodel!='':
            self.Nh  = self.params[nmu+psrs+agns+1]  #halo est
            self.Nhe = self.errs[nmu+psrs+agns+1]           #halo err                               #halo err


        #for it in range(len(self.params)):
        #    print np.sqrt(self.errs[it][it]),np.sqrt(self.errs2[it])

        """seps = np.arange(-2.0,2.1,0.1)
        py.figure(figsize=(16,16))
        rows = int(np.sqrt(len(self.params)))+1

        for it in range(len(self.params)):
            py.subplot(rows,rows,it+1)
            er = np.zeros(len(self.params))
            er[it] = np.sqrt(self.errs2[it])
            like = []
            for x in seps:
                modif = self.params+er*x
                tlike = self.minuit.fval-self.likelihood(modif)
                like.append(tlike)
            py.plot(seps,like)
        py.savefig('likes.png')"""

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
                cov = self.cov[it2][nmu+it1]                              #covariance of psf and number estimator
                Ne = self.Npje[it1]
                me = self.psfe[it2]
                ve = np.sqrt((dvdm*me)**2+(dvdN*Ne)**2)      #naive error propagation
                tvije.append(ve)

            self.vij.append(tvij)
            self.vije.append(tvije)

        self.vij = np.array(self.vij)                                        #PSR background number estimator
        self.vije = np.array(self.vije)                                      #PSR background number estimator errors

        return self.minuit.fval

    #######################################################################
    #                  Summary of Likelihood Analysis                     #
    #######################################################################
    ## printResults() - summary
    def printResults(self):
        print str(self)

    ######################################################################
    #      Likelihood function from (3) in paper                         #
    ######################################################################
    ## likelihood function
    #  @param params likelihood function parameters
    def likelihood(self,params,verb=False):#npij,bpij,naij,mi,Npj,Naj,Ni,Nh=0):

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
        Nh = params[nmu+psrs+agns+agns:nmu+psrs+agns+agns+1][0]
        acc = 0
        #set to true for slow,verbose output

        mi = np.append(mi,[1-np.sum(mi)])
        if verb:
            print '**************************************************************'
            print '--------------------------------------------------------------'
            print '                        Pulsars                               '
            print '--------------------------------------------------------------'
            print 'n\tb\tNi\tvi\tmod\tcont1\tcont2\tacc'

        ########################################
        #          first sum in (3)            #
        ########################################
        #loop over pulsars
        for it1,row in enumerate(npij):
            if verb:
                print '--------------------------------------------------------------'
                print '                        %s                               '%self.pulsars[it1][0]
                print '--------------------------------------------------------------'
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
                if lterm <0. or v<0.:
                    print lterm,v
                    return np.Infinity

                cont1 = -n*np.log(lterm) if n>0 else 0
                cont1 = cont1 + N*m + a*v
                acc = acc + cont1

                cont2 = -b*np.log(v) if b>0 else 0
                cont2 = cont2 + v
                acc = acc + cont2

                if verb:
                    print '%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f'%(n,b,N*m,a*v,lterm,cont1,cont2,acc)
                    t.sleep(0.25)
        
        if verb:
            print '--------------------------------------------------------------'
            print '                        AGN                                   '
            print '--------------------------------------------------------------'
            print 'n\tNi\tIi\tHi\tlterm\tcont1\tacc'
        ########################################
        #         second sum in (3)            #
        ########################################
        #loop over agn
        for it1,row in enumerate(naij):

            if verb:
                print '--------------------------------------------------------------'
                print '                        %s                               '%self.agnlist[it1]
                print '--------------------------------------------------------------'
            #loop over angular bins
            for it2,bin in enumerate(row):

                #make sure log term is proper
                lterm = Naj[it1]*mi[it2]+Ni[0]*self.iso[it2] + Nh*self.hmd[it2]
                if lterm<0.:
                    return np.Infinity

                cont1 = -bin*np.log(lterm) if bin>0 else 0
                cont1 = cont1 + Naj[it1]*mi[it2] + Ni[it1]*self.iso[it2] + Nh*self.hmd[it2]
                acc = acc + cont1

                if verb:
                    print '%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f\t%1.4f'%(bin,Naj[it1]*mi[it2],Ni[it1]*self.iso[it2],Nh*self.hmd[it2],lterm,cont1,acc)
                    t.sleep(0.25)
        return acc

    ######################################################################
    #      Find single or double PSF fits to fractions                   #
    ######################################################################
    ## finds sigma gamma and fraction
    # @param double fits a double PSF
    def fitpsf(self,double=False):
        psf = CALDBPsf(CALDBManager(irf='P6_v11_diff'))
        de = 0.45
        if double:
            if psf.newstyle:
                nc,nt,gc,gt,sc,st,w = psf.get_p(self.ebar,self.ctype)
                pars = [sc[0],gc[0],st[0],gt[0],nc[0]]
                lims = [[0,100],[1,100],[0,100],[1,100],[0.5-de,0.5+de]]
            else:
                gc,si,w = psf.get_p(self.ebar,self.ctype)
                pars = [si[0],gc[0]*1.2,si[0],gc[0]*0.8,0.5]
                lims = [[0,100],[1,100],[0,100],[1,100],[0.5-de,0.5+de]]
        else:
            if psf.newstyle:
                nc,nt,gc,gt,sc,st,w = psf.get_p(self.ebar,self.ctype)
                pars = [sc[0],gc[0]]
                lims = [[0,100],[1,100]]
            else:
                gc,si,w = psf.get_p(self.ebar,self.ctype)
                pars = [si[0],gc[0]]
                lims = [[0,100],[1,100]]
        tmin = Minuit(lambda x: self.psfchisq(x),pars,limits=lims,printMode=-1,up=1,strategy=1)
        tmin.minimize()
        print tmin.params
        print tmin.fval
        if not double:
            psf1 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=tmin.params)
            print psf1.rcl(0.68),psf1.rcl(0.95)
            print psf.inverse_integral(self.ebar,self.ctype,68.),psf.inverse_integral(self.ebar,self.ctype,95.)
        else:
            psf1 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=tmin.params[0:2])
            psf2 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=tmin.params[2:4])
            alph = tmin.params[4]
            #print psf1.rcl(0.68),psf1.rcl(0.95)
            print alph*psf1.rcl(0.68)+(1-alph)*psf2.rcl(0.68),alph*psf1.rcl(0.95)+(1-alph)*psf2.rcl(0.95)
            print psf.inverse_integral(self.ebar,self.ctype,68.),psf.inverse_integral(self.ebar,self.ctype,95.)
        return np.insert(tmin.params,0,self.ebar)


    ######################################################################
    #      Find single or double PSF fits to fractions                   #
    ######################################################################
    ## finds sigma gamma and fraction
    # @param double fits a double PSF
    def psfchisq(self,pars):
        psfs = []
        fints = []
        psf1 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=[pars[0],pars[1]])
        fint1 = psf1.integral(psf1.lims[0],psf1.lims[1])
        psfs.append(psf1)
        fints.append(fint1)
        alphas = [1.,0]
        if len(pars)>2:
            alphas = [pars[4],1-pars[4]]
            psf2 = PSF(lims=[min(self.angbins),max(self.angbins)],model_par=[pars[2],pars[3]])
            fint2 = psf1.integral(psf2.lims[0],psf2.lims[1])
            psfs.append(psf2)
            fints.append(fint2)
        tints = np.zeros(self.nbins)
        for it,psf in enumerate(psfs):
            cint = alphas[it]*np.array([psf.integral(self.angbins[it2],self.angbins[it2+1])/fints[it] for it2 in range(self.nbins)])
            tints = tints + cint
        chisqa = [((tints[it]-self.psf[it])/self.psfe[it])**2 if self.psfe[it]>0 else 0 for it in range(len(tints))]
        chisq = sum(chisqa)
        #print pars,chisqa,chisq
        #t.sleep(0.25)
        return chisq


    ######################################################################
    #      Makes plots of PSR, AGN fits, PSF and background residuals    #
    ######################################################################
    ## outputs a plot of the distributions
    # @param name output PNG filename
    def makeplot(self,name):
        py.ioff()
        py.figure(1,figsize=(32,16))
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
            p1 = py.errorbar(self.midpts,(hist)/self.areas,xerr=self.widths,yerr=np.sqrt(hist)/self.areas,marker='o',ls='None')
            p3 = py.errorbar(self.midpts,(self.Npj[it]*self.psf+self.pulsars[it][1]*self.vij[it])/self.areas,xerr=self.widths,yerr=np.sqrt((self.Npje[it]*self.psf)**2+(self.Npj[it]*self.psfe)**2)/self.areas,marker='o',ls='None')
            p2 = py.errorbar(self.midpts,(self.pofhists[it])/self.areas,xerr=self.widths,yerr=np.sqrt(self.pofhists[it])/self.areas,marker='o',ls='None')
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
        for it,hist in enumerate(self.agnhists):
            model = self.Naj[it]*self.psf+self.Ni[it]*self.iso + self.Nh*self.hmd
            #print self.Naj[it],self.psf,self.Ni[it],self.iso,self.Nh[0],self.hmd
            modelerrs = np.sqrt((self.Naj[it]*self.psf/sum(self.psf)*np.sqrt((self.psfe/self.psf)**2+(self.Naje[it]/self.Naj[it])**2))**2+(self.Nie[it]*self.iso)**2+(self.Nhe*self.hmd)**2)
            back = self.Ni[it]*self.iso
            backerrs = self.Nie[it]*self.iso
            p1 = py.errorbar(self.midpts,hist/self.areas,xerr=self.widths,yerr=np.sqrt(hist)/self.areas,marker='o',ls='None')
            p2 = py.errorbar(self.midpts,(model)/self.areas,xerr=self.widths,yerr=np.sqrt(model)/self.areas,marker='o',ls='None')

            names.append(self.agnlist[it]+' Data')
            pts.append(p1[0])
            names.append(self.agnlist[it]+' Model')
            pts.append(p2[0])

            if self.halomodel!='':
                #p4 = py.errorbar(self.midpts,self.Nh*self.hmd,xerr=self.widths,yerr=(self.Nhe*self.hmd)/self.areas,marker='o',ls='None')
                #names.append('Halo')
                #pts.append(p4[0])
                p5 = py.errorbar(self.midpts,(self.Nh*self.hmd+self.Ni[it]*self.iso)/self.areas,xerr=self.widths,marker='o',ls='None')
                names.append('Halo+Iso')
                pts.append(p5[0])
            else:
                p3 = py.errorbar(self.midpts,back/self.areas,xerr=self.widths,yerr=backerrs/self.areas,marker='o',ls='None')
                names.append(self.agnlist[it]+' Iso')
                pts.append(p3[0])
        py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$dN/d\theta^{2}$')
        py.xlim(min(self.midpts-self.widths)*(1-scale),max(self.midpts+self.widths)*(1+scale))
        if mi==1e40:
            mi = min(hist[amask]/self.areas[amask])
        py.ylim(0.25*mi,2*max(hist[amask]/self.areas[amask]))
        py.grid()
        py.legend(pts,names,loc=3)

        ############  PSF residuals plot  ###############
        ax = py.subplot(2,4,3)
        names = []
        pts = []
        ax.set_xscale("log", nonposx='clip')
        err = np.ones(self.nbins)
        p1 = py.errorbar(self.midpts,(self.psf-self.psfm)/self.psfm,xerr=self.widths,yerr=self.psfe/self.psfm,ls='None',marker='o')[0]
        cmask = self.psfe>0
        chisq = sum(((self.psf[cmask]-self.psfm[cmask])/self.psfe[cmask])**2)
        
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
                pt1 = py.errorbar(self.midpts,(psr-self.vij[it])/self.vij[it],xerr=self.widths,yerr=self.vije[it]/self.vij[it],ls='None',marker='o')[0]
                names.append(self.pulsars[it][0])
                pts.append(pt1)
                mask = psr>0
                up = max((psr[mask]-self.vij[it][mask]+self.vije[it][mask])/self.vij[it][mask])
                down = min((psr[mask]-self.vij[it][mask]-self.vije[it][mask])/self.vij[it][mask])
                ma = max(ma,up)
                ma = max(ma,abs(down))
                cmask = self.vije[it]>0
                chisq = sum(((psr[cmask]-self.vij[it][cmask])/self.vije[it][cmask])**2)
                tchisq = tchisq + '%s: %1.1f (%d)\n'%(self.pulsars[it][0],chisq,len(psr[cmask]))
            except:
                print 'Bad plotting' 
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
                pt1 = py.errorbar(self.midpts,(psr-model)/model,xerr=self.widths,yerr=modelerr/model,ls='None',marker='o')[0]
                names.append(self.pulsars[it][0])
                pts.append(pt1)
                mask = psr>0
                up = max((psr[mask]-model[mask]+modelerr[mask])/model[mask])
                down = min((psr[mask]-model[mask]-modelerr[mask])/model[mask])
                ma = max(ma,up)
                ma = max(ma,abs(down))
                cmask = modelerr>0
                chisq = sum(((psr[cmask]-model[cmask])/modelerr[cmask])**2)
                tchisq = tchisq + '%s: %1.1f (%d)\n'%(self.pulsars[it][0],chisq,len(psr[cmask]))
            except:
                print 'Bad plotting' 
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
            #try:
            model = self.Naj[it]*self.psf + self.Ni[it]*self.iso + self.Nh*self.hmd
            modelerr = np.sqrt((self.Naj[it]*self.psfe)**2 + (self.Naje[it]*self.psf)**2 + (self.Nie[it]*self.iso)**2 + (self.Nhe*self.hmd)**2)
            pt1 = py.errorbar(self.midpts,(agn-model)/model,xerr=self.widths,yerr=modelerr/model,ls='None',marker='o')[0]
            names.append(self.agnlist[it])
            pts.append(pt1)
            mask = agn>0
            up = max((agn[mask]-model[mask]+modelerr[mask])/model[mask])
            down = min((agn[mask]-model[mask]-modelerr[mask])/model[mask])
            ma = max(ma,up)
            ma = max(ma,abs(down))
            cmask = modelerr>0
            chisq = sum(((agn[cmask]-model[cmask])/modelerr[cmask])**2)
            tchisq = tchisq + '%s: %1.1f (%d)\n'%(self.agnlist[it],chisq,len(agn[cmask]))
            #except:
            #print 'Bad plotting' 
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
            print 'Unphysical Solution: %1.4f'%sterm

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
            print 'Unphysical Solution: %1.4f'%sterm

        #calculate gradient of background estimator analytically
        grad1 = -N-a*N+(4*a*(1+a)*b*N+2*(m*N-a*(b+n-m*N))*(N-a*(-N)))/(2*np.sqrt(sterm))    #psf derivative
        grad2 = -m-a*m+(4*a*(1+a)*b*m+2*(m*N-a*(b+n-m*N))*(m-a*(-m)))/(2*np.sqrt(sterm))    #number estimator derivative
        v = np.array([grad1,grad2])
        v = v/(2.*a*(1+a))
        return v

    ######################################################################
    #         Gradient of maximum likelihood from (3)                    #
    ######################################################################
    ## The gradient of the maximum likelhood estimator
    # @param params likelihood parameters defined by likelihood function
    # @param verb verbose slow output of gradient calculation
    def gradient(self,params,verb=False):

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
        Nh = params[nmu+psrs+agns+agns:nmu+psrs+agns+agns+1]

        mun = 1-np.sum(mi)
        mi = np.append(mi,[mun])
        grad = []
        if verb:
            print 'Obs\tmod\tfact\tnum\tacc'
            print '----------------------'

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
                    print npij[it2][it],denom,npij[it2][nmu],denomn,bpij[it2][it],backg,acc
                    t.sleep(0.25)

            #loop over AGN
            for it2 in range(agns):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh*self.hmd[it]                                             #ith model
                denomn = Naj[it2]*mi[nmu]+Ni[it2]*self.iso[nmu]+Nh*self.hmd[nmu]                                         #nth model
                fact0 = 0 if naij[it2][it]==0 else Naj[it2]*(naij[it2][it]/(denom) - 1.)                                    #ith contribution
                fact1 = 0 if naij[it2][nmu]==0 else -Naj[it2]*(naij[it2][nmu]/denomn - 1.)                                  #nth contribution

                #catch bad gradients
                if (denom <=0 and naij[it2][it]>0) or (denomn<=0 and naij[it2][nmu]>0):
                    flag=True
                acc = acc + fact0 + fact1
                if verb:
                    print naij[it2][it],denom,naij[it2][nmu],denomn,acc
                    t.sleep(0.25)

            if flag:
                grad.append(-np.Infinity)
            else:
                grad.append(-acc.item())
            if verb:
                print '----------------------'

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
                    print npij[it2][it],denom,mi[it],alpha*grad1,acc
                    t.sleep(0.25)
            if flag:
                grad.append(-np.Infinity)
            else:
                grad.append(-acc.item())
            if verb:
                print '----------------------'

        #AGN number estimator gradient
        for it2 in range(agns):
            acc = 0
            flag = False
            for it in range(self.nbins):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh*self.hmd[it]
                fact = 0 if naij[it2][it]==0 else mi[it]*(naij[it2][it]/denom - 1.)
                if (denom <=0 and naij[it2][it]>0):
                    flag=True
                acc = acc + fact
                if verb:
                    print naij[it2][it],denom,acc
                    t.sleep(0.25)
            if flag:
                grad.append(-np.Infinity)
            else:
                grad.append(-acc.item())
            if verb:
                print '----------------------'

        #Isotropic number estimator gradient for AGN
        for it2 in range(agns):
            acc = 0
            flag = False
            for it in range(self.nbins):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh*self.hmd[it]
                fact = 0 if naij[it2][it]==0 else self.iso[it]*(naij[it2][it]/denom - 1.)
                if (denom <=0 and naij[it2][it]>0):
                    flag=True
                acc = acc + fact
                if verb:
                    print naij[it2][it],denom,acc
                    t.sleep(0.25)
            if flag:
                grad.append(-np.Infinity)
            else:
                grad.append(-acc.item())
            if verb:
                print '----------------------'
        

        #Halo number estimator gradient for AGN
        acc = 0
        flag = False
        for it2 in range(agns):
            for it in range(self.nbins):
                denom = Naj[it2]*mi[it]+Ni[it2]*self.iso[it]+Nh*self.hmd[it]
                fact = 0 if naij[it2][it]==0 else self.hmd[it]*(naij[it2][it]/denom - 1.)
                if (denom <=0 and naij[it2][it]>0):
                    flag=True
                acc = acc + fact
                if verb:
                    print naij[it2][it],denom,acc
                    t.sleep(0.25)
        if flag:
            grad.append(-np.Infinity)
        else:
            grad.append(-acc.item())
        if verb:
            print '----------------------'
        
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
        err1 = so.fmin_powell(lambda x: abs(self.likelihood(self.minuit.params+x[0]*eig)-self.minuit.fval-0.5),
                              [self.minuit.params[num]*0.01],full_output=1,disp=disp)
        err1 = abs(err1[0])
        err2 = so.fmin_powell(lambda x: abs(self.likelihood(self.minuit.params-x[0]*eig)-self.minuit.fval-0.5),
                              [self.minuit.params[num]*0.01],full_output=1,disp=disp)
        err2 = abs(err2[0])

        #try to catch badly formed likelihood surfaces
        if self.likelihood(self.minuit.params-err1*eig)==np.Infinity:
            return err2*err2
        if self.likelihood(self.minuit.params+err2*eig)==np.Infinity:
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
        rows = int(np.sqrt(len(self.minuit.params)))+1
        #go through all parameters
        #find the quadratic form of the likelihood surface
        #determine the increase in the variance from covariance of parameters
        for it in range(len(self.minuit.params)):

            if it!=num:
                #py.subplot(rows,rows,it+1)
                eigx = np.zeros(len(self.minuit.params))
                eigx[num]=err
                eigy = np.zeros(len(self.minuit.params))
                eigy[it]=np.sqrt(self.errs2[it])

                #calulate likelihood along ring around maximum (1-sigma in likelihood)
                px = (self.likelihood(self.minuit.params+eigx)-self.minuit.fval)
                pxpy = (self.likelihood(self.minuit.params+(eigx+eigy)/rt)-self.minuit.fval)
                py1 = (self.likelihood(self.minuit.params+eigy)-self.minuit.fval)
                mxpy = (self.likelihood(self.minuit.params+(-eigx+eigy)/rt)-self.minuit.fval)
                mx = (self.likelihood(self.minuit.params-eigx)-self.minuit.fval)
                mxmy = (self.likelihood(self.minuit.params+(-eigx-eigy)/rt)-self.minuit.fval)
                my = (self.likelihood(self.minuit.params-eigy)-self.minuit.fval)
                pxmy = (self.likelihood(self.minuit.params+(eigx-eigy)/rt)-self.minuit.fval)
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

                    #print xmin,ymin
                    #t.sleep(0.25)
                except:
                    pass
                    #print 'Caught poorly formed quad surface'
        #py.savefig('likes%d.png'%num)

        return (ma*err)**2

    def tspval(self):
        pval = 1-sst.stats.chisqprob(self.TS,2)/2.
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
def test(bins=8,ctype=0,emin=1000,emax=1778,days=30,irf='P6_v3_diff',maxr=-1,sel='[0]',agnlis=['agn-psf-study-bright'],double=False,ctlim=[0.4,1.0]):
    psf = CALDBPsf(CALDBManager(irf=irf))
    ebar = np.sqrt(emin*emax)
    psrs = ''

    if maxr<0:
        maxr = psf.inverse_integral(ebar,ctype,99.5)*1.5             #use 1.5 times the 99.5% containment as the maximum distance
    cl = CombinedLike(irf=irf,mode=-1,pulsars = eval('pulsars'+sel+''),agnlist=agnlis,verbose=True,ctmin=ctlim[0],ctmax=ctlim[1])
    cl.loadphotons(0,maxr,emin,emax,239557417,239517417+days*86400,ctype)
    cl.bindata(bins)
    f0 = cl.fit()
    params = cp.copy(cl.fitpsf(double))
    print str(cl)
    for psr in cl.pulsars:
        psrs = psrs + '%s_'%psr[0]
    cl.makeplot('figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%s%s'%(emin,emax,ctype,maxr,bins,psrs,('').join(cl.agnlist)))
    # cl
    return params,cl

def test2(bins=8,ctype=0,emin=1000,emax=1778,days=730,irf='P6_v3_diff',maxr=-1,sel='[0]',agnlis='agn-psf-study-bright',model='Gaussian'):
    psf = CALDBPsf(CALDBManager(irf=irf))
    ebar = np.sqrt(emin*emax)
    if maxr<0:
        maxr = psf.inverse_integral(ebar,ctype,99.5)*1.5             #use 1.5 times the 99.5% containment as the maximum distance
    cl = CombinedLike(irf=irf,mode=-1,pulsars = eval('pulsars'+sel+''),agnlist=[agnlis],verbose=False)
    cl.loadphotons(0,maxr,emin,emax,239557417,239517417+days*86400,ctype)
    cl.bindata(bins)
    f0 = cl.fit()
    print cl
    fitpars = cl.fitpsf()
    psrs = ''
    psfphotons = cl.Naj[0]
    for psr in cl.pulsars:
        psrs = psrs + '%s_'%psr[0]
    r05 = psf.inverse_integral(ebar,ctype,5.)
    r34 = psf.inverse_integral(ebar,ctype,34.)
    r68 = psf.inverse_integral(ebar,ctype,68.)
    r95 = psf.inverse_integral(ebar,ctype,95.)
    r99 = psf.inverse_integral(ebar,ctype,99.)
    #testr = (r68+r99)/2.
    nps = 10
    pts = np.arange(0,nps+1,1)
    pts = (pts*(r99-r05)/nps+r05)
    npars = len(cl.params)

    uplims2 = []
    uplims1 = []
    detect = []
    TSs = []
    of = open('figures/uplimemi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%s%s_%s.txt'%(emin,emax,ctype,maxr,bins,psrs,cl.agnlist[0],model),'w')
    print >>of,'%1.3f'%f0
    print >>of,'size(rad)\tmaxL\tNha \tcl68\tcl95\tTS'
    for pt in pts:
        cl.fit(halomodel=model,haloparams=[pt/rd,ebar,ctype,fitpars[1]/rd,fitpars[2]])
        print cl
        lmax = cl.minuit.fval
        cl.TS = 2*(f0-lmax)
        TSs.append(cl.TS)
        #norm = si.quad(lambda x: np.exp(lmax-cl.profile(npars-1,x)),0,cl.Naj[0])[0]
        #print norm

        total = cl.Naj[0]+cl.Nh
        bestp = so.fmin_powell(lambda x: abs(cl.profile(npars-1,x[0])-lmax-0.5) if x[0]>cl.Nh else 1e40,[cl.Nh+total/10.],disp=0,full_output=1)

        fac = bestp[0]#cl.Nh/30. if cl.Nhe>cl.Nh else cl.Nhe
        #xr = np.array([x for x in (cl.Nh+cl.Naj[0]+np.arange(-10,11,1)*fac) if (x>0) ])
        #if len(xr)<10:
        xr = np.arange(0,11,1)
        if (cl.Nh+(-5)*(fac-cl.Nh)/2.)<0:
            xr = fac*xr
        else:
            xr = cl.Nh+(xr-5)*(fac-cl.Nh)/2.
        fx = np.array([lmax-cl.profile(npars-1,x) for x in xr])
        print xr
        print fx
        cl1 = BayesianLimit(xr,fx).getLimit(0.32)
        cl2 = BayesianLimit(xr,fx).getLimit(0.05)
        detect.append(cl.Nh/(psfphotons))
        uplims2.append(cl2[0]/(psfphotons))
        uplims1.append(cl1[0]/(psfphotons))
        print cl2[0],psfphotons,cl2[0]/(psfphotons)
        #py.clf()
        #py.plot(xr,np.exp(fx),'ro')
        #py.savefig('figures/proftest%1.4f.png'%pt)
        cl.makeplot('figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s'%(emin,emax,ctype,maxr,bins,pt,psrs,cl.agnlist[0],model))
        print >>of,'%1.5f\t%1.3f\t%1.1f\t%1.1f\t%1.1f\t%1.3f'%(pt/rd,lmax,cl.Nh,cl1,cl2,cl.TS)
        #sig1 = si.quad(lambda x: np.exp(lmax-cl.profile(npars-1,x))/norm,0,cl.Nh+cl.Nhe)[0]
        #print sig1
        #sig2 = si.quad(lambda x: np.exp(lmax-cl.profile(npars-1,x))/norm,0,cl.Nh+2*cl.Nhe)[0]
        #print sig2
        #sig3 = si.quad(lambda x: np.exp(lmax-cl.profile(npars-1,x))/norm,0,cl.Nh+3*cl.Nhe)[0]
        #print sig3
        #print norm,sig1,sig2,sig3,cl.Nh,cl.Nhe

    py.figure(10,figsize=(8,8))
    py.clf()
    maxts = max(TSs)*1.1
    p2=py.plot(pts,np.array(uplims2)*maxts/1.1,'v')
    p1=py.plot(pts,np.array(uplims1)*maxts/1.1,'v')
    p5=py.plot(pts,np.array(detect)*maxts/1.1,'o')
    p6=py.plot(pts,TSs,'d')
    p3=py.plot([r68,r68],[0,maxts],'r--')
    p4=py.plot([r95,r95],[0,maxts],'r-.')
    p7=py.plot([0,max(pts)*1.1],[1,1],'b--')
    p8=py.plot([0,max(pts)*1.1],[4,4],'b-.')
    py.grid()
    py.xlabel('%s halo size (deg)'%model)
    py.ylabel('%s halo fraction'%model)
    py.ylim(0,maxts)
    py.xlim(0,max(pts)*1.1)
    prop = mpl.font_manager.FontProperties(size=9) 
    py.legend((p5,p1,p2,p3,p4,p6,p7,p8),(r'hfract',r'68% CL',r'95% CL','%s R68'%irf,'%s R95'%irf,'TS',r'$1 \sigma$',r'$2 \sigma$'),bbox_to_anchor=(1.1, 1.0),prop=prop)
    py.savefig('figures/uplimemi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%s%s_%s.png'%(emin,emax,ctype,maxr,bins,psrs,cl.agnlist[0],model))
    of.close()
    #for pt in pts:
    #    tval = cl.profile(npars-1,pt)
    #    print pt,(lmax-tval)
    """cl.makeplot('figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%s%s'%(emin,emax,ctype,maxr,bins,psrs,cl.agnlist[0]))
    f1=[]
    cl.printResults()
    halorange = np.arange(0.05,maxr,maxr/10.)
    likes = []
    nest = []
    mi = 1e40
    midx = 0
    model = 'Gaussian'
    for num,it in enumerate(halorange):
        likes.append(cl.fit(halomodel=model,haloparams=[it/rd,ebar,ctype]))
        nest.append(cl.Nh)
        print it,nest[num],likes[num]
        if likes[num]<mi:
            mi = likes[num]
            midx = num

    #minu = Minuit(lambda x:cl.fit(halomodel=model,haloparams=[x[0]/rd,ebar,ctype]),[halorange[midx]],limits=[[0,maxr]],steps=[0.0001],printMode=0,strategy=0)
    #minu.minimize()
    #minu.errors()
    bestp = so.fmin_powell(lambda x:cl.fit(halomodel=model,haloparams=[x[0]/rd,ebar,ctype]),(halorange[midx]),full_output=1)
    cl.printResults()
    print 'Halo size was: %1.3f'%bestp[0]
    print 'TS was: %1.2f'%(2*(f0-bestp[1]))
    cl.TS = 2*(f0-bestp[1])
    #for it in (np.arange(1.,10.,1.)/10.):
    #it=bestp[0]
    #f1.append(cl.fit(halomodel='Halo',haloparams=[it/rd]))
    cl.likelihood(cl.minuit.params,True)
    cl.makeplot('figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s'%(emin,emax,ctype,maxr,bins,bestp[0],psrs,cl.agnlist[0],model))"""

def runall(ctype,bins=12):
    emins = [100,178,316,562,1000,1778,3162,5623,10000,17783,31623]
    emaxs = [178,316,562,1000,1778,3162,5623,10000,17783,31623,100000]
    ctypes = [0,1]
    irfs = ['P6_v11_diff']#,'P6_v3_diff']
    sels = ['[0:2]']
    #models=['Halo','Gaussian']
    #lists = ['agn_redshift2_lo_bzb','agn_redshift2_hi_bzb']#'tevbzb','1es0229p200','1es0347-121','1es1101-232']
    lists = ['agn-psf-study-bright']
    #for model in models:
    ff = open('figures7c/output%d.txt'%ctype,'w')
    for it,emin in enumerate(emins):
        #for ctype in ctypes:
        for sel in sels:
            for irf in irfs:
                if emin<10000:
                    if (emin+ctype*179)<1778:
                        lis = []
                    else:
                        lis = lists
                    pars = test(bins=bins,ctype=ctype,emin=emin,emax=emaxs[it],days=730,irf=irf,sel=sel,agnlis=lis,double=False)
                else:
                    sl = StackLoader(lis=lists[0],irf=irf,srcdir=srcdir,useft2s=False)
                    sl.files = glob.glob(agndir+'*.fits')
                    sl.loadphotons(0,3,emin,emaxs[it],0,239557417+730*86400,ctype)
                    sl.solvepsf()
                    pars = [sl.ebar,sl.sigma*rd,sl.gamma]#,sl.sigma2*rd,sl.gamma2,sl.Npsf/(sl.Npsf+sl.Npsf2)]
                print >>ff,'%1.6f\t%1.6f\t%1.6f'%(pars[0],pars[1],pars[2])
    ff.close()
    
def runconvolution(model=['CDisk','CHalo']):
    machines = 'tev01 tev02 tev03 tev04 tev05 tev06 tev07 tev08 tev09 tev10 tev11'.split()
    engines = 2
    setup_string = 'import uw.stacklike.binned as ub;reload(ub);from uw.stacklike.binned import *'
    emins = [1000,1778,3162,5623,10000,17783,31623]#[100,178,316,562,
    emaxs = [1778,3162,5623,10000,17783,31623,100000]#178,316,562,1000,
    tasks = ['ub.test2(bins=12,ctype=%d,emin=%d,emax=%d,days=730,irf=\'P6_v3_diff\',maxr=-1,sel=\'[0:2]\',agnlis=\'%s\',model=\'%s\')'%(y,emins[x],emaxs[x],lis,mod) for x in range(len(emins)) for y in range(2) for lis in ['agn_redshift2_lo_bzb','agn_redshift2_hi_bzb'] for mod in model]
    print tasks
    ua.setup_mec(engines=engines,machines=machines)
    t.sleep(15)
    logfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/python/mec.log','w')
    at = ua.AssignTasks(setup_string,tasks,log=logfile,timelimit=10000,progress_bar=True,ignore_exception=True)
    at(30)
    ua.kill_mec()

def MCrun():
    for i in range(100):
        cl.ponhists = sst.poisson.rvs(t1,1)
        cl.pofhists = sst.poisson.rvs(t2,1)
        cl.agnhists = sst.poisson.rvs(t3,1)

        f0 = cl.fit(halomodel='')
        f1 = cl.fit(halomodel=model,haloparams=[pt/rd,ebar,ctype,fitpars[1]/rd,fitpars[2]])
        TS = 2*(f0-f1)
        cl.TS=TS
        TSs.append(TS)
        cl.makeplot('figures/emi%1.0f_ema%1.0f_ec%1.0f_roi%1.2f_bins%1.0f_%1.1f%s%s_%s%s'%(emin,emax,ctype,maxr,bins,pt,psrs,cl.agnlist[0],model,i))
        print TS