"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/stacklike/stacklike.py,v 1.19 2012/07/17 19:14:24 mar0 Exp $
author: M.Roth <mar0@u.washington.edu>
"""

import numpy as np
import pylab as py
import skymaps as s
import pointlike as pl
import astropy.io.fits as pf
import glob as glob
import scipy.integrate as si
import scipy.optimize as so
import os as os
from uw.like import pypsf
from uw.stacklike.angularmodels import *
from uw.stacklike.CLHEP import HepRotation,Hep3Vector,Photon
from uw.stacklike.stcaldb import IrfLoader
import uw.thb_roi.roi_factory as uf
import copy as cp
import uw.like.SpatialModels as us
import uw.utilities.convolution as uc
from uw.utilities.minuit import Minuit
import matplotlib.font_manager
from uw.like.pycaldb import CALDBManager
from time import sleep
import copy as cp
from uw.stacklike import dataset

#################################################        STACKLIKE         ####################################################
#
#                           Stacklike is an analysis suite which models angular distributions of photons
#                           by stacking sources. The angular distributions are then fit to the models
#                           by maximizing the extended likelihood with respect model parameters. This
#                           includes the number of photons associated with each submodel (eg PSF parameters
#                           and # of photons from the PSF). The submodels are defined in angularmodels.py,
#                           and are combined into a composite model where the overall likelihood is calculated.
#                           This makes it natural to add additional model components with different angular
#                           distributions.
#
#                           files:
#                           stacklike.py     - this file
#                           angularmodels.py - contains all models, composite models, and fitters
#                           CLHEP.py         - contains implementation of matrix and vector classes for transforming 
#                                              photon positions 
#                           
#


################################################# HOWTO MAKE SOURCE LISTS #####################################################
#
#  simply make a new text file: 'test.txt'
#
#  the first line of the file will be ignored so feel free to make a header:
#
#         #name    ra    dec
#         source1  128.8379  -45.1763
#         source2  0.0000    0.0000
#         .........................
#         (und so weiter)
#



################################################## ENVIRONMENT SETUP        ###################################################
debug=False
rd = 180./np.pi        #to radians
monthlist = ['aug2008','sep2008','oct2008','nov2008','dec2008','jan2009','feb2009','mar2009','apr2009'\
    ,'may2009','jun2009','jul2009','aug2009','sep2009','oct2009','nov2009','dec2009','jan2010','feb2010','mar2010'\
    ,'apr2010','may2010','jun2010','jul2010']

environ = 'TEV'

if environ=='FERMITS':
    basedir = r'd:\common\mar0\stacklike/'
    CALDBdir = r'd:\fermi/CALDB/v1r1/CALDB/data/glast/lat'
    datadir = r'd:\fermi\data\flight/'                             #directory for FT1 files
    ft2dir = r'd:\fermi\data\flight/'                              #directory for FT2 files
    srcdir = r'd:\common\mar0\sourcelists/'                        #directory for source lists
if environ=='TEV':
    basedir = r'/phys/groups/tev/scratch1/users/Fermi/'
    CALDBdir = basedir+'packages/ScienceTools-09-20-00/irfs/caldb/CALDB/data/glast/lat'
    datadir = basedir+'mar0/data/7.3src/'                             #directory for FT1 files
    ft2dir = basedir+r'mar0/data/'                              #directory for FT2 files
    srcdir = r'/phys/users/mar0/sourcelists/'                        #directory for source lists
files = []                                                     #default list of FT1 file names (minus '-ft1.fits')
ft2s = []                                                      #default list of FT2 file names (minus '-ft2.fits')
for month in monthlist:
    files.append('%s'%month)
    ft2s.append('%s'%month)
pulsar=True


###################################################  START STACKLOADER CLASS ##################################################

##  StackLoader class
#   
#   loads a set of photons with instrument coordinate information and stacks them based on a list of sources
#
#   usage:
#
#   sl = Stackloader(args1)
#   sl.loadphotons(args2)
#   
#
#   StackLoader contains methods to solve alignment, PSF, and halo parameters
#   ***************************************************************************
#   solverot() - outputs the boresight alignment in arcseconds for Rx(x)*Ry(y)*Rz(z)
#
#   solvepsf() - fits the parameters of a single King function in a uniform backround
#
#   solvehalo() - fits the parameters of a gaussain halo with a PSF and uniform backround
#
#   solvedoublepsf() - fits the parameters of two King functions in a uniform backround
#   !!!Warning: large degeneracy in on-orbit data, only works with very low uniform backround!!!
class StackLoader(object):


    ##  Constructor - sets up sources and pointing history file
    #   @param lis text file with all the source names, ra's, and dec's
    #   @param tev University of Washington flag for TEV machines
    #   @param rot optional rotation information input beforehand (in radians)
    #   @param files list of filenames in the form ['name1'] corresponding to 'name1-ft1.fits' and 'name1-ft2.fits'
    #   @param CALDBdr CALDBdir pointlike dir (should be '......./glast/lat/')
    #   @param ft2dr ft2 file directory
    #   @param datadr ft1 file directory
    #   @param srcdr source list directory (ascii file with 'name    ra     dec', first line is skipped)
    #   @param irf response function of the form 'P?_v?'
    def __init__(self,**kwargs):
        self.name=''
        self.rot=[0,0,0]
        self.binfile=''
        self.irf = 'P6_v3_diff'
        self.ctmin=0.4
        self.ctmax=1.0
        self.quiet=True
        
        self.__dict__.update(kwargs)

        self.dsel = eval('dataset.%s'%(self.name))
        
        self.irf = self.dsel.irf
        self.rot = self.dsel.rot

        if not self.quiet:
            print ('')
            print ('**********************************************************')
            print ('*                                                        *')
            print ('*                       STACKLIKE                        *')
            print ('*                                                        *')
            print ('**********************************************************')
            print ('Using irf: %s'%self.irf)
            print ('Using list: %s.txt'%self.name)
            print ('Using boresight alignment (in arcsec): Rx=%1.0f Ry=%1.0f Rz=%1.0f'%(self.rot[0]*rd*3600,self.rot[1]*rd*3600,self.rot[2]*rd*3600))

        self.rotang = self.rot
        self.rot = HepRotation(self.rot,False)
        self.atb = []
        self.ds = []
        self.excess = []
        self.diffs =[]
        self.aeff = []
        self.photons = []
        self.offsrc = []
        self.photoncount=0
        self.Npsf =0
        self.Nback =0
        self.ebar=0
        self.cls=0
        self.srcs = self.dsel.srcs
        self.files = self.dsel.ft1s
        self.ft2s = self.dsel.ft2s
        self.useft2s = self.dsel.useft2s
        self.phasecut = self.dsel.phase

    ##  loads photons from FT1 file near sources  
    #   @param minroi min ROI in degrees
    #   @param maxroi max ROI in degrees
    #   @param emin minimum energy in MeV
    #   @param emax maximum energy in MeV
    #   @param start start time (MET)
    #   @param stop end time (MET)
    #   @param cls conversion type: 0=front,1=back,-1=all
    def loadphotons(self,minroi,maxroi,emin,emax,start,stop,cls):
        self.start = start
        self.stop = stop
        self.atb = []
        self.ebar=self.ebar*len(self.photons)
        self.emin=emin
        self.emax=emax
        self.minroi = minroi
        self.maxroi = maxroi
        self.cls = cls
        self.eave = np.sqrt(self.emin*self.emax)
        if self.dsel.binfile is not None:
            print ('Applying masks to data')
            print ('**********************************************************')
            print ('%1.0f < Energy < %1.0f'%(self.emin,self.emax))
            print ('Event class = %1.0f'%(cls))
            print ('**********************************************************')

            self.bpd = s.BinnedPhotonData(self.datadir+self.binfile)
            self.ebar = 0
            for bnd in self.bpd:
                if cls ==0:
                    catchbad = bnd.event_class()==0 or bnd.event_class()==-2147483648
                else:
                    catchbad = bnd.event_class()==1
                if bnd.emax()>self.emin and bnd.emin()<self.emax and (catchbad or cls==-1):
                    ebar = np.sqrt(bnd.emin()*bnd.emax())
                    for src in self.srcs:
                        wsdl = s.WeightedSkyDirList(bnd,src,self.maxroi/rd)
                        sds = len(wsdl)
                        for wsd in wsdl:
                            self.photons.append(Photon(wsd.ra(),wsd.dec(),ebar,0,cls,[],[],src,weight=wsd.weight()))
                            self.ebar = self.ebar + wsd.weight()*ebar
            pcnts = sum([x.weight for x in self.photons])
            self.photonscount = self.photoncount+pcnts
            if len(self.photons)==0:
                print ('No photons!')
                raise 'No photons!'
            else:
                self.ebar = self.ebar/pcnts
            print ('%d photons remain'%pcnts)
            print ('**********************************************************')
        else:
            print ('Applying masks to data')
            print ('**********************************************************')
            print ('%1.0f < Energy (MeV) < %1.0f'%(self.emin,self.emax))
            print ('%1.2f < Separation (deg) < %1.2f'%(self.minroi,self.maxroi))
            print ('%1.0f < Time < %1.0f'%(start,stop))
            print ('Event class = %1.0f'%(cls))
            print ('%1.2f < costh < %1.2f'%(self.ctmin,self.ctmax))
            print ('**********************************************************')

            """#go through each fits file, mask unwanted events, and setup tables
            tcuts = np.array([0,0,0,0,0,0,0,0,0])
            for j,ff in enumerate(self.files):
                if not self.quiet:
                    print (ff)
                tff = pf.open(ff)
                ttb = tff[1].data
                ttb,ttcut = self.mask(ttb,self.srcs,emin,emax,start,stop,self.maxroi,cls)
                tcuts=tcuts+ttcut
                if len(ttb)>0:
                    self.atb.append(ttb)
                else:
                    if self.useft2s:
                        del self.ft2s[j]
                    del self.files[j]
                    del ttb

            print ('Photon pruning information')
            print ('%d photons available'%tcuts[0])
            print ('---------------------------------------')
            print ('%d cut by time range'%tcuts[1])
            print ('%d cut by energy range'%tcuts[2])
            print ('%d cut by instrument theta'%tcuts[3])
            print ('%d cut by zenith angle'%tcuts[4])
            print ('%d cut by conversion type'%tcuts[5])
            print ('%d cut by proximity to sources'%tcuts[6])
            print ('%d cut by phase'%tcuts[7])
            print ('---------------------------------------')
            print ('%d photons remain'%tcuts[8])
            print ('**********************************************************')"""


            #go through each table and contruct Photon objects
            for j,ff in enumerate(self.files):
                tff = pf.open(ff)
                tb = tff[1].data
                hdr = tff[1].header
                if float(hdr['TSTART'])>stop or float(hdr['TSTOP'])<start:
                    if not self.quiet:
                        print ('Skipped %s [%d,%d], [%d,%d]: Out of time range'%(ff,hdr['TSTART'],hdr['TSTOP'],start,stop))
                    continue
                try:
                    tb.field('THETA');tb.field('RA');tb.field('TIME');
                except:
                    continue
                tb,ttcut = self.mask(tb,self.srcs,emin,emax,start,stop,self.maxroi,cls)
                try:
                    tb.field('CTBCLASSLEVEL')
                    pass7=False
                except:
                    pass7=True
                if not self.quiet:
                    print ('Examining %d events'%len(tb))
                if self.useft2s:
                    if not self.quiet:
                        print ('    *Loading pointing history')
                    try:
                        phist = pl.PointingHistory(self.ft2s[j])
                    except:
                        print ('Failed pointing history file: %s'%self.ft2s[j])
                        phist = None
                else:
                    phist=None
                if not self.quiet:
                    print ('    *Loading events')
                try:
                    gti = s.Gti(self.files[j])
                    print ('Loaded GTI: %d - %d'%(gti.minValue(),gti.maxValue()))
                except:
                    print ('Failed GTI in file: %s'%self.ft2s[j])
                    gti = None
                #iterate over events for photons

                if gti is not None and ((phist is not None and self.useft2s) or not self.useft2s):
                    pcut = 0
                    for k in range(len(tb)):
                        event = tb[k]
                        s_ra,s_dec=float(event.field('RA')),float(event.field('DEC'))
                        sd = s.SkyDir(s_ra,s_dec)#Hep3Vector([s_ra/rd,s_dec/rd])
                        rsrc = self.srcs[0]
                        diff = 1e40
                        if pass7:
                            if (event.field('EVENT_CLASS') & self.dsel.eventclass) == 0:
                                pcut+=1
                                continue
                        #associate with a particular source
                        for src in self.srcs:
                            diff = sd.difference(src)*rd
                        
                            #if photons is not too far away from any source
                            #add it to the list
                            if diff<self.maxroi and diff>self.minroi:
                                time = event.field('TIME')
                                xv = []
                                zv = []
                                if self.useft2s:
                                    pi = phist(time)
                                    xv = Hep3Vector([])
                                    zv = Hep3Vector([])
                                    xax = pi.xAxis()
                                    zax = pi.zAxis()
                                    zen = pi.zenith()
                                    xv(xax)
                                    zv(zax)

                                accept = gti.accept(time)
                                if accept:
                                    photon = Photon(s_ra,s_dec,event.field('ENERGY'),time,event.field('CONVERSION_TYPE'),xv,zv,src,ct=np.cos(event.field('THETA')/rd))
                                    #print ('%1.2e'%(np.cos(zv.dir().difference(sd))-photon.ct))
                                    #t.sleep(0.1)
                                    self.ebar = self.ebar+event.field('ENERGY')
                                    self.photons.append(photon)
                                    self.photoncount = self.photoncount + 1
                print (pcut)
                if self.useft2s:
                    del phist #free up memory from pointing history object
                tff.close()
                del tb
                del gti
            #del self.atb
            if len(self.photons)>0:
                #print ('%d of %d were in a GTI'%(len(self.photons),tcuts[8]))
                self.ebar = self.ebar/len(self.photons)  #calculate mean energy of photons
            else:
                print ('No Photons!')

                #raise 'No Photons!'

    ## Bins all of the angular separations together in separation
    def bindata(self):
        self.getds()
        nbins = 2*max(1,int(np.sqrt(self.photoncount)))
        print (nbins)
        hist = np.histogram(np.log10(self.ds),bins=nbins)
        diffs = np.sqrt((10**(2*hist[1][1:])+10**(2*hist[1][:len(hist[1])-1]))/2.)*rd   #average separation of each bin is determined by area
        
        #free up memory
        try:
            del self.photons
        except:
            pass
        try:
            del self.ds
        except:
            pass
        try:
            del self.atb
        except:
            pass
        self.photons=[]
        self.ds=[]

        #single source
        self.srcs=[s.SkyDir(0,0)]
        for it2,nums in enumerate(hist[0]):
            if nums>0:
                self.photons.append(Photon(diffs[it2],0,self.ebar,0,self.cls,[],[],self.srcs[0],weight=nums))
        self.bin=True
        self.photoncounts = sum([photon.weight for photon in self.photons])
        self.getds()

    ## returns number of photons above isotropic component
    #  @param exp number estimator for the isotropic component
    def solvenonback(self,exp):
        self.getds()
        backd = rd*rd*exp/(self.maxroi*self.maxroi-self.minroi*self.minroi)
        hexp = np.histogram(self.ds,bins=200)
        ntot = 0
        for it,nbin in enumerate(hexp[0]):
            dA = (hexp[1][it+1]*hexp[1][it+1]-hexp[1][it]*hexp[1][it])
            ntot = ntot + nbin - backd*dA
        return ntot

    def estback(self,frac):
        #Try to estimate the number of Isotropic photons
        #by looking at density near tails
        self.getds()
        area = 1-frac*frac
        density = len(self.ds[self.ds>(frac*self.maxroi/rd)])
        denback = density*(self.maxroi*self.maxroi-self.minroi*self.minroi)/((self.maxroi**2)*area)
        return denback

    def excessdist(self,perc=0.75):
        back = self.estback(perc)
        self.getds()
        self.ds = np.sort(self.ds)
        excess = np.array([(i+1)-((ds*rd)**2)*back/(self.maxroi*self.maxroi-self.minroi*self.minroi) for i,ds in enumerate(self.ds)])
        self.excess = excess/max(excess)
        return excess/max(excess)

    def rcontainext(self,frac,perc=0.75):
        if True:#len(self.excess)==0 or len(self.ds)==0:
            self.excessdist(perc)
        best = 1.0
        for it,exc in enumerate(self.excess):
            if abs(exc-frac)<best:
                dis = self.ds[it]
                best=abs(exc-frac)
        return dis*rd

    def modelprofile(self,which,val):
        pars = cp.copy(self.cm.minuit.params)
        pars[which]=val
        fixed = [i==which for i in range(len(pars))]
        prof = Minuit(lambda x: self.cm.extlikelihood(x,self.photons),pars,fixed=fixed,printMode=-1)
        prof.minimize()
        return prof.fval

    ## estimates Isotropic component from PSF
    #  @param func optional custom background component to fit [[angular separations],[counts]]
    #  @param free free or fix model number estimators
    #  @param exp specify number estimators
    def solveback(self,func=[[],[]],free=[True,True],exp=[],model_par=[]):

        #denback = self.estback(0.75)
        #if denback>self.photoncount:
        denback =self.photoncount/2
        
        # setup default or aligned PSF
        if model_par==[]:
            psf = self.getpsf()
            if self.useft2s:
                psf2 = PSFAlign(lims=[self.minroi/rd,self.maxroi/rd],model_par=self.rotang,free=[False,False,False])
                psf2.PSF = cp.copy(psf)
                psf = cp.copy(psf2)

        # setup PSF or aligned PSF with specified parameters
        else:
            if len(model_par)==2:
                psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=model_par,free=[False for x in range(len(model_par))])
            else:
                psf = PSFDouble(lims=[self.minroi/rd,self.maxroi/rd],model_par=model_par,free=[False for x in range(len(model_par))])
            if self.useft2s:
                psf2 = PSFAlign(lims=[self.minroi/rd,self.maxroi/rd],model_par=self.rotang,free=[False,False,False])
                psf2.PSF = cp.copy(psf)
                psf = cp.copy(psf2)

        # setup the models
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)

        #if custom background is specified, add it
        if len(func[0])!=0:
            bck2 = Custom(lims=[self.minroi/rd,self.maxroi/rd],func=func)
            cm.addModel(bck2)

        #check to see if number estimators are specified
        if exp==[]: 
            if len(cm.models)==2:
                exp=[self.photoncount-denback,denback]
            else:
                exp=[self.photoncount-denback,0.,denback]
        if len(free)!=len(cm.models):
            free=[True,True,True]

        #maximize the likelihood
        cm.fit(self.photons,mode=0,quiet=self.quiet,free=free,exp=exp)

        #get parameters and errors
        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Nbacke = np.sqrt(self.errs[1][1])
        phs = (self.Npsf+self.Nback)
        tr = sum([self.errs[i][i] for i in range(len(self.errs))])
        f1 = self.Npsf/phs
        f1e = f1*np.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        f2 = self.Nback/phs
        f2e = f2*np.sqrt(self.errs[1][1]/(self.Nback**2)+tr/(phs**2))
        nonback = self.solvenonback(self.Nback)
        if not self.quiet:
            print ('**********************************************************')
            print ('Npsf  = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100))
            print ('Nback = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f2*100,f2e*100))
            print ('Nsource = %1.0f'%nonback)
            print ('**********************************************************')
        if len(cm.models)==2:
            if len(psf.model_par)==2:
                self.il=cm.extlikelihood([self.Npsf,self.Nback,psf.model_par[0],psf.model_par[1]],self.photons)
            if len(psf.model_par)==3:
                self.il=cm.extlikelihood([self.Npsf,self.Nback,psf.model_par[0],psf.model_par[1],psf.model_par[2]],self.photons)
            if len(psf.model_par)==5:
                self.il=cm.extlikelihood([self.Npsf,self.Nback,psf.model_par[0],psf.model_par[1],psf.model_par[2],psf.model_par[3],psf.model_par[4]],self.photons)
        else:
            self.Ncust = cm.minuit.params[2]
            self.Ncuste = np.sqrt(self.errs[2][2])
            self.il=cm.extlikelihood([self.Npsf,self.Nback,self.Ncust,psf.model_par[0],psf.model_par[1]],self.photons)
        self.cm = cm
        return cm.minuit.fval

    def rcontainlike(self,ct=[]):
        if len(self.photons)==0:
            return 0
        if ct==[]:
            psf = self.getpsf()
        else:
            psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.001,2.25],free=[True,True])
            psf.fromcontain([ct[0],ct[1]],[0.68,0.95])

        return self.solveback(model_par=psf.model_par,exp=[self.Npsf,self.Nback])

    ## returns the TS for a point source in a uniform background
    #  @param ct optionally define the PSF by [R68,R95] in radians, default is the self.irf
    def TS(self,ct=[]):

        if len(self.photons)==0:
            return 0
        if ct==[]:
            psf = self.getpsf()
        else:
            psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.001,2.25],free=[True,True])
            psf.fromcontain([ct[0],ct[1]],[0.68,0.95])

        self.solveback()#model_par=psf.model_par)
        cm = CompositeModel()

        back = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm.addModel(psf)
        cm.addModel(back)
        f0 = cm.extlikelihood([0,self.Npsf+self.Nback,psf.model_par[0],psf.model_par[1]],self.photons)
        f1 = cm.extlikelihood([self.Npsf,self.Nback,psf.model_par[0],psf.model_par[1]],self.photons)
        return -2*(f1-f0)

    ## Finds boresight alignment solution
    def solverot(self):

        #estimate Isotropic component
        self.solveback()

        psfa = PSFAlign(lims=[self.minroi/rd,self.maxroi/rd],free=[True,True,True],ebar=self.ebar,ec=self.cls,irf=self.irf)
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm = CompositeModel()
        cm.addModel(psfa)
        cm.addModel(bck)
        cm.fit(self.photons,free=[True,True],exp=[self.Npsf,self.Nback],mode=1,quiet=self.quiet)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Nbacke = np.sqrt(self.errs[1][1])
        self.params=[x*rd*3600 for x in cm.minuit.params[2:]]
        self.errors=[np.sqrt(self.errs[x+2][x+2])*rd*3600 for x in range(3)]
        self.il=cm.extlikelihood([self.Npsf,self.Nback,0,0,0],self.photons)
        for x in range(3):
            print ('R%d %3.0f +/-%3.0f arcsec'%(x+1,self.params[x],self.errors[x]))
        try:
            TS = -2*(fl-self.il)
            print ('Significance of rotation was TS = %d'%TS)
        except:
            print ('Cannot determine improvement')
        print ('Expected %d photons, got %d ( %d (%d) + %d (%d) )'%(len(self.photons),int(self.Npsf+self.Nback),int(self.Npsf),int(np.sqrt(self.Npsfe)),int(self.Nback),int(np.sqrt(self.Nbacke))))
        print ('Called likelihood %d times'%cm.calls)
        self.cm=cm

    ## Finds boresight alignment solution
    def solvefish(self,freepsf=False):

        #estimate Isotropic component
        self.solveback()

        psfa = PSFFish(lims=[self.minroi/rd,self.maxroi/rd],free=[freepsf,freepsf,True],ebar=self.ebar,ec=self.cls,irf=self.irf)
        #print (psfa.limits)
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm = CompositeModel()
        cm.addModel(psfa)
        cm.addModel(bck)
        cm.fit(self.photons,free=[True,True],exp=[self.Npsf,self.Nback],mode=1,quiet=self.quiet)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Nbacke = np.sqrt(self.errs[1][1])
        self.params=[x*rd*3600 for x in cm.minuit.params[4:]]
        self.errors=[np.sqrt(self.errs[x+2][x+2])*rd*3600 for x in range(1)]
        self.il=cm.extlikelihood([self.Npsf,self.Nback,psfa.model_par[0],psfa.model_par[1],0],self.photons)
        for x in range(1):
            print ('dTh %3.0f +/-%3.0f arcsec'%(self.params[x],self.errors[x]))
        try:
            TS = -2*(fl-self.il)
            print ('Significance of rotation was TS = %d'%TS)
        except:
            print ('Cannot determine improvement')
        print ('Expected %d photons, got %d ( %d (%d) + %d (%d) )'%(len(self.photons),int(self.Npsf+self.Nback),int(self.Npsf),int(np.sqrt(self.Npsfe)),int(self.Nback),int(np.sqrt(self.Nbacke))))
        print ('Called likelihood %d times'%cm.calls)
        self.cm=cm

    ## tries to solve parameters for a single King function with Isotropic background
    #  @param func optional custom background component to fit [[angular separations],[counts]]
    #  @param free free or fix model number estimators
    #  @param exp specify number estimators
    #  @param theta fit theta dependence
    def solvepsf(self,func=[[],[]],free=[True,True],exp=[],theta=False,outfile=None,model_par=[]):

        #estimate Isotropic component
        self.solveback(model_par=model_par)

        cpars = cp.copy(model_par)
        #solve for PSF parameters
        psf = self.getpsf(True,theta,False)
        
        #print (model_par)
        if not model_par==[]:
            psf.model_par=model_par
            #print (psf.model_par)


        sigma = psf.model_par[0]
        gamma = psf.model_par[1]
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)

        #if custom background is specified, add it
        if len(func[0])!=0:
            bck2 = Custom(lims=[self.minroi/rd,self.maxroi/rd],func=func)
            cm.addModel(bck2)

        #check to see if number estimators are specified
        if exp==[]: 
            if len(cm.models)==2:
                exp=[self.Npsf,self.Nback]
            else:
                exp=[self.Npsf,0.,self.Nback]
                free=[True,True,True]

        #check for nearby point sources, add to model
        if self.offsrc!=[]:
            for sr in self.offsrc:
                exp.append(1.)
                free.append(True)
                off = OffPsf(lims=bck.lims,model_par=psf.model_par,off=sr)
                for sr2 in self.srcs:
                    diff = sr.difference(sr2)*rd
                    if diff<self.maxroi:
                        off.addsrc(sr2)
                cm.addModel(off)

        #fit all free parameters
        cm.fit(self.photons,free=free,exp=exp,mode=1,quiet=self.quiet)

        #get parameters and errors
        nm = len(cm.models)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Nbacke = np.sqrt(self.errs[1][1])
        if nm==3:
            self.Ncust = cm.minuit.params[2]
            self.Ncuste = np.sqrt(self.errs[2][2])
        self.sigma = cm.minuit.params[nm]
        self.sigmae = np.sqrt(self.errs[nm][nm])
        self.gamma = cm.minuit.params[nm+1]
        self.gammae = np.sqrt(self.errs[nm+1][nm+1])
        if len(psf.model_par)>2:
            self.sigma2 = cm.minuit.params[nm+2]
            self.sigma2e = np.sqrt(self.errs[nm+2][nm+2])
            self.gamma2 = cm.minuit.params[nm+3]
            self.gamma2e = np.sqrt(self.errs[nm+3][nm+3])
            self.frac = cm.minuit.params[nm+4]
            self.frace = np.sqrt(self.errs[nm+4][nm+4])
        if theta:
            self.ms = cm.minuit.params[nm+2]
            self.mse = np.sqrt(self.errs[nm+2][nm+2])
            self.mg = cm.minuit.params[nm+3]
            self.mge = np.sqrt(self.errs[nm+3][nm+3])
        self.cov = self.errs[nm][nm+1]
        if len(cm.models)==2:
            if theta:
                self.il=cm.extlikelihood([self.Npsf,self.Nback,self.sigma,self.gamma,0.,0.],self.photons)
            else:
                if len(psf.model_par)>2:
                    self.il=cm.extlikelihood([self.Npsf,self.Nback,self.sigma,self.gamma,self.sigma2,self.gamma2,self.frac],self.photons)
                else:
                    self.il=cm.extlikelihood([self.Npsf,self.Nback,self.sigma,self.gamma],self.photons)
        else:
            if theta:
                self.il=cm.extlikelihood([self.Npsf,self.Nback,0.,sigma,gamma,0.,0.],self.photons)
            else:
                self.il=cm.extlikelihood([self.Npsf,self.Nback,0.,sigma,gamma],self.photons)
        self.r68 = cm.models[0].rcl(0.68)
        self.r95 = cm.models[0].rcl(0.95)
        #if len(psf.model_par)>2:
        #    self.r68e = self.r68*np.sqrt((self.sigmae/self.sigma)**2+(self.sigma2e/self.sigma2)**2)#cm.models[0].clerr(0.68,self.sigmae,self.gammae,self.cov)
        #    self.r95e = self.r95*np.sqrt((self.sigmae/self.sigma)**2+(self.sigma2e/self.sigma2)**2)#cm.models[0].clerr(0.95,self.sigmae,self.gammae,self.cov)
        #else:
        self.r68e = cm.models[0].clerr(0.68,self.sigmae,self.gammae,self.cov)
        self.r95e = cm.models[0].clerr(0.95,self.sigmae,self.gammae,self.cov)
        if (self.r68+self.r68e)>self.r95:
            self.r68e = self.r95-self.r68
        if (self.r95-self.r95e)<self.r68:
            self.r95e = self.r95-self.r68
        r68o = cm.models[0].recl(0.68,sigma,gamma)
        r95o = cm.models[0].recl(0.95,sigma,gamma)

        phs = (self.Npsf+self.Nback)
        tr = sum([self.errs[i][i] for i in range(len(self.errs))])
        f1 = self.Npsf/phs
        f1e = f1*np.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        f2 = self.Nback/phs
        f2e = f2*np.sqrt(self.errs[1][1]/(self.Nback**2)+tr/(phs**2))
        print ('**********************************************************')
        print ('Npsf  = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100))
        print ('Nback = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f2*100,f2e*100))
        print ('      Dens  = %1.0f/deg**2'%(self.Nback/(self.maxroi**2-self.minroi**2)))
        print ('Sigma = %1.3f [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.sigma*rd,self.sigmae/self.sigma,self.irf,self.sigma/sigma))
        print ('Gamma = %1.2f  [1 +/- %1.2f]        Ratio to %s: (%1.2f)'%(self.gamma,self.gammae/self.gamma,self.irf,self.gamma/gamma))
        print ('R68   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r68*rd,self.r68e/self.r68,self.irf,self.r68/r68o))
        print ('R95   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r95*rd,self.r95e/self.r95,self.irf,self.r95/r95o))
        if theta:
            print ('ms = %1.4f [1 +/- %1.4f]'%(self.ms,self.mse/abs(cm.minuit.params[nm+2])))
            print ('mg = %1.4f [1 +/- %1.4f]'%(self.mg,self.mge/abs(cm.minuit.params[nm+3])))
        self.Ts = -2*(fl-self.il)
        print ('Significance of psf change was TS = %d'%self.Ts)
        if outfile is not None:
            print ('**********************************************************', file=outfile)
            print ('Npsf  = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100), file=outfile)
            print ('Nback = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f2*100,f2e*100), file=outfile)
            print ('      Dens  = %1.0f/deg**2'%(self.Nback/(self.maxroi**2-self.minroi**2)), file=outfile)
            print ('Sigma = %1.4f [1 +/- %1.4f] (deg)  Ratio to %s: (%1.2f)'%(self.sigma*rd,self.sigmae/self.sigma,self.irf,self.sigma/sigma), file=outfile)
            print ('Gamma = %1.4f  [1 +/- %1.4f]        Ratio to %s: (%1.2f)'%(self.gamma,self.gammae/self.gamma,self.irf,self.gamma/gamma), file=outfile)
            print ('R68   = %1.4f  [1 +/- %1.4f] (deg)  Ratio to %s: (%1.2f)'%(self.r68*rd,self.r68e/self.r68,self.irf,self.r68/r68o), file=outfile)
            print ('R95   = %1.4f  [1 +/- %1.4f] (deg)  Ratio to %s: (%1.2f)'%(self.r95*rd,self.r95e/self.r95,self.irf,self.r95/r95o), file=outfile)
            if theta:
                print ('ms = %1.4f [1 +/- %1.4f]'%(cm.minuit.params[nm+2],np.sqrt(self.errs[nm+2][nm+2])/abs(cm.minuit.params[nm+2])), file=outfile)
                print ('mg = %1.4f [1 +/- %1.4f]'%(cm.minuit.params[nm+3],np.sqrt(self.errs[nm+3][nm+3])/abs(cm.minuit.params[nm+3])), file=outfile)
            self.Ts = -2*(fl-self.il)
            print ('Significance of psf change was TS = %d'%self.Ts, file=outfile)
        self.cm=cm
        if theta:
            for ct in np.linspace(0.,1.,24):
                tsig,tgam = ((1-ct)*self.ms+1)*self.sigma*rd,((1-ct)*self.mg+1)*self.gamma
                r68,r95 = cm.models[0].recl(0.68,tsig,tgam),cm.models[0].recl(0.95,tsig,tgam)
                print ('%1.2f %1.3f %1.2f %1.2f %1.2f'%(ct,tsig,tgam,r68,r95))

    # tries to fit a halo component on top of a PSF defined by 'irf' in a Isotropic background
    #  @param cust [sigma,gamma] of PSF, where sigma is in radians
    def solvehalo(self,cust=[],opt=True):

        #estimate Isotropic component
        self.solvepsf()

        #solve for Halo model component
        if cust==[]:
            psf = self.getpsf(opt)
        else:
            psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=cust,free=[opt,opt])
        halo = Halo(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.2/rd])#PLCutoff(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.2/rd,-2.,self.ebar,self.cls])
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])

        #setup model
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(halo)
        cm.addModel(bck)


        #fit all free parameters
        cm.fit(self.photons,free=[True,True,True],exp=[self.Npsf/2.,self.Npsf/2.,self.Nback],mode=1,quiet=self.quiet)

        #get parameters and errors
        self.errs = cm.minuit.errors()
        self.Npsf = cm.minuit.params[0]
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Nhalo = cm.minuit.params[1]
        self.Nhaloe = np.sqrt(self.errs[1][1])
        self.Nback = cm.minuit.params[2]
        self.Nbacke = np.sqrt(self.errs[2][2])
        self.theta = cm.minuit.params[5]
        self.thetae = np.sqrt(self.errs[5][5])/self.theta
        self.alph = cm.minuit.params[6]
        self.alphe = np.sqrt(self.errs[6][6])/self.alph
        self.frac = self.Nhalo*100./(self.Nhalo+self.Npsf)
        ferr = self.frac*np.sqrt(self.Nhaloe**2/(self.Nhalo**2)+self.Npsfe**2/(self.Npsf**2))
        self.frace = ferr/self.frac
        phs = (self.Npsf+self.Nback+self.Nhalo)
        tr = sum([self.errs[i][i] for i in range(len(self.errs))])
        f1 = self.Npsf/phs
        f1e = f1*np.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        f2 = self.Nhalo/phs
        f2e = f2*np.sqrt(self.errs[1][1]/(self.Nhalo**2)+tr/(phs**2))
        f3 = self.Nback/phs
        f3e = f3*np.sqrt(self.errs[2][2]/(self.Nback**2)+tr/(phs**2))
        print ('**********************************************************')
        print ('Npsf  = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100))
        print ('Nhalo = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Nhalo,self.Nhaloe/self.Nhalo,f2*100,f2e*100))
        print ('Nback = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f3*100,f3e*100))
        print ('Halo width was %1.3f [1 +/- %1.2f] deg'%(self.theta*rd,self.thetae))
        #print ('Halo PL was %1.3f [1 +/- %1.2f] deg'%(self.alph,self.alphe))
        print ('Halo fraction was %1.0f [1 +/- %1.2f]'%(self.frac,self.frace))

        fl = cm.minuit.fval
        null = cm.extlikelihood([self.Npsf,0.,self.Nback,self.sigma,self.gamma,self.theta],self.photons)
        print ('Significance of Halo was: %1.0f'%(-2*(fl-null)))
        self.cm=cm

    ## tries to fit two King functions in an Isotropic background
    #  @param func optional custom background component to fit [[angular separations],[counts]]
    #  @param free free or fix model number estimators
    #  @param exp specify number estimators
    def solvedoublepsf(self,func=[[],[]],free=[True,True],exp=[]):
        
        #estimate uniform Isotropic component
        self.solveback()

        #solve two King function parameters
        psf = self.getpsf(True)
        sigma = psf.model_par[0]
        gamma = psf.model_par[1]
        psf2 = PSFDouble(lims=[self.minroi/rd,self.maxroi/rd],model_par=[psf.model_par[0],psf.model_par[1],1.2*psf.model_par[0],0.7*psf.model_par[1],0.5]) #try to separate psfs initially
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])

        #setup models
        cm = CompositeModel()
        #cm.addModel(psf)
        cm.addModel(psf2)
        cm.addModel(bck)

        #check for custom background
        if len(func[0])!=0:
            bck2 = Custom(lims=[self.minroi/rd,self.maxroi/rd],func=func)
            cm.addModel(bck2)

        #check for specified estimators
        if exp==[]: 
            if cm.models==3:
                exp=[self.Npsf,self.Nback]
                free = [True,True]
            else:
                exp=[self.Npsf,self.Nback,0.]
                free=[True,True,True]

        #check for nearby sources, add to model
        if self.offsrc!=[]:
            for sr in self.offsrc:
                exp.append(1.)
                free.append(True)
                off = OffPsf(lims=bck.lims,model_par=psf.model_par,off=sr)
                for sr2 in self.srcs:
                    diff = sr.difference(sr2)*rd
                    if diff<self.maxroi:
                        off.addsrc(sr2)
                cm.addModel(off)

        #maximize likelihood
        cm.fit(self.photons,free=free,exp=exp,mode=1,quiet=self.quiet)

        #get parameters and errors
        nm = len(cm.models)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        #self.Npsf2 = cm.minuit.params[1]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = np.sqrt(self.errs[0][0])
        #self.Npsf2e = np.sqrt(self.errs[1][1])
        self.Nbacke = np.sqrt(self.errs[1][1])
        if nm==4:
            self.Ncust = cm.minuit.params[2]
            self.Ncuste = np.sqrt(self.errs[2][2])
        self.sigma = cm.minuit.params[nm]
        self.sigmae = np.sqrt(self.errs[nm][nm])
        self.gamma = cm.minuit.params[nm+1]
        self.gammae = np.sqrt(self.errs[nm+1][nm+1])
        self.sigma2 = cm.minuit.params[nm+2]
        self.sigmae2 = np.sqrt(self.errs[nm+2][nm+2])
        self.gamma2 = cm.minuit.params[nm+3]
        self.gammae2 = np.sqrt(self.errs[nm+3][nm+3])
        self.frac = cm.minuit.params[nm+4]
        self.frace = np.sqrt(self.errs[nm+4][nm+4])
        if nm==3:
            self.il=cm.extlikelihood([self.Npsf,self.Nback,sigma,gamma,sigma,gamma,0.5],self.photons)
        else:
            self.il=cm.extlikelihood([self.Npsf,self.Nback,0.,sigma,gamma,sigma,gamma,0.5],self.photons)

        if self.gamma>gamma:
            gsign = '+'
        else:
            gsign = '-'
        if self.sigma>sigma:
            ssign = '+'
        else:
            ssign = '-'


        if self.gamma2>gamma:
            gsign2 = '+'
        else:
            gsign2 = '-'
        if self.sigma2>sigma:
            ssign2 = '+'
        else:
            ssign2 = '-'
        phs = (self.Npsf+self.Nback)#+self.Npsf2)
        tr = sum([self.errs[i][i] for i in range(len(self.errs))])
        f1 = self.Npsf/phs
        f1e = f1*np.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        #f2 = self.Npsf2/phs
        #f2e = f2*np.sqrt(self.errs[1][1]/(self.Npsf2**2)+tr/(phs**2))
        f3 = self.Nback/phs
        f3e = f3*np.sqrt(self.errs[1][1]/(self.Nback**2)+tr/(phs**2))

        def integrate(d1,d2,s,g):
            return ((1+(d1/s)**2/(2*g))**(-g+1)-(1+(d2/s)**2/(2*g))**(-g+1))*(s**2)

        def cintegrate(d1,d2,s1,g1,n1,s2,g2,n2):
            return n1*integrate(d1,d2,s1,g1)+n2*integrate(d1,d2,s2,g2)

        def findcontain(frac,s1,g1,n1,s2,g2,n2):
            fint = cintegrate(0,1e40,s1,g1,n1,s2,g2,n2)
            best = so.fmin_powell(lambda x: (cintegrate(0,x[0],s1,g1,n1,s2,g2,n2)/fint-frac)**2,[s1],full_output=1,disp=1)
            return best[0]

        self.r68 = findcontain(0.68,self.sigma*rd,self.gamma,self.Npsf,self.sigma2*rd,self.gamma2,(1-self.frac)*self.Npsf)
        self.r95 = findcontain(0.95,self.sigma*rd,self.gamma,self.Npsf,self.sigma2*rd,self.gamma2,(1-self.frac)*self.Npsf)

        print ('***************************************')
        print ('Npsf   = %1.0f  [1 +/- %1.2f]     Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100))
        #print ('Npsf2  = %1.0f  [1 +/- %1.2f]     Fraction: %1.0f +/- %1.0f'%(self.Npsf2,self.Npsf2e/self.Npsf2,f2*100,f2e*100))
        print ('Nback  = %1.0f  [1 +/- %1.2f]     Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f3*100,f3e*100))
        print ('Sigma  = %1.3f [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.sigma*rd,self.sigmae/self.sigma,ssign,100*abs((self.sigma-sigma)/sigma)))
        print ('Gamma  = %1.2f  [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.gamma,self.gammae/self.gamma,gsign,100*abs((self.gamma-gamma)/gamma)))
        print ('Sigma2 = %1.3f [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.sigma2*rd,self.sigmae2/self.sigma2,ssign2,100*abs((self.sigma2-sigma)/sigma)))
        print ('Gamma2 = %1.2f  [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.gamma2,self.gammae2/self.gamma2,gsign2,100*abs((self.gamma2-gamma)/gamma)))
        #print ('frac   = %1.3f  [])
        print ('R68 = %1.3f'%(self.r68))
        print ('R95 = %1.3f'%(self.r95))
        TS = -2*(fl-self.il)
        print ('Significance of psf change was TS = %d'%TS)
        self.cm=cm

    ## solvediffuse - use the diffuse background information
    #  @param fit optimize PSF
    #  @param exp2 number estimators
    #  @param free2 for offset point sources, free or fixed number estimators
    #  @param ptfree free or fix the central PSF number estimator (only estimate diffuse/isotropic)
    #  @param cust custom [sigma,gamma] for PSF
    def solvediffuse(self,fit=False,exp2=[],free2=[],ptfree=True,cust=[]):
        #Try to estimate the number of Isotropic photons
        #by looking at density near tails
        ct = ['front','back']
        self.getds()
        frac = 0.98
        area = 1-frac*frac
        density = len(self.ds[self.ds>(frac*self.maxroi/rd)])
        denback = density*(self.maxroi*self.maxroi-self.minroi*self.minroi)/((self.maxroi**2)*area)
        if denback>len(self.ds):
            denback = len(self.ds)
        exp=[len(self.ds)-denback,denback/2.,denback/2.]

        diffusefile = basedir+'data/galprop/ring_21month_P6v11.fits'
        s.EffectiveArea.set_CALDB(basedir+'packages/ScienceTools-09-19-00/irfs/caldb/CALDB/data/glast/lat')
        efa = 'P6_v11_diff_%s'%(ct[self.cls])
        ltc = basedir+'data/2yr/2years_lt.fits'

        #setup models: PSF, diffuse, and isotropic
        if cust==[]:
            psf = self.getpsf(opt=fit)
        else:
            psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=cust,free=[fit,fit])
        sigma = psf.model_par[0]
        gamma = psf.model_par[1]
        bck = Diffuse(lims=[self.minroi/rd,self.maxroi/rd],diff=diffusefile,ebar=self.ebar,ltc=ltc,ea=efa)
        bck2 = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        offs = []

        #check for nearby sources
        for off in self.offsrc:
            toff = OffPsf(lims=[self.minroi/rd,self.maxroi/rd],model_par=[psf.model_par[0],psf.model_par[1]],off=off)
            offs.append(toff)

        #add sources to background
        for src in self.srcs:
            bck.addsrc(src)
            for tof in offs:
                tof.addsrc(src)
        
        #setup models
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)
        cm.addModel(bck2)
        for tof in offs:
            cm.addModel(tof)
        nm = len(cm.models)

        #METHOD1 - no nearby point sources
        if nm<=3:
            print (exp)

            #maximize likelihood
            cm.fit(self.photons,mode=0,quiet=self.quiet,free=[True,True,True],exp=exp)

            #get parameters and errors
            self.Npsf = cm.minuit.params[0]
            self.Nbacki = cm.minuit.params[1]
            self.Nbackd = cm.minuit.params[2]
            self.Nback = self.Nbacki+self.Nbackd
            self.errs = cm.minuit.errors()
            self.Npsfe = np.sqrt(self.errs[0][0])
            self.Nbacke = np.sqrt(self.errs[1][1]+self.errs[2][2])
            self.cm = cm
            self.sigma = cm.minuit.params[nm]
            self.sigmae = np.sqrt(self.errs[nm][nm])
            self.gamma = cm.minuit.params[nm+1]
            self.gammae = np.sqrt(self.errs[nm+1][nm+1])
            self.cov = self.errs[nm][nm+1]
            fl = cm.minuit.fval
            self.il=cm.extlikelihood([self.Npsf,self.Nbackd,self.Nbacki,sigma,gamma],self.photons)
            self.r68 = cm.models[0].rcl(0.68)
            self.r95 = cm.models[0].rcl(0.95)
            self.r68e = cm.models[0].clerr(0.68,self.sigmae,self.gammae,self.cov)
            self.r95e = cm.models[0].clerr(0.95,self.sigmae,self.gammae,self.cov)
            if (self.r68+self.r68e)>self.r95:
                self.r68e = self.r95-self.r68
            if (self.r95-self.r95e)<self.r68:
                self.r95e = self.r95-self.r68
            r68o = cm.models[0].recl(0.68,sigma,gamma)
            r95o = cm.models[0].recl(0.95,sigma,gamma)

            phs = (self.Npsf+self.Nback)
            tr = sum([self.errs[i][i] for i in range(len(self.errs))])
            f1 = self.Npsf/phs
            f1e = f1*np.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
            f2 = self.Nback/phs
            f2e = f2*np.sqrt(self.errs[1][1]/(self.Nback**2)+tr/(phs**2))
            print ('**********************************************************')
            print ('Npsf  = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100))
            print ('Nback = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f2*100,f2e*100))
            print ('Sigma = %1.3f [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.sigma*rd,self.sigmae/self.sigma,self.irf,self.sigma/sigma))
            print ('Gamma = %1.2f  [1 +/- %1.2f]        Ratio to %s: (%1.2f)'%(self.gamma,self.gammae/self.gamma,self.irf,self.gamma/gamma))
            print ('R68   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r68*rd,self.r68e/self.r68,self.irf,self.r68/r68o))
            print ('R95   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r95*rd,self.r95e/self.r95,self.irf,self.r95/r95o))
            TS = -2*(fl-self.il)
            print ('Significance of psf change was TS = %d'%TS)

        #METHOD2 - nearby point sources
        else:
            left = len(self.ds)-denback

            #get number estimators for nearby sources
            #if exp == []:
            #    exp = [left/(nm-2.) for x in range(nm)]
            #    free = [True for x in range(nm)]
            #    exp[1]=denback/2.
            #    exp[2]=denback/2.
            #    cm.fit(self.photons,mode=0,quiet=self.quiet,free=free,exp=exp)
            
            #number estimators are specified
            #else:
            #tfile = open('/phys/users/mar0/figures/%s_%d_%d.txt'%(self.lis,self.ebar,self.cls),'w')
            
            #calculate photons left for offset sources
            npsfs = sum(exp2)
            left = len(self.ds)-npsfs
            exp2.insert(1,left)
            exp2.insert(2,0)
            texp = cp.copy(exp2)

            #setup free estimators
            free2.insert(0,ptfree)
            free2.insert(1,True)
            free2.insert(2,True)
            print (free2,exp2)

            #maximize likelihood
            cm.fit(self.photons,mode=0,quiet=self.quiet,free=free2,exp=exp2)
            exp2 = cm.minuit.params


            #print ('pred\tfit\tdiff',file=tfile)
            #for it5 in range(len(texp)):
            #    print ('%1.4f\t%1.4f\t%1.4f'%(texp[it5],exp2[it5],(exp2[it5]-texp[it5])), file=tfile)
            #tfile.close()


            #get parameters
            print (cm.minuit.params)
            errs = cm.minuit.errors()
            pars = len(cm.minuit.params)
            sigma, gamma = cm.minuit.params[pars-2],cm.minuit.params[pars-1]
            sige, game, cov = np.sqrt(errs[pars-2][pars-2]),np.sqrt(errs[pars-1][pars-1]),errs[pars-2][pars-1]
            r68n,r95n = cm.models[0].rcl(0.68),cm.models[0].rcl(0.95)
            r68e,r95e = cm.models[0].clerr(0.68,sige,game,cov),cm.models[0].clerr(0.95,sige,game,cov)
            print (r68n*rd,r68e*rd)
            print (r95n*rd,r95e*rd)


            #for i in range(len(cm.minuit.params)):
            #    outline = ''
            #    for j in range(len(cm.minuit.params)):
            #        outline = outline + '%1.6f\t'%errs[i][j]
            #    print (outline)
            #self.cm = cm


    def solvehalolims(self,frac):
        rf = uf.ROIfactory(binfile=basedir+'data/2yr/2years_4bpd.fits',ltcube=basedir+'data/2yr/2years_lt.fits',irf='P6_v11_diff')
        roi = rf('vela')
        for band in roi.bands:
            if band.ec==self.cls and self.ebar>band.emin and self.ebar<band.emax:
                eb = band
        self.solveback()
        psf = self.getpsf(False)
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        #cm = CompositeModel()
        #cm.addModel(psf)
        #cm.addModel(bck)
        #cm.fit(self.photons,free=[True,True],exp=[self.Npsf,self.Nback],quiet=self.quiet,mode=1)
        psf2 = psf#PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=psf.model_par,free=[False,False])
        initial = self.cm.minuit.fval

        uplim=[]
        thetas = np.arange(-16,6,1)
        thetas = thetas/8.
        thetas = 10**(thetas)
        py.ioff()
        py.figure(figsize=(8,8))
        uthetas = []
        for theta in thetas:
            if theta>psf2.model_par[0]*rd and theta<self.maxroi:
                halo = CHalo(lims=[self.minroi/rd,self.maxroi/rd],model_par=[theta/rd],free=[False],energy=self.ebar,ct=self.cls)
                cm = CompositeModel()
                cm.addModel(psf2)
                cm.addModel(bck)
                cm.addModel(halo)
                cm.fit(self.photons,free=[True,True,True],exp=[self.Npsf,self.Nback,0],quiet=self.quiet,mode=0)
                #self.il=cm.minuit.fval
                #self.errors=cm.minuit.errors()
                self.Nhalo=cm.minuit.params[2]
                def like(x):
                    if x<0 or x>self.Npsf:
                        return 0.
                    val = np.exp(cm.minuit.fval-cm.extlikelihood([self.Npsf-x,self.Nback,x,psf2.model_par[0],psf2.model_par[1]],self.photons,True))
                    #print (val)
                    return val
                
                def integ(x):
                    val = si.quad(lambda y:like(y),0,x)[0]
                    return val

                def contain(x,fnt,ct,limit=False):
                    ft = integ(x)
                    if ft/fnt>ct:
                        return 1.
                    val = abs(ft/fnt-ct)
                    print (x,val)
                    return val

                fint = integ(self.Npsf)
                
                best = so.fmin_powell(lambda x:contain(x,fint,frac),self.Nhalo,maxiter=2,full_output=1,disp=1)
                upl = best[0]#self.Nhalo+np.sqrt(self.errors[2][2])*2
                uplim.append(upl/(self.Npsf))
                uthetas.append(theta)
                print (upl)
                cm.nest=[self.Npsf-best[0],self.Nback,best[0]]
                self.cm=cm
                self.makeplot('/phys/users/mar0/figures/uplims_%1.3f_%1.0f_%s_%s_%s.png'%(theta,self.ebar,self.lis,self.irf,halo.name),fig=1)
        return uthetas,uplim

    ## Makes a plot of all of the models and the angular distribution of photons
    #  @param name filename of output file
    #  @param fig fignum of pylab figure, default is 1
    #  @param bin number of bins for unbinned data
    #  @param scale makes the x-scale log/linear  ['log','linear']
    #  @param jac either counts/deg or counts/(deg**2) ['square','linear']
    #  @param xs x-scale either deg or deg**2 ['linear','square']
    #  @param outhist creates an output histogram of all models and binned data
    #  @param stats put statistics on right hand side of the plot
    def makeplot(self,name='test',fig=-1,bin=25.,scale='log',jac='square',xs='linear',xlims=[],outhist=False,stats=False):
        #calculate angular separations
        self.getds()

        if not xlims==[]:
            self.ds = self.ds[(self.ds>(xlims[0]/rd))&(self.ds<(xlims[1]/rd))]

        #plotting setup
        bins = bin                       #angular bins

        #figure out scales depending on scale/jacobian parameters
        if scale == 'log':
            if xs=='square':
                jmin = min(self.ds)**2
                jmax = max(self.ds)**2
            else:
                jmin = min(self.ds) 
                jmax = max(self.ds)
            mi = np.log10(jmin)
            ma = np.log10(jmax)
        else:
            if xs=='square':
                jmin = (self.minroi/rd)**2#min(self.ds)**2
                jmax = (self.maxroi/rd)**2#max(self.ds)**2
            else:
                jmin = self.minroi/rd#min(self.ds) 
                jmax = self.maxroi/rd#max(self.ds)
            mi = jmin
            ma = jmax
        d2 = ma-mi
        be = np.arange(mi,ma,d2/bins)
        if be[len(be)-1]-ma<0:
            be2 = np.append(be,ma)#be[len(be)-1]+d2/bins)
        else:
            be2 = np.append(be,ma+d2/bins)
        if scale =='log':
            if xs=='square':
                be3 = np.sqrt(10**(be2))
                be4 = (10**(be+d2/bins/2))*rd*rd
            else:
                be3 = 10**(be2)
                be4 = (10**(be+d2/bins/2))*rd
        else:
            if xs=='square':
                be3 = np.sqrt(be2)
                be4 = (be+d2/bins/2)*rd*rd
            else:
                be3 = be2
                be4 = (be+d2/bins/2)*rd
        py.ioff()
        if fig ==-1:
            py.figure(figsize=(8,8))
        else:
            py.clf()
        if scale == 'log':
            if xs=='square':
                binning = np.log10(self.ds*self.ds)
            else:
                binning = np.log10(self.ds)
        else:
            if xs=='square':
                binning = self.ds*self.ds
            else:
                binning = self.ds
        hist = py.hist(binning,bins=be2,fc='None')
        py.clf()
        fits=[]
        hists=[]
        errs=[]
        names=[]                      #names of the various plots
        ctn=[]
        names.append('data')
        
        setoff = True
        #retrieve all model names
        for model in range(len(self.cm.models)):
            fits.append([])
            if self.cm.models[model].name=='offpsf':
                if setoff:
                    names.append('offpsf')
                    setoff=False
            else:
                names.append(self.cm.models[model].name)
            if self.cm.models[model].name=='psf' or self.cm.models[model].name=='psf2':
                names.append('psf r68')
                names.append('psf r95')
                if xs=='square':
                    ctn.append((self.cm.models[model].rcl(0.68)*rd)**2)
                    ctn.append((self.cm.models[model].rcl(0.95)*rd)**2)
                else:
                    ctn.append(self.cm.models[model].rcl(0.68)*rd)
                    ctn.append(self.cm.models[model].rcl(0.95)*rd)


        fits.append([])
        names.append('total')

        #fill all model plots
        for it,ba in enumerate(be):
            left=be3[it]
            right=be3[it+1]
            if jac=='square':
                area = (right**2-left**2)/(rd**2)
            else:
                area = (right-left)/rd
            for it2,model in enumerate(self.cm.models):
                fits[it2].append(self.cm.nest[it2]*model.integral(left,right)/model.integral(self.minroi/rd,self.maxroi/rd)/area)
            fits[len(fits)-1].append(self.cm.integral(left,right,self.minroi/rd,self.maxroi/rd)/area)
            hists.append(hist[0][it]/area)
            errs.append(max(min(np.sqrt(hist[0][it])/area,(hist[0][it]/area)*(1.-1.e-15)),1e-40))

        pts=[]
        pmin,pmax=0.01*min(fits[len(fits)-1]),100*max(fits[len(fits)-1])
        if stats:
            py.subplot(2,2,1)
        else:
            py.subplot(2,1,1)
        #plot histogrammed events
        py.errorbar(be4,hists,yerr=errs,ls='None',marker='None',ecolor='blue',capsize=80./bin)
        p1 = py.step(be4,hists,where='mid')
        pts.append(p1[0])

        pcts=0

        setoff = True
        #plot all models
        for model in range(len(self.cm.models)):
            p1 = py.plot(be4,fits[model],self.cm.models[model].mark,label=self.cm.models[model].name)
            if self.cm.models[model].name=='offpsf':
                if setoff:
                    pts.append(p1)
                    setoff=False
            else:
                pts.append(p1)
            if self.cm.models[model].name=='psf' or self.cm.models[model].name=='psf2':
                p1 = py.plot([ctn[2*pcts],ctn[2*pcts]],[pmin,pmax],'-.',label=self.cm.models[model].name+'68')
                pts.append(p1)
                p1 = py.plot([ctn[2*pcts+1],ctn[2*pcts+1]],[pmin,pmax],'-.',label=self.cm.models[model].name+'95')
                pts.append(p1)
                pcts=pcts+1
        p1 = py.plot(be4,fits[len(fits)-1],self.cm.mark,label='total')
        pts.append(p1)

        #calculate chisq
        self.chisq=0.
        total = fits[len(fits)-1]
        for it,bin in enumerate(hists):
            if errs[it]>1e-35:
                self.chisq = self.chisq + ((bin-total[it])/errs[it])**2

        prop = matplotlib.font_manager.FontProperties(size=9) 
        #finish up plotting
        py.legend(bbox_to_anchor=(1.1, 1.0),prop=prop)
        if scale=='log':
            py.loglog()
        else:
            py.semilogy()
        py.ylim(pmin,pmax)
        if scale=='log':
            py.xlim(0.8*min(be4),max(be4)*1.2)
        else:
            py.xlim(0,max(be4)*1.01)
        if xs=='square':
            py.xlabel(r'$\theta^2\/(\rm{deg^2})$')
        else:
            py.xlabel(r'$\theta\/(\rm{deg})$')
        if jac=='square':
            py.ylabel(r'$dN/d\theta^2$')
        else:
            py.ylabel(r'$dN/d\theta$')
        py.grid()
        clas = ['front','back']
        py.suptitle('%s: Energy(MeV): [%1.0f,%1.0f]\nClass: %s    Costh: [%1.1f,%1.1f]\nChisq: %1.1f / %1.0f dof'%(self.name,self.emin,self.emax,clas[self.cls],self.ctmin,self.ctmax,self.chisq,len(be2)-1))
        if stats:
            py.subplot(2,2,3)
        else:
            py.subplot(2,1,2)
        if scale=='log':
            py.semilogx()
        py.plot([0.8*min(be4),max(be4)*1.2],[0,0],'r-')
        py.errorbar(be4,(np.array(hists)-np.array(total))/np.array(errs),yerr=1,marker='o',ls='None')
        if scale=='log':
            py.xlim(0.8*min(be4),max(be4)*1.2)
        else:
            py.xlim(0,max(be4)*1.01)
        py.ylim(-5,5)
        py.grid()
        if xs=='square':
            py.xlabel(r'$\theta^2\/(\rm{deg^2})$')
        else:
            py.xlabel(r'$\theta\/(\rm{deg})$')
        py.ylabel(r'$\frac{\rm{Data}-\rm{Model}}{\sigma}$')
        head = ''
        for md in self.cm.models:
            head = head+'N'+md.name+'\n'
        for md in self.cm.models:
            head = head+md.header.replace('\t','\n')
        data = ''
        terrs = self.errs
        for it,pm in enumerate(self.cm.minuit.params):
            if it<len(self.cm.models):
                data = data + '%1.0f (%1.0f)\n'%(pm,np.sqrt(terrs[it][it]))
            else:
                data = data +'%1.6f (%1.6f)\n'%(pm,np.sqrt(terrs[it][it]))
        if stats:
            py.figtext(0.55,0.5,head)
            py.figtext(0.68,0.5,data)
        py.savefig(name+'.png')

        if outhist:
            #dt = np.dtype([('name',np.str_),('values',np.float32,(2,))])
            output = []
            headertb = []
            output.append(be4)
            headertb.append('xbins')
            output.append(hists)
            headertb.append('density')
            output.append(errs)
            headertb.append('densunc')
            for itfit,fit in enumerate(fits):
                if itfit<(len(fits)-1):
                    output.append(fit)
                    headertb.append(self.cm.models[itfit].name)
                else:
                    output.append(fit)
                    headertb.append('TOTAL')
            output.append(1.*np.array(hist[0]))
            headertb.append('hist')
            #print (output)
            #outa = np.array(output,dtype=dt)
            np.save(name+'tb.npy',output)
            np.save(name+'hdr.npy',headertb)

    ## helper function - cuts unused data
    #  @param table 'EVENTS' FITS table from FT1 data
    #  @param srcs list of skymaps::SkyDirs of sources to stack
    #  @param emin minimum energy (MeV)
    #  @param emax maximum energy (MeV)
    #  @param start minimum MET
    #  @param stop maximum MET
    #  @param rad ROI in degrees around sources
    #  @param cls conversion type: 0=front,1=back,-1=all
    def mask(self,table,srcs,emin,emax,start,stop,rad,cls):
        cuts = [0,0,0,0,0,0,0]
        total = len(table)
        tc = len(table)

        #make time cuts
        msk = (table.field('TIME')>start) & (table.field('TIME')<stop)
        table = table[msk]
        tc = tc - len(table)
        cuts[0] = tc
        tc = len(table)
        if len(table)>0:
            #make energy cuts
            msk = (table.field('ENERGY')>emin) & (table.field('ENERGY')<emax)
            table = table[msk]
            tc = tc - len(table)
            cuts[1] = tc
            tc = len(table)
            if len(table)>0:
                #make instrument theta cuts
                msk = (table.field('THETA')<(np.arccos(self.ctmin)*rd)) & (table.field('THETA')>(np.arccos(self.ctmax)*rd))
                table = table[msk]
                tc = tc - len(table)
                cuts[2] = tc
                tc = len(table)
                if len(table)>0:
                    #make zenith angle cuts to remove limb photons
                    msk = table.field('ZENITH_ANGLE')<90.
                    table = table[msk]
                    tc = tc - len(table)
                    cuts[3] = tc
                    tc = len(table)
                    if len(table)>0:
                        #optionally cut out conversion type
                        if cls!=-1:
                            msk = table.field('CONVERSION_TYPE')==cls
                            table = table[msk]
                            tc = tc - len(table)
                            cuts[4] = tc
                            tc = len(table)
                        #make ROI cuts
                        if len(table)>0:
                            msk = self.dirmask(srcs[0],rad,table)
                            for j in range(len(srcs)-1):
                                msk = msk | self.dirmask(srcs[j+1],rad,table)
                            table = table[msk]
                            tc = tc - len(table)
                            cuts[5] = tc
                            tc = len(table)
                        #optionally cut on pulse phase for pulsars
                        if self.phasecut!=[] and len(table)>0:
                            msk = []
                            for it in range(len(self.phasecut)/2):
                                print (self.phasecut[2*it],self.phasecut[2*it+1])
                                tmsk = (table.field('PULSE_PHASE')>self.phasecut[2*it]) & (table.field('PULSE_PHASE')<self.phasecut[2*it+1])
                                msk = tmsk if it==0 else msk | tmsk
                            table = table[msk]
                            tc = tc - len(table)
                            cuts[6] = tc

        tc = len(table)
        #display number of photons cut from each step and finally the number of photons examined
        if not self.quiet:
            print ('TOTAL: %d    TIMEC: %d    ENERGY: %d    THETA: %d    ZENITH: %d    ECLASS: %d    POS: %d    PHASE: %d    EXAM: %d'%(total,cuts[0],cuts[1],cuts[2],cuts[3],cuts[4],cuts[5],cuts[6],tc))
        return table,np.array([total,cuts[0],cuts[1],cuts[2],cuts[3],cuts[4],cuts[5],cuts[6],tc])
    

    ## direction mask - masks area around sources to speed execution
    #  @param sd skymaps::SkyDir of source location
    #  @param deg radius in degrees
    #  @param tb 'EVENTS' FITS table from FT1 file
    def dirmask(self,sd,deg,tb):
        ra,dec=sd.ra(),sd.dec()

        #does the region overlap the poles?
        if deg+abs(dec)>90:
            if dec<0:
                ldec = dec+deg
                mask = tb.field('DEC')<ldec
            else:
                ldec = dec-deg
                mask = tb.field('DEC')>ldec
        else:
            if dec<0:
                cdec = np.cos((dec-deg)/rd)
                ramin = ra-deg/cdec
                ramax = ra+deg/cdec
                decmin = dec-deg
                decmax = dec+deg
            else:
                cdec = np.cos((dec+deg)/rd)
                ramin = ra-deg/cdec
                ramax = ra+deg/cdec
                decmin = dec-deg
                decmax = dec+deg
                   #ramin < 0 go to 360+ramin                              racut                                         deccut
            mask = ( ((tb.field('RA')>360+ramin)&(ramin<0))  |  ((tb.field('RA')>ramin)&(tb.field('RA')<ramax))  )\
            & ( (tb.field('DEC')>decmin)&(tb.field('DEC')<decmax) )
        return mask

    ##  sets scaled angular deviations for rotated events
    #   @param x1 x-rotation
    #   @param y1 y-rotation
    #   @param z1 z-rotation
    #def getus(self,x1,y1,z1):
    #    self.us=[]
    #    x,y,z=x1/rd/3600.,y1/rd/3600.,z1/rd/3600.
    #    rot = HepRotation([x,y,z],False)
    #    for photon in self.photons:
    #        sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
    #        umax = (self.maxroi/sig/rd)
    #        umax = umax*umax/2
    #        u = photon.diff(rot)/sig
    #        u = u*u/2
    #        self.us.append(u)

    #sets angular deviations
    def getds(self):
        if self.ds==[]:
            for photon in self.photons:
                if self.dsel.binfile is not None:
                    u = photon.srcdiff()
                    for k in range(int(photon.weight)):
                        self.ds.append(u)
                else:
                    if self.useft2s:
                        u = photon.diff(self.rot)
                    else:
                        u = photon.srcdiff()
                    self.ds.append(u)
            self.ds = np.array(self.ds)

    def saveds(self,fname='ds.npy'):
        if self.ds==[]:
            self.getds()
        np.save(fname,self.ds)

    def loadds(self,fname):
        self.ds = np.load(fname)
        self.srcs = [s.SkyDir(0,0)]
        for ang in self.ds:
            photon = Photon(ang*rd,0,0,0,0,[],[],self.srcs[0])
            self.photons.append(photon)
            self.photoncount = self.photoncount+1
        self.maxroi = max(self.ds)*rd
        self.minroi = min(self.ds)*rd
        self.bin=False

    def spickle(self,fname):
        import pickle
        self.getds()
        try:
            del self.photons
        except:
            pass
        try:
            del self.cm
        except:
            pass
        try:
            del self.psf
        except:
            pass
        try:
            del self.irfloader
        except:
            pass
        cfile = open(fname,'w')
        pickle.dump(self,cfile)
        cfile.close()

    ## returns a PSF object based on containment specified by self.irf
    #  @param opt frees or fixes PSF parameters
    def getpsf(self,opt=False,theta=False,double=False):
        dmask = [True,True,False,False,False]
        self.psf = pypsf.CALDBPsf(CALDBManager(irf=self.irf))
        self.irfloader = IrfLoader(self.irf)
        psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.001,2.25],free=[opt,opt])
        pars = self.irfloader.params(self.ebar,self.cls)[:-1]
        R68 = self.irfloader.rcontain(pars,0.68)
        R95 = self.irfloader.rcontain(pars,0.95)
        psf.fromcontain([R68,R95],[0.68,0.95])
        return psf
        pars = self.irfloader.params(self.ebar,self.cls)
        pars = [pars[0],pars[1],pars[3],pars[4],pars[2] if double else 1.0]
        if double:
            psf = PSFDouble(lims=[self.minroi/rd,self.maxroi/rd],model_par=pars,free=[(opt and dmask[x]) or double for x in range(len(pars))])
        if theta:
            psf = PSFTheta(lims=psf.lims,model_par=[psf.model_par[0],psf.model_par[1],0.,0.],free=[opt,opt,True,True])
        return psf

    def writebpd(self,bpd=8):
        binfile = s.BinnedPhotonData(bpd)
        for photon in self.photons:
            if self.bin:
                u = photon.srcdiff()
            else:
                if self.useft2s:
                    u = photon.diff(self.rot)
                else:
                    u = photon.srcdiff()
            u = u*rd
            binfile.addPhoton(s.Photon(s.SkyDir(u,0),float(photon.energy),float(photon.time),int(self.cls)))
        binfile.write('%s_%1.0f_%1.0f_%1.0f.fits'%(self.lis,self.cls,self.emin,self.emax))


################################################### END ALIGNMENT CLASS ###########################################

################################################### TEST METHODS  #################################################


##   alignment test program
##   runs on a small set of flight data and uses Crab, Geminga, and Vela as sources
##   if everything is ok, should return 0
def test():
    plr = ''#os.environ['POINTLIKEROOT']
    fdir = plr+basedir+'mar0/pointlikedev/uw/stacklike/boresighttest/'
    al = StackLoader(lis='cgv',files=['test'],datadir=fdir,ft2dir=fdir,srcdir=fdir,quiet=False,irf='P6_v8_diff')
    al.loadphotons(0,4,1000,2e4,0,999999999,0)
    al.solverot()
    al.makeplot('aligntest')
    ret = 0
    if abs(al.params[0]+37)>10:
        ret = ret + 1
    if abs(al.params[1]-5)>10:
        ret = ret + 2
    if abs(al.params[2]+69)>10:
        ret = ret + 4
    if abs(al.errors[0]-54)>10:
        ret = ret+8
    if abs(al.errors[1]-51)>10:
        ret = ret+16
    if abs(al.errors[2]-84)>10:
        ret = ret +32
    return ret

##   psf parameter test program
##   runs on a small set of flight data and uses Crab, Geminga, and Vela as sources
##   if everything is ok, should return 0
def test2(bins=25.):
    plr = ''#os.environ['POINTLIKEROOT']
    fdir = plr+basedir+'mar0/pointlikedev/uw/stacklike/boresighttest/'
    os.system('cd %s'%fdir)
    al = StackLoader(lis='cgv',files=['test'],datadir=fdir,ft2dir=fdir,srcdir=fdir,quiet=False,useft2s=False)
    al.loadphotons(0,10,1000,1770,0,999999999,0)
    al.solvepsf()
    ret=0
    if (al.sigma*rd-0.149)>1e-3:
        ret = ret + 1
    if (al.gamma-1.63)>1e-2:
        ret = ret + 2 
    al.makeplot('psftest',bin=bins,scale='log',jac='square',xs='square')
    return ret
    

##   background parameter test program
##   runs on a small set of flight data and uses Crab, Geminga, and Vela as sources
##   if everything is ok, should return 0
def test3(bins=25.,stats=True):
    plr = ''#os.environ['POINTLIKEROOT']
    fdir = plr+basedir+'mar0/pointlikedev/uw/stacklike/boresighttest/'
    os.system('cd %s'%fdir)
    al = StackLoader(lis='cgv',files=['test'],datadir=fdir,ft2dir=fdir,srcdir=fdir,quiet=False,useft2s=False)
    al.loadphotons(0,10,1000,1770,0,999999999,0)
    al.solveback()
    ret=0
    al.makeplot('backtest',bin=bins,scale='log',jac='square',xs='square',stats=stats)
    return ret

##   background parameter test program
##   runs on a small set of flight data and uses Crab, Geminga, and Vela as sources
##   if everything is ok, should return 0
def test4(bins=25.,stats=True):
    plr = ''#os.environ['POINTLIKEROOT']
    fdir = plr+basedir+'mar0/pointlikedev/uw/stacklike/boresighttest/'
    os.system('cd %s'%fdir)
    al = StackLoader(lis='cgv',files=['test'],datadir=fdir,ft2dir=fdir,srcdir=fdir,quiet=False,useft2s=False)
    al.loadphotons(0,10,1000,1770,0,999999999,0)
    al.solveback()
    ret=0
    al.makeplot('backtest',bin=bins,scale='log',jac='square',xs='square',stats=stats)
    return ret
###################################################################################################