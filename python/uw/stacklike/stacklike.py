"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/stacklike/stacklike.py,v 1.6 2010/10/21 22:35:35 mar0 Exp $
author: M.Roth <mar0@u.washington.edu>
"""

import numpy as np
import pylab as py
import skymaps as s
import pointlike as pl
import pyfits as pf
import glob as glob
import scipy.integrate as si
import scipy.optimize as so
import os as os
from uw.like import pypsf
from uw.stacklike.angularmodels import *
from uw.stacklike.CLHEP import HepRotation,Hep3Vector,Photon
import matplotlib.font_manager


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
CALDBdir = r'd:\fermi/CALDB/v1r1/CALDB/data/glast/lat'
datadir = r'd:\fermi\data\flight/'                             #directory for FT1 files
ft2dir = r'd:\fermi\data\flight/'                              #directory for FT2 files
srcdir = r'd:\common\mar0\sourcelists/'                        #directory for source lists
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
        self.lis='strong'
        self.rot=[0,0,0]
        self.files=files
        self.CALDBdir=CALDBdir
        self.datadir=datadir
        self.srcdir=srcdir
        self.ft2dir=ft2dir
        self.binfile=''
        self.irf = 'P6_v3_diff'
        self.ctmin=0.4
        self.ctmax=1.0
        self.quiet=True
        self.firstlight = False
        self.tev = False
        self.bin=False
        self.useft2s=True
        self.__dict__.update(kwargs)

        print ''
        print '**********************************************************'
        print '*                                                        *'
        print '*                       STACKLIKE                        *'
        print '*                                                        *'
        print '**********************************************************'
        print 'Using irf: %s'%self.irf
        print 'Using list: %s.txt'%self.lis
        print 'Using boresight alignment (in arcsec): Rx=%1.0f Ry=%1.0f Rz=%1.0f'%(self.rot[0]*rd*3600,self.rot[1]*rd*3600,self.rot[2]*rd*3600)

        self.rot = HepRotation(self.rot,False)

        if self.tev:
            self.CALDBdir = r'/phys/groups/tev/scratch1/users/Fermi/CALDB/v1r1/CALDB/data/glast/lat/'
            self.datadir = r'/phys/groups/tev/scratch1/users/Fermi/data/flight/'
            self.ft2dir = r'/phys/users/mar0/python/'
            self.srcdir = r'/phys/users/mar0/python/'
            self.logfile = open(r'/phys/users/mar0/python/alignment_log.txt','w')

        if not self.firstlight:
            ftemp = [self.datadir + fil+'-ft1.fits' for fil in self.files]
            if self.useft2s:
                self.ft2s = [self.ft2dir + fil+'-ft2.fits' for fil in self.files]
            else:
                self.ft2s = []
            self.files=ftemp
        else:
            self.ft2dir = r'd:\common\mar0\data\firstlight/'
            self.files = [self.ft2dir+'jul2008-ft1.fits']
            self.ft2s = [self.ft2dir+'jul2008-ft2.fits']

        if self.binfile != '':
            self.bin = True

        #s.IParams.set_CALDB(self.CALDBdir)
        #s.IParams.init(self.irf)
        self.atb = []

        self.psf = pypsf.CALDBPsf(irf=self.irf)
        self.ds = []
        self.diffs =[]
        self.aeff = []
        self.photons = []
        self.offsrc = []
        self.photoncount=0
        self.Npsf =0
        self.Nback =0
        self.ebar=0
        self.cls=0
        self.loadsrclist(self.lis)

    def loadsrclist(self,fname):
        self.srcs = []
        sf = file(self.srcdir+fname+'.txt')
        header = sf.readline()
        for lines in sf:
            line = lines.split()
            ra = float(line[1])
            dec = float(line[2])
            sd = s.SkyDir(ra,dec)
            self.srcs.append(sd)

    ##  loads photons from FT1 file near sources  
    #   @param minroi min ROI in degrees
    #   @param maxroi max ROI in degrees
    #   @param emin minimum energy in MeV
    #   @param emax maximum energy in MeV
    #   @param start start time (MET)
    #   @param stop end time (MET)
    #   @param cls conversion type: 0=front,1=back,-1=all
    def loadphotons(self,minroi,maxroi,emin,emax,start,stop,cls):
        self.atb = []
        self.ebar=self.ebar*len(self.photons)
        self.emin=emin
        self.emax=emax
        self.minroi = minroi
        self.maxroi = maxroi
        self.cls = cls
        self.eave = np.sqrt(self.emin*self.emax)
        if self.bin:
            print 'Applying masks to data'
            print '**********************************************************'
            print '%1.0f < Energy < %1.0f'%(self.emin,self.emax)
            print 'Event class = %1.0f'%(cls)
            print '**********************************************************'

            self.bpd = s.BinnedPhotonData(self.datadir+self.binfile)
            self.ebar = 0
            for bnd in self.bpd:
                if bnd.emin()<self.eave and bnd.emax()>self.eave and (bnd.event_class()==cls or cls==-1):
                    ebar = np.sqrt(bnd.emin()*bnd.emax())
                    for src in self.srcs:
                        wsdl = s.WeightedSkyDirList(bnd,src,self.maxroi/rd)
                        sds = len(wsdl)
                        for wsd in wsdl:
                            self.photons.append(Photon(wsd.ra(),wsd.dec(),ebar,0,bnd.event_class(),[],[],src,wsd.weight()))
                            self.ebar = self.ebar + wsd.weight()*ebar
            pcnts = sum([x.weight for x in self.photons])
            self.photonscount = self.photoncount+pcnts
            if len(self.photons)==0:
                print 'No photons!'
                raise 'No photons!'
            else:
                self.ebar = self.ebar/pcnts
            print '%d photons remain'%pcnts
            print '**********************************************************'
        else:
            print 'Applying masks to data'
            print '**********************************************************'
            print '%1.0f < Energy (MeV) < %1.0f'%(self.emin,self.emax)
            print '%1.2f < Separation (deg) < %1.2f'%(self.minroi,self.maxroi)
            print '%1.0f < Time < %1.0f'%(start,stop)
            print 'Event class = %1.0f'%(cls)
            print '%1.2f < costh < %1.2f'%(self.ctmin,self.ctmax)
            print '**********************************************************'

            #go through each fits file, mask unwanted events, and setup tables
            tcuts = np.array([0,0,0,0,0,0,0,0])
            for ff in self.files:
                if not self.quiet:
                    print ff
                tff = pf.open(ff)
                ttb = tff[1].data
                ttb,ttcut = self.mask(ttb,self.srcs,emin,emax,start,stop,self.maxroi,cls)
                tcuts=tcuts+ttcut
                self.atb.append(ttb)
                tff.close()
            print 'Photon pruning information'
            print '%d photons available'%tcuts[0]
            print '---------------------------------------'
            print '%d cut by time range'%tcuts[1]
            print '%d cut by energy range'%tcuts[2]
            print '%d cut by instrument theta'%tcuts[3]
            print '%d cut by zenith angle'%tcuts[4]
            print '%d cut by conversion type'%tcuts[5]
            print '%d cut by proximity to sources'%tcuts[6]
            print '---------------------------------------'
            print '%d photons remain'%tcuts[7]
            print '**********************************************************'


            #go through each table and contruct Photon objects
            for j,tb in enumerate(self.atb):
                if not self.quiet:
                    print 'Examining %d events'%len(tb)
                if len(tb)>0:
                    if not self.quiet:
                        if self.useft2s:
                            print '    *Loading pointing history'
                            phist = pl.PointingHistory(self.ft2s[j])
                    if not self.quiet:
                        print '    *Loading events'
                    
                    #iterate over events for photons
                    for k in range(len(tb)):
                        event = tb[k]
                        sd = s.SkyDir(float(event.field('RA')),float(event.field('DEC')))
                        rsrc = self.srcs[0]
                        diff = 1e40

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
                                photon = Photon(sd.ra(),sd.dec(),event.field('ENERGY'),time,event.field('CONVERSION_TYPE'),xv,zv,src)
                                self.ebar = self.ebar+event.field('ENERGY')
                                self.photons.append(photon)
                                self.photoncount = self.photoncount + 1
                    if self.useft2s:
                        del phist #free up memory from pointing history object
            if len(self.photons)>0:
                self.ebar = self.ebar/len(self.photons)  #calculate mean energy of photons
            else:
                print 'No Photons!'
                raise 'No Photons!'

    ## estimates Isotropic component from PSF
    #  @param func optional custom background component to fit [[angular separations],[counts]]
    #  @param free free or fix model number estimators
    #  @param exp specify number estimators
    def solveback(self,func=[[],[]],free=[True,True],exp=[]):

        #Try to estimate the number of Isotropic photons
        #by looking at density near tails
        self.getds()
        frac = 0.98
        area = 1-frac*frac
        density = len(self.ds[self.ds>(frac*self.maxroi/rd)])
        denback = density*(self.maxroi*self.maxroi-self.minroi*self.minroi)/((self.maxroi**2)*area)

        if denback>len(self.ds):
            denback = len(self.ds)
        psf = self.getpsf()
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
                exp=[len(self.photons)-denback,denback]
            else:
                exp=[len(self.photons)-denback,0.,denback]
        if len(free)!=len(cm.models):
            free=[True,True,True]
        cm.fit(self.photons,mode=0,quiet=self.quiet,free=free,exp=exp)
        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Nbacke = np.sqrt(self.errs[1][1])
        if len(cm.models)==2:
            self.il=cm.extlikelihood([self.Npsf,self.Nback,psf.model_par[0],psf.model_par[1]],self.photons)
        else:
            self.Ncust = cm.minuit.params[2]
            self.Ncuste = np.sqrt(self.errs[2][2])
            self.il=cm.extlikelihood([self.Npsf,self.Nback,self.Ncust,psf.model_par[0],psf.model_par[1]],self.photons)
        self.cm = cm

    ## returns the TS for a point source in a uniform background
    #  @param ct optionally define the PSF by [R68,R95] in radians, default is the self.irf
    def TS(self,ct=[]):
        if len(self.photons)==0:
            return 0
        if (self.Npsf+self.Nback)==0:
            self.solveback()
        cm = CompositeModel()
        if ct==[]:
            psf = self.getpsf()
        else:
            psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.001,2.25],free=[True,True])
            psf.fromcontain(ct[0],ct[1],0.68,0.95)
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
            print 'R%d %3.0f +/-%3.0f arcsec'%(x+1,self.params[x],self.errors[x])
        try:
            TS = -2*(fl-self.il)
            print 'Significance of rotation was TS = %d'%TS
        except:
            print 'Cannot determine improvement'
        print 'Expected %d photons, got %d ( %d (%d) + %d (%d) )'%(len(self.photons),int(self.Npsf+self.Nback),int(self.Npsf),int(np.sqrt(self.Npsfe)),int(self.Nback),int(np.sqrt(self.Nbacke)))
        print 'Called likelihood %d times'%cm.calls
        self.cm=cm

    ## tries to solve parameters for a single King function with Isotropic background
    #  @param func optional custom background component to fit [[angular separations],[counts]]
    #  @param free free or fix model number estimators
    #  @param exp specify number estimators
    def solvepsf(self,func=[[],[]],free=[True,True],exp=[]):

        #estimate Isotropic component
        self.solveback()

        #solve for PSF parameters
        psf = self.getpsf(True)
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

        #fit all free parameters
        cm.fit(self.photons,free=free,exp=exp,mode=1,quiet=self.quiet)
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
        self.cov = self.errs[nm][nm+1]
        if len(cm.models)==2:
            self.il=cm.extlikelihood([self.Npsf,self.Nback,sigma,gamma],self.photons)
        else:
            self.il=cm.extlikelihood([self.Npsf,self.Nback,0.,sigma,gamma],self.photons)
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
        print '**********************************************************'
        print 'Npsf  = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100)
        print 'Nback = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f2*100,f2e*100)
        print 'Sigma = %1.3f [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.sigma*rd,self.sigmae/self.sigma,self.irf,self.sigma/sigma)
        print 'Gamma = %1.2f  [1 +/- %1.2f]        Ratio to %s: (%1.2f)'%(self.gamma,self.gammae/self.gamma,self.irf,self.gamma/gamma)
        print 'R68   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r68*rd,self.r68e/self.r68,self.irf,self.r68/r68o)
        print 'R95   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r95*rd,self.r95e/self.r95,self.irf,self.r95/r95o)
        TS = -2*(fl-self.il)
        print 'Significance of psf change was TS = %d'%TS
        self.cm=cm

    # tries to fit a halo component on top of a PSF defined by 'irf' in a Isotropic background
    #  @param cust [sigma,gamma] of PSF, where sigma is in radians
    def solvehalo(self,cust=[]):

        #estimate Isotropic component
        self.solveback()

        #solve for Halo model component while freezing PSF parameters

        if cust==[]:
            psf = self.getpsf()
        else:
            psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=cust,free=[False,False])
        halo = Halo(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.2/rd])
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(halo)
        cm.addModel(bck)
        cm.fit(self.photons,free=[True,True,True],exp=[self.Npsf/2.,self.Npsf/2.,self.Nback],mode=1,quiet=self.quiet)
        self.errs = cm.minuit.errors()
        self.Npsf = cm.minuit.params[0]
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Nhalo = cm.minuit.params[1]
        self.Nhaloe = np.sqrt(self.errs[1][1])
        self.Nback = cm.minuit.params[2]
        self.Nbacke = np.sqrt(self.errs[2][2])
        self.theta = cm.minuit.params[5]
        self.thetae = np.sqrt(self.errs[5][5])/self.theta
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
        print '**********************************************************'
        print 'Npsf  = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100)
        print 'Nhalo = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Nhalo,self.Nhaloe/self.Nhalo,f2*100,f2e*100)
        print 'Nback = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f3*100,f3e*100)
        print 'Halo width was %1.3f [1 +/- %1.2f] deg'%(self.theta*rd,self.thetae)
        print 'Halo fraction was %1.0f [1 +/- %1.2f]'%(self.frac,self.frace)
        self.cm=cm

    ## tries to fit two King functions in an Isotropic background
    #  @param func optional custom background component to fit [[angular separations],[counts]]
    #  @param free free or fix model number estimators
    #  @param exp specify number estimators
    def solvedoublepsf(self,func=[[],[]],free=[True,True,True],exp=[]):
        
        #estimate uniform Isotropic component
        self.solveback()

        #solve two King function parameters
        psf = self.getpsf(True)
        sigma = psf.model_par[0]
        gamma = psf.model_par[1]
        psf2 = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=[1.2*psf.model_par[0],0.7*psf.model_par[1]]) #try to separate psfs initially
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(psf2)
        cm.addModel(bck)
        if len(func[0])!=0:
            bck2 = Custom(lims=[self.minroi/rd,self.maxroi/rd],func=func)
            cm.addModel(bck2)
        if exp==[]: 
            if cm.models==3:
                exp=[self.Npsf/2.,self.Npsf/2.,self.Nback]
                free = [True,True,True]
            else:
                exp=[self.Npsf/2.,self.Npsf/2.,self.Nback,0.]
                free=[True,True,True,True]
        cm.fit(self.photons,free=free,exp=exp,mode=1,quiet=self.quiet)
        nm = len(cm.models)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        self.Npsf2 = cm.minuit.params[1]
        self.Nback = cm.minuit.params[2]
        self.errs = cm.minuit.errors()
        self.Npsfe = np.sqrt(self.errs[0][0])
        self.Npsf2e = np.sqrt(self.errs[1][1])
        self.Nbacke = np.sqrt(self.errs[2][2])
        if nm==4:
            self.Ncust = cm.minuit.params[3]
            self.Ncuste = np.sqrt(self.errs[3][3])
        self.sigma = cm.minuit.params[nm]
        self.sigmae = np.sqrt(self.errs[nm][nm])
        self.gamma = cm.minuit.params[nm+1]
        self.gammae = np.sqrt(self.errs[nm+1][nm+1])
        self.sigma2 = cm.minuit.params[nm+2]
        self.sigmae2 = np.sqrt(self.errs[nm+2][nm+2])
        self.gamma2 = cm.minuit.params[nm+3]
        self.gammae2 = np.sqrt(self.errs[nm+3][nm+3])
        if nm==3:
            self.il=cm.extlikelihood([self.Npsf,self.Npsf2,self.Nback,sigma,gamma,sigma,gamma],self.photons)
        else:
            self.il=cm.extlikelihood([self.Npsf,self.Npsf2,self.Nback,0.,sigma,gamma,sigma,gamma],self.photons)

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
        phs = (self.Npsf+self.Nback+self.Npsf2)
        tr = sum([self.errs[i][i] for i in range(len(self.errs))])
        f1 = self.Npsf/phs
        f1e = f1*np.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        f2 = self.Npsf2/phs
        f2e = f2*np.sqrt(self.errs[1][1]/(self.Npsf2**2)+tr/(phs**2))
        f3 = self.Nback/phs
        f3e = f3*np.sqrt(self.errs[2][2]/(self.Nback**2)+tr/(phs**2))
        print '***************************************'
        print 'Npsf   = %1.0f  [1 +/- %1.2f]     Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100)
        print 'Npsf2  = %1.0f  [1 +/- %1.2f]     Fraction: %1.0f +/- %1.0f'%(self.Npsf2,self.Npsf2e/self.Npsf2,f2*100,f2e*100)
        print 'Nback  = %1.0f  [1 +/- %1.2f]     Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f3*100,f3e*100)
        print 'Sigma  = %1.3f [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.sigma*rd,self.sigmae/self.sigma,ssign,100*abs((self.sigma-sigma)/sigma))
        print 'Gamma  = %1.2f  [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.gamma,self.gammae/self.gamma,gsign,100*abs((self.gamma-gamma)/gamma))
        print 'Sigma2 = %1.3f [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.sigma2*rd,self.sigmae2/self.sigma2,ssign2,100*abs((self.sigma2-sigma)/sigma))
        print 'Gamma2 = %1.2f  [1 +/- %1.2f]      Fractional Change: (%s%1.0f)'%(self.gamma2,self.gammae2/self.gamma2,gsign2,100*abs((self.gamma2-gamma)/gamma))
        TS = -2*(fl-self.il)
        print 'Significance of psf change was TS = %d'%TS
        self.cm=cm

    def solvediffuse(self,fit=False):
        #Try to estimate the number of Isotropic photons
        #by looking at density near tails
        self.getds()
        frac = 0.98
        area = 1-frac*frac
        density = len(self.ds[self.ds>(frac*self.maxroi/rd)])
        denback = density*(self.maxroi*self.maxroi-self.minroi*self.minroi)/((self.maxroi**2)*area)

        if denback>len(self.ds):
            denback = len(self.ds)
        exp=[len(self.ds)-denback,denback/2.,denback/2.]
        #exp = [len(self.ds)/2.,len(self.ds)/2.,0]
        diffusefile = r'd:\fermi\data\galprop\ring_21month_v1.fits'

        psf = self.getpsf(opt=fit)
        sigma = psf.model_par[0]
        gamma = psf.model_par[1]
        bck = Diffuse(lims=[self.minroi/rd,self.maxroi/rd],diff=diffusefile,ebar=self.ebar)
        bck2 = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        offs = []
        for off in self.offsrc:
            toff = OffPsf(lims=[self.minroi/rd,self.maxroi/rd],model_par=[psf.model_par[0],psf.model_par[1]],off=off)
            offs.append(toff)
        #bck3 = Gaussian(lims=[self.minroi/rd,self.maxroi/rd],model_par=[1./rd])
        for src in self.srcs:
            bck.addsrc(src)
            for tof in offs:
                tof.addsrc(src)
        
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)
        cm.addModel(bck2)
        for tof in offs:
            cm.addModel(tof)
        nm = len(cm.models)
        if nm<=3:
            cm.fit(self.photons,mode=0,quiet=self.quiet,free=[True,True,True],exp=exp)
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
            print '**********************************************************'
            print 'Npsf  = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100)
            print 'Nback = %1.0f [1 +/- %1.2f]      Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f2*100,f2e*100)
            print 'Sigma = %1.3f [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.sigma*rd,self.sigmae/self.sigma,self.irf,self.sigma/sigma)
            print 'Gamma = %1.2f  [1 +/- %1.2f]        Ratio to %s: (%1.2f)'%(self.gamma,self.gammae/self.gamma,self.irf,self.gamma/gamma)
            print 'R68   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r68*rd,self.r68e/self.r68,self.irf,self.r68/r68o)
            print 'R95   = %1.2f  [1 +/- %1.2f] (deg)  Ratio to %s: (%1.2f)'%(self.r95*rd,self.r95e/self.r95,self.irf,self.r95/r95o)
            TS = -2*(fl-self.il)
            print 'Significance of psf change was TS = %d'%TS
        else:
            left = len(self.ds)-denback
            exp = [left/(nm-2.) for x in range(nm)]
            free = [True for x in range(nm)]
            exp[1]=denback/2.
            exp[2]=denback/2.
            cm.fit(self.photons,mode=0,quiet=self.quiet,free=free,exp=exp)
            self.cm = cm


    def solvehalolims(self,frac):
        self.solveback()
        psf = self.getpsf(True)
        bck = Isotropic(lims=[self.minroi/rd,self.maxroi/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)
        cm.fit(self.photons,free=[True,True],exp=[self.Npsf,self.Nback],quiet=self.quiet,mode=1)
        psf2 = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=psf.model_par,free=[False,False])
        initial = cm.minuit.fval

        uplim=[]
        thetas = np.arange(-16,6,1)
        thetas = thetas/8.
        thetas = 10**(thetas)
        py.ioff()
        py.figure(figsize=(8,8))
        for theta in thetas:
            halo = Halo(lims=[self.minroi/rd,self.maxroi/rd],model_par=[theta/rd],free=[False])
            cm = CompositeModel()
            cm.addModel(psf2)
            cm.addModel(bck)
            cm.addModel(halo)
            cm.fit(self.photons,free=[True,True,False],exp=[self.Npsf,self.Nback,0],quiet=self.quiet,mode=0)
            #self.il=cm.minuit.fval
            #self.errors=cm.minuit.errors()
            #self.Nhalo=cm.minuit.params[2]
            def like(x):
                if x<0 or x>self.Npsf:
                    return 0.
                val = np.exp(cm.minuit.fval-cm.extlikelihood([self.Npsf-x,self.Nback,x,psf2.model_par[0],psf2.model_par[1],halo.model_par[0]],self.photons,True))
                #print val
                return val
            
            def integ(x):
                val = si.quad(lambda y:like(y),0,x)[0]
                return val

            def contain(x,fnt,ct,limit=False):
                ft = integ(x)
                if ft/fnt>ct:
                    return 1.
                val = abs(ft/fnt-ct)
                print x,val
                return val

            fint = integ(self.Npsf)
            
            best = so.fmin_powell(lambda x:contain(x,fint,frac),0,maxiter=2,full_output=1,disp=1)
            upl = best[0]#self.Nhalo+np.sqrt(self.errors[2][2])*2
            uplim.append(upl)
            print upl
            cm.nest=[self.Npsf-best[0],self.Nback,best[0]]
            self.cm=cm
            self.makeplot('uplims_%1.3f_%1.0f_%s_%s.png'%(theta,self.ebar,self.lis,self.irf),fig=1)
        return thetas,uplim

    ## Makes a plot of all of the models and the angular distribution of photons
    #  @param name filename of output file
    #  @param fig fignum of pylab figure, default is 1
    #  @param bin number of bins for unbinned data
    #  @param scale makes the x-scale log/linear  ['log','linear']
    #  @param jac either counts/deg or counts/(deg**2) ['square','linear']
    #  @param xs x-scale either deg or deg**2 ['linear','square']
    def makeplot(self,name='test.png',fig=-1,bin=25.,scale='log',jac='square',xs='linear',xlims=[]):
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
                jmin = min(self.ds)**2
                jmax = max(self.ds)**2
            else:
                jmin = min(self.ds) 
                jmax = max(self.ds)
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
        

        #retrieve all model names
        for model in range(len(self.cm.models)):
            fits.append([])
            names.append(self.cm.models[model].name)
            if self.cm.models[model].name=='psf':
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
        py.subplot(2,1,1)
        #plot histogrammed events
        py.errorbar(be4,hists,yerr=errs,ls='None',marker='None',ecolor='blue',capsize=80./bin)
        p1 = py.step(be4,hists,where='mid')
        pts.append(p1[0])

        pcts=0

        #plot all models
        for model in range(len(self.cm.models)):
            p1 = py.plot(be4,fits[model],self.cm.models[model].mark)
            pts.append(p1)
            if self.cm.models[model].name=='psf':
                p1 = py.plot([ctn[2*pcts],ctn[2*pcts]],[pmin,pmax],'--')
                pts.append(p1)
                p1 = py.plot([ctn[2*pcts+1],ctn[2*pcts+1]],[pmin,pmax],'--')
                pts.append(p1)
                pcts=pcts+1
        p1 = py.plot(be4,fits[len(fits)-1],self.cm.mark)
        pts.append(p1)

        #calculate chisq
        self.chisq=0.
        total = fits[len(fits)-1]
        for it,bin in enumerate(hists):
            if errs[it]>1e-35:
                self.chisq = self.chisq + ((bin-total[it])/errs[it])**2

        prop = matplotlib.font_manager.FontProperties(size=9) 
        #finish up plotting
        py.legend(pts,names,bbox_to_anchor=(1.1, 1.25),prop=prop)
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
        py.suptitle('%s: Energy(MeV): [%1.0f,%1.0f]\nClass: %s    Costh: [%1.1f,%1.1f]\nChisq: %1.1f / %1.0f dof'%(self.lis,self.emin,self.emax,clas[self.cls],self.ctmin,self.ctmax,self.chisq,len(be2)-1))
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
        py.savefig(name)
        


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
        cuts = [0,0,0,0,0,0]
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
                    msk = table.field('ZENITH_ANGLE')<105
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
        #display number of photons cut from each step and finally the number of photons examined
        if not self.quiet:
            print 'TOTAL: %d    TIMEC: %d    ENERGY: %d    THETA: %d    ZENITH: %d    ECLASS: %d    POS: %d    EXAM: %d'%(total,cuts[0],cuts[1],cuts[2],cuts[3],cuts[4],cuts[5],tc)
        return table,np.array([total,cuts[0],cuts[1],cuts[2],cuts[3],cuts[4],cuts[5],tc])
    

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
                if self.bin:
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

    ## returns a PSF object based on containment specified by self.irf
    #  @param opt frees or fixes PSF parameters
    def getpsf(self,opt=False):
        self.psf = pypsf.CALDBPsf(irf=self.irf)
        if self.cls==-1:
            r68f=self.psf.inverse_integral(self.ebar,0,68.)/rd
            r95f=self.psf.inverse_integral(self.ebar,0,95.)/rd
            r68b=self.psf.inverse_integral(self.ebar,1,68.)/rd
            r95b=self.psf.inverse_integral(self.ebar,1,95.)/rd
            r68=np.sqrt(r68f*r68b)
            r95=np.sqrt(r95f*r95b)
        else:
            r68=self.psf.inverse_integral(self.ebar,int(self.photons[0].event_class),68.)/rd
            r95=self.psf.inverse_integral(self.ebar,int(self.photons[0].event_class),95.)/rd
        psf = PSF(lims=[self.minroi/rd,self.maxroi/rd],model_par=[0.001,2.25],free=[opt,opt])
        factor = 2.
        psf.fromcontain(1./factor,r95/r68/factor,0.68,0.95)
        psf.model_par[0]=psf.model_par[0]*r68*factor
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
    plr = os.environ['POINTLIKEROOT']
    fdir = plr+'/python/uw/stacklike/boresighttest/'
    al = StackLoader(lis='cgv',files=['test'],datadir=fdir,ft2dir=fdir,srcdir=fdir,quiet=False,irf='P6_v8')
    al.loadphotons(0,4,1000,2e6,0,999999999,0)
    al.solverot()
    al.makeplot('aligntest.png')
    ret = 0
    if abs(al.params[0]+31)>10:
        ret = ret + 1
    if abs(al.params[1]-9)>10:
        ret = ret + 2
    if abs(al.params[2]+45)>10:
        ret = ret + 4
    if abs(al.errors[0]-60)>10:
        ret = ret+8
    if abs(al.errors[1]-59)>10:
        ret = ret+16
    if abs(al.errors[2]-96)>10:
        ret = ret +32
    return ret

##   psf parameter test program
##   runs on a small set of flight data and uses Crab, Geminga, and Vela as sources
##   if everything is ok, should return 0
def test2(bins=25.):
    plr = os.environ['POINTLIKEROOT']
    fdir = plr+'/python/uw/stacklike/boresighttest/'
    os.system('cd %s'%fdir)
    al = StackLoader(lis='cgv',files=['test'],datadir=fdir,ft2dir=fdir,srcdir=fdir,quiet=False,useft2s=False)
    al.loadphotons(0,10,1000,1770,0,999999999,0)
    al.solvepsf()
    ret=0
    if (al.sigma*rd-0.149)>1e-3:
        ret = ret + 1
    if (al.gamma-1.63)>1e-2:
        ret = ret + 2 
    al.makeplot('psftest.png',bin=bins,scale='log',jac='square',xs='square')
    return ret
    
###################################################################################################