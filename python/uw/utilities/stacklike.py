"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/stacklike.py,v 1.1 2010/08/02 20:48:10 mar0 Exp $
author: M.Roth <mar0@u.washington.edu>
"""

import numpy as N
import pylab as p
import skymaps as s
import pointlike as pl
import pyfits as pf
import glob as glob
import scipy.integrate as si
import os as os
from uw.utilities.angularmodels import CompositeModel,PSF,PSFAlign,Backg,Halo
from uw.utilities.CLHEP import HepRotation,Hep3Vector,Photon

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
rd = 180./N.pi        #to radians
monthlist = ['aug2008','sep2008','oct2008','nov2008','dec2008','jan2009','feb2009','mar2009','apr2009'\
    ,'may2009','jun2009','jul2009','aug2009','sep2009','oct2009','nov2009','dec2009','jan2010','feb2010','mar2010'\
    ,'apr2010','may2010','jun2010','jul2010']
CALDBdir = r'y:\fermi/CALDB/v1r1/CALDB/data/glast/lat'
datadir = r'y:\fermi\data\flight/'                             #directory for FT1 files
ft2dir = r'y:\fermi\data\flight/'                              #directory for FT2 files
srcdir = r'y:\common\mar0\sourcelists/'                        #directory for source lists
files = []                                                     #default list of FT1 file names (minus '-ft1.fits')
ft2s = []                                                      #default list of FT2 file names (minus '-ft2.fits')
for month in monthlist:
    files.append('%s'%month)
    ft2s.append('%s'%month)

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
#   solvepsf() - fits the parameters of a single King function in a uniform background
#
#   solvehalo() - fits the parameters of a gaussain halo with a PSF and uniform background
#
#   solvedoublepsf() - fits the parameters of two King functions in a uniform background
#   !!!Warning: large degeneracy in on-orbit data, only works with very low uniform background!!!
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
    def __init__(self,lis='strong',tev=False,rot=[0,0,0],files=files,CALDBdr=CALDBdir,datadr=datadir,ft2dr=ft2dir,srcdr=srcdir,irf = 'P6_v3'):
        
        self.firstlight = False
        self.rot = HepRotation(rot,False)
        print 'Using list: %s.txt'%lis
        print 'Using boresight alignment (in arcsec): Rx=%1.0f Ry=%1.0f Rz=%1.0f'%(rot[0]*rd*3600,rot[1]*rd*3600,rot[2]*rd*3600)
        self.tev = tev
        self.CALDBdir = CALDBdr
        self.datadir = datadr
        self.srcdir = srcdr
        self.ctmin=0.4
        self.ctmax=1.0
        self.irf=irf

        if tev:
            self.CALDBdir = r'/phys/groups/tev/scratch1/users/Fermi/CALDB/v1r1/CALDB/data/glast/lat/'
            self.datadir = r'/phys/groups/tev/scratch1/users/Fermi/data/flight/'
            self.ft2dir = r'/phys/users/mar0/python/'
            self.srcdir = r'/phys/users/mar0/python/'
            self.logfile = open(r'/phys/users/mar0/python/alignment_log.txt','w')

        if not self.firstlight:
            self.ft2dir = ft2dr
            self.files = [self.datadir + fil+'-ft1.fits' for fil in files]
            self.ft2s = [self.ft2dir + fil+'-ft2.fits' for fil in files]
        else:
            self.ft2dir = r'd:\common\mar0\data\firstlight/'
            self.files = [self.ft2dir+'jul2008-ft1.fits']
            self.ft2s = [self.ft2dir+'jul2008-ft2.fits']


        s.IParams.set_CALDB(self.CALDBdir)
        s.IParams.init(self.irf)
        self.atb = []
        self.srcs = []
        sf = file(self.srcdir+lis+'.txt')
        header = sf.readline()
        for lines in sf:
            line = lines.split()
            ra = float(line[1])
            dec = float(line[2])
            sd = s.SkyDir(ra,dec)
            self.srcs.append(sd)

    ##  loads photons from FT1 file near sources  
    #   @param rad ROI in degrees
    #   @param emin minimum energy in MeV
    #   @param emax maximum energy in MeV
    #   @param start start time (MET)
    #   @param stop end time (MET)
    #   @param cls conversion type: 0=front,1=back,-1=all
    def loadphotons(self,rad,emin,emax,start,stop,cls):
        self.emin=emin
        self.emax=emax
        self.ebar=0
        self.rad = rad
        self.cls = cls
        self.us = []
        self.atb = []
        self.diffs =[]
        self.aeff = []
        
        print 'Applying masks to data'
        print '*********************************'
        print '%1.0f < Energy < %1.0f'%(self.emin,self.emax)
        print '%1.2f degrees from source'%(self.rad)
        print '%1.0f < Time < %1.0f'%(start,stop)
        print 'Event class = %1.0f'%(cls)
        print '%1.2f < costh < %1.2f'%(self.ctmin,self.ctmax)
        print '*********************************'

        #go through each fits file, mask unwanted events, and setup tables
        for ff in self.files:
            print ff
            tff = pf.open(ff)
            ttb = tff[1].data
            ttb = self.mask(ttb,self.srcs,emin,emax,start,stop,rad,cls)
            self.atb.append(ttb)
            tff.close()
        self.photons = []

        #go through each table and contruct Photon objects
        for j,tb in enumerate(self.atb):
            print 'Examining %d events'%len(tb)
            if len(tb)>0:
                print '    *Loading pointing history'
                phist = pl.PointingHistory(self.ft2s[j])
                print '    *Loading events'
                
                #iterate over events for photons
                for k in range(len(tb)):
                    event = tb[k]
                    sd = s.SkyDir(float(event.field('RA')),float(event.field('DEC')))
                    rsrc = self.srcs[0]
                    diff = 1e40

                    #associate with a particular source
                    for src in self.srcs:
                        tdiff = sd.difference(src)*rd
                        if tdiff<diff:
                            rsrc=src
                            diff=tdiff
                    
                    #if photons is not too far away from any source
                    #add it to the list
                    if diff<rad:
                        time = event.field('TIME')
                        pi = phist(time)
                        xax = pi.xAxis()
                        zax = pi.zAxis()
                        zen = pi.zenith()
                        xv = Hep3Vector([])
                        zv = Hep3Vector([])
                        xv(xax)
                        zv(zax)
                        photon = Photon(sd.ra(),sd.dec(),event.field('ENERGY'),time,event.field('EVENT_CLASS'),xv,zv,rsrc)
                        self.ebar = self.ebar+event.field('ENERGY')
                        self.photons.append(photon)

                del phist #free up memory from pointing history object
        self.ebar = self.ebar/len(self.photons)  #calculate mean energy of photons

    ## estimates uniform background component from PSF
    def solveback(self):
        #Try to estimate the number of background photons
        sigma=s.IParams.sigma(self.ebar,int(self.photons[0].event_class))
        gamma=s.IParams.gamma(self.ebar,int(self.photons[0].event_class))
        psf = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma],free=[False,False])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)
        cm.fit(self.photons,mode=0)

        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = N.sqrt(self.errs[0][0])
        self.Nbacke = N.sqrt(self.errs[1][1])
        
    ## Finds boresight alignment solution
    def solverot(self):
        psfa = PSFAlign(lims=[0,self.rad/rd],free=[True,True,True],ebar=self.ebar)
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psfa)
        cm.addModel(bck)
        cm.fit(self.photons)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = N.sqrt(self.errs[0][0])
        self.Nbacke = N.sqrt(self.errs[1][1])
        self.params=[x*rd*3600 for x in cm.minuit.params[2:]]
        self.errors=[N.sqrt(self.errs[x+2][x+2])*rd*3600 for x in range(3)]
        self.il=cm.extlikelihood([self.Npsf,self.Nback,0,0,0],self.photons)
        for x in range(3):
            print 'R%d %3.0f +/-%3.0f arcsec'%(x+1,self.params[x],self.errors[x])
        try:
            TS = -2*(fl-self.il)
            print 'Significance of rotation was TS = %d'%TS
        except:
            print 'Cannot determine improvement'
        print 'Expected %d photons, got %d ( %d (%d) + %d (%d) )'%(len(self.photons),int(self.Npsf+self.Nback),int(self.Npsf),int(N.sqrt(self.Npsfe)),int(self.Nback),int(N.sqrt(self.Nbacke)))
        print 'Called likelihood %d times'%cm.calls
        self.cm=cm

    ## tries to solve parameters for a single King function in a uniform background
    def solvepsf(self):

        #estimate uniform background component
        self.solveback()

        #solve for PSF parameters
        sigma=s.IParams.sigma(self.ebar,int(self.photons[0].event_class))
        gamma=s.IParams.gamma(self.ebar,int(self.photons[0].event_class))
        psf = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)
        cm.fit(self.photons,free=[True,True],exp=[self.Npsf,self.Nback],mode=1)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        self.Nback = cm.minuit.params[1]
        self.errs = cm.minuit.errors()
        self.Npsfe = N.sqrt(self.errs[0][0])
        self.Nbacke = N.sqrt(self.errs[1][1])

        self.sigma = cm.minuit.params[2]
        self.sigmae = N.sqrt(self.errs[2][2])
        self.gamma = cm.minuit.params[3]
        self.gammae = N.sqrt(self.errs[3][3])
        self.cov = self.errs[2][3]

        self.il=cm.extlikelihood([self.Npsf,self.Nback,sigma,gamma],self.photons)

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
        f1e = f1*N.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        f2 = self.Nback/phs
        f2e = f2*N.sqrt(self.errs[1][1]/(self.Nback**2)+tr/(phs**2))
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

    # tries to fit a halo component on top of a PSF defined by 'irf' in a uniform background
    def solvehalo(self):

        #estimate uniform background component
        self.solveback()

        #solve for Halo model component while freezing PSF parameters
        sigma=s.IParams.sigma(self.ebar,int(self.photons[0].event_class))
        gamma=s.IParams.gamma(self.ebar,int(self.photons[0].event_class))
        psf = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma],free=[False,False])
        halo = Halo(lims=[0,self.rad/rd],model_par=[0.2/rd])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(halo)
        cm.addModel(bck)
        cm.fit(self.photons,free=[True,True,True],exp=[self.Npsf/2.,self.Npsf/2.,self.Nback],mode=1)
        self.errs = cm.minuit.errors()
        self.Npsf = cm.minuit.params[0]
        self.Npsfe = N.sqrt(self.errs[0][0])
        self.Nhalo = cm.minuit.params[1]
        self.Nhaloe = N.sqrt(self.errs[1][1])
        self.Nback = cm.minuit.params[2]
        self.Nbacke = N.sqrt(self.errs[2][2])
        self.theta = cm.minuit.params[5]
        self.frac = self.Nhalo*100./(self.Nhalo+self.Npsf)
        ferr = self.frac*N.sqrt(self.Nhaloe/(self.Nhalo**2)+self.Npsfe/(self.Npsf**2))
        phs = (self.Npsf+self.Nback+self.Nhalo)
        tr = sum([self.errs[i][i] for i in range(len(self.errs))])
        f1 = self.Npsf/phs
        f1e = f1*N.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        f2 = self.Nhalo/phs
        f2e = f2*N.sqrt(self.errs[1][1]/(self.Nhalo**2)+tr/(phs**2))
        f3 = self.Nback/phs
        f3e = f3*N.sqrt(self.errs[2][2]/(self.Nback**2)+tr/(phs**2))
        print '**********************************************************'
        print 'Npsf  = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Npsf,self.Npsfe/self.Npsf,f1*100,f1e*100)
        print 'Nhalo = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Nhalo,self.Nhaloe/self.Nhalo,f2*100,f2e*100)
        print 'Nback = %1.0f [1 +/- %1.2f]  Fraction: %1.0f +/- %1.0f'%(self.Nback,self.Nbacke/self.Nback,f3*100,f3e*100)
        print 'Halo width was %1.3f [1 +/- %1.2f] deg'%(self.theta*rd,N.sqrt(self.errs[5][5])/self.theta)
        print 'Halo fraction was %1.0f [1 +/- %1.2f]'%(self.frac,ferr/self.frac)
        self.cm=cm

    ## tries to fit two King functions in a uniform background
    def solvedoublepsf(self):
        
        #estimate uniform background component
        self.solveback()

        #solve two King function parameters
        sigma=s.IParams.sigma(self.ebar,int(self.photons[0].event_class))
        gamma=s.IParams.gamma(self.ebar,int(self.photons[0].event_class))
        psf = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma])
        psf2 = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(psf)
        cm.addModel(bck)
        cm.fit(self.photons,free=[True,True,True],exp=[self.Npsf/2.,self.Npsf/2.,self.Nback],mode=1)
        fl = cm.minuit.fval
        self.Npsf = cm.minuit.params[0]
        self.Npsf2 = cm.minuit.params[1]
        self.Nback = cm.minuit.params[2]
        self.errs = cm.minuit.errors()
        self.Npsfe = N.sqrt(self.errs[0][0])
        self.Npsf2e = N.sqrt(self.errs[1][1])
        self.Nbacke = N.sqrt(self.errs[2][2])

        self.sigma = cm.minuit.params[3]
        self.sigmae = N.sqrt(self.errs[3][3])
        self.gamma = cm.minuit.params[4]
        self.gammae = N.sqrt(self.errs[4][4])
        self.sigma2 = cm.minuit.params[5]
        self.sigmae2 = N.sqrt(self.errs[5][5])
        self.gamma2 = cm.minuit.params[6]
        self.gammae2 = N.sqrt(self.errs[6][6])

        self.il=cm.extlikelihood([self.Npsf,self.Npsf2,self.Nback,sigma,gamma,sigma,gamma],self.photons)

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
        f1e = f1*N.sqrt(self.errs[0][0]/(self.Npsf**2)+tr/(phs**2))
        f2 = self.Npsf2/phs
        f2e = f2*N.sqrt(self.errs[1][1]/(self.Npsf2**2)+tr/(phs**2))
        f3 = self.Nback/phs
        f3e = f3*N.sqrt(self.errs[2][2]/(self.Nback**2)+tr/(phs**2))
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

    ## Makes a plot of all of the models and the angular distribution of photons
    #  @param name filename of output file
    def makeplot(self,name='test.png',fig=-1):
        #calculate angular separations
        self.getds()

        #plotting setup
        bins = 25.                       #angular bins
        mi = N.log10((min(self.ds))**2)
        ma = N.log10((max(self.ds))**2)
        d2 = ma-mi
        be = N.arange(mi,ma,d2/bins)
        be2 = N.append(be,be[len(be)-1]+d2/bins)
        be3 = N.sqrt(10**(be2))
        be4 = 10**(be+d2/bins/2)*rd*rd
        p.ioff()
        if fig ==-1:
            p.figure(figsize=(8,8))
        else:
            p.clf()
        hist = p.hist(N.log10(self.ds*self.ds),bins=be2,fc='None')
        p.clf()
        fits=[]
        hists=[]
        errs=[]
        names=[]                      #names of the various plots
        names.append('data')
        
        #retrieve all model names
        for model in range(len(self.cm.models)):
            fits.append([])
            names.append(self.cm.models[model].name)
        fits.append([])
        names.append('total')

        #fill all model plots
        for it,ba in enumerate(be):
            left=be3[it]
            right=be3[it+1]
            for it2,model in enumerate(self.cm.models):
                fits[it2].append(self.cm.nest[it2]*model.integral(left,right)/model.integral(0,self.rad/rd)/(right**2-left**2)/(rd**2))
            fits[len(fits)-1].append(self.cm.integral(left,right,0,self.rad/rd)/(right**2-left**2)/(rd**2))
            hists.append(hist[0][it]/(right**2-left**2)/(rd**2))
            errs.append(N.sqrt(hist[0][it])/(right**2-left**2)/(rd**2))
        pts=[]

        #plot histogrammed events
        p1 = p.errorbar(be4,hists,yerr=errs,ls='None',marker='o')
        pts.append(p1[0])

        #plot all models
        for ad in fits:
            p1 = p.plot(be4,ad,'-')
            pts.append(p1)

        #finish up plotting
        p.legend(pts,names)
        ax = p.gca()
        ax.set_yscale("log", nonposy='clip')
        p.semilogx()
        p.ylim(0.01*min(fits[len(fits)-1]),100*max(fits[len(fits)-1]))
        p.xlim(0.8*min(be4),max(be4)*1.2)
        p.xlabel(r'$\theta^2\/(\rm{deg^2})$')
        p.ylabel(r'$dN/d\theta^2$')
        p.grid()
        p.savefig(name)
        


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
                msk = (table.field('THETA')<(N.arccos(self.ctmin)*rd)) & (table.field('THETA')>(N.arccos(self.ctmax)*rd))
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
        print 'TOTAL: %d    TIMEC: %d    ENERGY: %d    THETA: %d    ZENITH: %d    ECLASS: %d    POS: %d    EXAM: %d'%(total,cuts[0],cuts[1],cuts[2],cuts[3],cuts[4],cuts[5],tc)
        return table
    

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
                cdec = N.cos((dec-deg)/rd)
                ramin = ra-deg/cdec
                ramax = ra+deg/cdec
                decmin = dec-deg
                decmax = dec+deg
            else:
                cdec = N.cos((dec+deg)/rd)
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
    def getus(self,x1,y1,z1):
        self.us=[]
        x,y,z=x1/rd/3600.,y1/rd/3600.,z1/rd/3600.
        rot = HepRotation([x,y,z],False)
        for photon in self.photons:
            sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
            umax = (self.rad/sig/rd)
            umax = umax*umax/2
            u = photon.diff(rot)/sig
            u = u*u/2
            self.us.append(u)

    #sets angular deviations
    def getds(self):
        self.ds=[]
        for photon in self.photons:
            u = photon.diff(self.rot)
            self.ds.append(u)
        self.ds = N.array(self.ds)

        self.hist = N.histogram(N.log(self.ds),bins=200,new=True)
        self.hist = [self.hist[0],self.hist[1],(self.hist[1][1]-self.hist[1][0])/2.]


################################################### END ALIGNMENT CLASS ###########################################

################################################### TEST METHODS  #################################################


##   alignment test program
##   runs on a small set of flight data and uses Crab, Geminga, and Vela as sources
##   if everything is ok, should return 0
def test():
    plr = os.environ['POINTLIKEROOT']
    fdir = plr+'/python/uw/utilities/boresighttest/'
    al = StackLoader(lis='cgv',files=['test'],datadr=fdir,ft2dr=fdir,srcdr=fdir)
    al.loadphotons(1,1000,2e6,0,999999999,-1)
    al.solverot()
    al.makeplot('aligntest.png')
    ret = 0
    if (al.params[0]+77.6158515)>10:
        ret = ret + 1
    if (al.params[1]+73.5283828322)>10:
        ret = ret + 2
    if (al.params[2]-67.346498744)>10:
        ret = ret + 4
    if (al.errors[0]-60.9043937028)>10:
        ret = ret+8
    if (al.errors[1]-59.0699490493)>10:
        ret = ret+16
    if (al.errors[2]-91.2188846332)>10:
        ret = ret +32
    return ret

##   psf parameter test program
##   runs on a small set of flight data and uses Crab, Geminga, and Vela as sources
##   if everything is ok, should return 0
def test2():
    plr = os.environ['POINTLIKEROOT']
    fdir = plr+'/python/uw/utilities/boresighttest/'
    os.system('cd %s'%fdir)
    al = StackLoader(lis='cgv',files=['test'],datadr=fdir,ft2dr=fdir,srcdr=fdir)
    al.loadphotons(10,1000,1770,0,999999999,0)
    al.solvepsf()
    ret=0
    if (al.sigma*rd-0.140)>1e-3:
        ret = ret + 1
    if (al.gamma-1.53)>1e-2:
        ret = ret + 2 
    al.makeplot('psftest.png')
    return ret
    
###################################################################################################