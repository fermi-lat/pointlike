"""
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/utilities/boresight.py,v 1.0 2010/07/29 13:53:17 mar0 Exp $
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

debug=False
rd = 180./N.pi
monthlist = ['aug2008','sep2008','oct2008','nov2008','dec2008','jan2009','feb2009','mar2009','apr2009'\
    ,'may2009','jun2009','jul2009','aug2009','sep2009','oct2009','nov2009','dec2009','jan2010','feb2010','mar2010'\
    ,'apr2010','may2010','jun2010']
CALDBdir = r'd:\fermi/CALDB/v1r1/CALDB/data/glast/lat/'
datadir = r'd:\fermi\data\flight/'
ft2dir = r'd:\fermi\data\flight/'
srcdir = r'd:\common\mar0\sourcelists/'
files = []
ft2s = []
for month in monthlist:
    files.append('%s'%month)
    ft2s.append('%s'%month)

###################################################  START STACKLOADER CLASS ##################################################    

##  Stackloader class
#   
#   loads a set of photons 
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
    def __init__(self,lis='strong',tev=False,rot=[0,0,0],files=files,CALDBdr=CALDBdir,datadr=datadir,ft2dr=ft2dir,srcdr=srcdir,test=False):
        
        self.firstlight = False
        self.rot = HepRotation(rot,False)
        self.rot.echo()
        self.tev = tev
        self.CALDBdir = CALDBdr
        self.datadir = datadr
        self.srcdir = srcdr
        self.ctmin=0.4
        self.ctmax=1.0

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
        s.IParams.init('P6_v3')
        self.atb = []
        self.srcs = []
        sf = file(self.srcdir+lis+'.txt')
        header = sf.readline()
        for lines in sf:
            line = lines.split()
            ra = float(line[1])
            dec = float(line[2])
            self.srcs.append(s.SkyDir(ra,dec))

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
        for ff in self.files:
            print ff
            tff = pf.open(ff)
            ttb = tff[1].data
            ttb = self.mask(ttb,self.srcs,emin,emax,start,stop,rad,cls)
            self.atb.append(ttb)
            tff.close()
        self.photons = []

        for j,tb in enumerate(self.atb):
            print 'Examining %d events'%len(tb)
            if len(tb)>0:
                print '    *Loading pointing history'
                phist = pl.PointingHistory(self.ft2s[j])
                print '    *Loading events'
                for k in range(len(tb)):
                    event = tb[k]
                    sd = s.SkyDir(float(event.field('RA')),float(event.field('DEC')))
                    rsrc = self.srcs[0]
                    diff = 1e40
                    for src in self.srcs:
                        tdiff = sd.difference(src)*rd
                        if tdiff<diff:
                            rsrc=src
                            diff=tdiff
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
                        if debug:
                            bpd.addPhoton(pl.Photon(sd,float(event.field('ENERGY')),float(time),int(event.field('EVENT_CLASS'))))
                del phist
        self.aeff = N.array(self.aeff).transpose()
        self.ebar = self.ebar/len(self.photons)

    #Find the uniform background component
    def solveback(self):
        alphas = [0.5]
        psfa = PSFAlign(lims=[0,self.rad/rd])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()

        #Model is a point source described 
        #by a psf in a uniform background
        cm.addModel(psfa)
        cm.addModel(bck)
        cm.fit(self.photons)
        Npsf = cm.minuit.params[0]
        Npsfe = cm.minuit.errors()[0][0]
        Nback = cm.minuit.params[1]
        Nbacke = cm.minuit.errors()[1][1]
        self.Npsf = Npsf
        self.Nback = Nback
        self.alpha = Npsf/(Npsf+Nback)
        self.sigalph = self.alpha*N.sqrt(Npsfe/(Npsf**2)+(Npsfe+Nbacke)/((Npsf+Nback)**2))
        self.il = cm.minuit.fval
        print 'Expected %d photons, got %d ( %d (%d) + %d (%d) )'%(len(self.photons),int(Npsf+Nback),int(Npsf),int(N.sqrt(Npsfe)),int(Nback),int(N.sqrt(Nbacke)))
        print self.alpha,self.sigalph
        self.cm=cm

    # Finds boresight alignment solution
    def solverot(self):
        psfa = PSFAlign(lims=[0,self.rad/rd],free=[True,True,True])
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


    def solvepsf(self):
        sigma=s.IParams.sigma(self.ebar,int(self.photons[0].event_class))
        gamma=s.IParams.gamma(self.ebar,int(self.photons[0].event_class))
        psf = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(bck)
        cm.fit(self.photons)
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

        self.il=cm.extlikelihood([self.Npsf,self.Nback,sigma,gamma],self.photons)

        if self.gamma>gamma:
            gsign = '+'
        else:
            gsign = '-'
        if self.sigma>sigma:
            ssign = '+'
        else:
            ssign = '-'
        print '**********************************************************'
        print 'Sigma = %1.3f [1 +/- %1.2f] (%s%1.0f)'%(self.sigma*rd,self.sigmae/self.sigma,ssign,100*abs((self.sigma-sigma)/sigma))
        print 'Gamma = %1.2f  [1 +/- %1.2f] (%s%1.0f)'%(self.gamma,self.gammae/self.gamma,gsign,100*abs((self.gamma-gamma)/gamma))
        TS = -2*(fl-self.il)
        print 'Significance of psf change was TS = %d'%TS
        self.cm=cm

    def solvehalo(self):
        sigma=s.IParams.sigma(self.ebar,int(self.photons[0].event_class))
        gamma=s.IParams.gamma(self.ebar,int(self.photons[0].event_class))
        psf = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma],free=[False,False])
        halo = Halo(lims=[0,self.rad/rd],model_par=[0.5/rd])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(halo)
        cm.addModel(bck)
        cm.fit(self.photons)
        print 'Halo width was %1.3f deg'%(cm.minuit.params[5]*rd)
        self.cm=cm

    ## Makes a plot of all of the models and the angular distribution of photons
    #  @param name filename of output file
    def makeplot(self,name='test.png'):
        self.getds()
        bins = 20.
        mi = N.log10((min(self.ds))**2)
        ma = N.log10((max(self.ds))**2)
        d2 = ma-mi
        be = N.arange(mi,ma,d2/bins)
        be2 = N.arange(mi,(1.-1./bins)*ma,d2/bins)
        be3 = N.sqrt(10**(be2))
        p.ioff()
        p.figure(figsize=(8,8))
        hist = p.hist(N.log10(self.ds*self.ds),bins=be2,fc='None')
        p.clf()
        fits=[]
        hists=[]
        errs=[]
        names=[]
        names.append('data')
        for model in range(len(self.cm.models)):
            fits.append([])
            names.append(self.cm.models[model].name)
        fits.append([])
        names.append('total')
        for it,ba in enumerate(be):
            left=be3[it]
            right=be3[it+1]
            for it2,model in enumerate(self.cm.models):
                fits[it2].append(self.cm.nest[it2]*model.integral(left,right)/model.integral(0,self.rad/rd)/(right**2-left**2)/(rd**2))
            fits[len(fits)-1].append(self.cm.integral(left,right,0,self.rad/rd)/(right**2-left**2)/(rd**2))
            hists.append(hist[0][it]/(right**2-left**2)/(rd**2))
            errs.append(N.sqrt(hist[0][it])/(right**2-left**2)/(rd**2))
        pts=[]

        p1 = p.errorbar(be+d2/bins/2+2*N.log10(rd),hists,yerr=errs,ls='None',marker='o')
        pts.append(p1[0])
        for ad in fits:
            p1 = p.plot(be+d2/bins/2+2*N.log10(rd),ad,'-')
            pts.append(p1)

        p.legend(pts,names)

        ax = p.gca()
        ax.set_yscale("log", nonposy='clip')
        p.ylim(0.1*min(hists),10*max(hists))
        p.xlim(min(be+2*N.log10(rd)),max(be+d2/bins+2*N.log10(rd)))
        p.xlabel(r'$log(\theta^2)\/[log(deg^2)]$')
        p.ylabel(r'$dN/d\Omega$')
        #p.ylim(0.5,max(hist[0])*2.)
        p.grid()
        p.savefig(name)
        
    ## tries to fit two King functions
    def solvedoublepsf(self):
        sigma=s.IParams.sigma(self.ebar,int(self.photons[0].event_class))
        gamma=s.IParams.gamma(self.ebar,int(self.photons[0].event_class))
        psf = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma])
        psf2 = PSF(lims=[0,self.rad/rd],model_par=[sigma,gamma])
        bck = Backg(lims=[0,self.rad/rd])
        cm = CompositeModel()
        cm.addModel(psf)
        cm.addModel(psf)
        cm.addModel(bck)
        cm.fit(self.photons)
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
        print '***************************************'
        print 'Sigma = %1.3f [1 +/- %1.2f] (%s%1.0f)'%(self.sigma*rd,self.sigmae/self.sigma,ssign,100*abs((self.sigma-sigma)/sigma))
        print 'Gamma = %1.2f  [1 +/- %1.2f] (%s%1.0f)'%(self.gamma,self.gammae/self.gamma,gsign,100*abs((self.gamma-gamma)/gamma))
        print 'Sigma2 = %1.3f [1 +/- %1.2f] (%s%1.0f)'%(self.sigma2*rd,self.sigmae2/self.sigma2,ssign2,100*abs((self.sigma2-sigma)/sigma))
        print 'Gamma2 = %1.2f  [1 +/- %1.2f] (%s%1.0f)'%(self.gamma2,self.gammae2/self.gamma2,gsign,100*abs((self.gamma2-gamma)/gamma))
        TS = -2*(fl-self.il)
        print 'Significance of psf change was TS = %d'%TS
        self.cm=cm

    #helper function - cuts unused data
    def mask(self,table,srcs,emin,emax,start,stop,rad,cls):
        cuts = [0,0,0,0,0,0]
        total = len(table)
        tc = len(table)
        msk = (table.field('TIME')>start) & (table.field('TIME')<stop)
        table = table[msk]
        tc = tc - len(table)
        cuts[0] = tc
        tc = len(table)
        if len(table)>0:
            msk = (table.field('ENERGY')>emin) & (table.field('ENERGY')<emax)
            table = table[msk]
            tc = tc - len(table)
            cuts[1] = tc
            tc = len(table)
            if len(table)>0:
                msk = (table.field('THETA')<(N.arccos(self.ctmin)*rd)) & (table.field('THETA')>(N.arccos(self.ctmax)*rd))
                table = table[msk]
                tc = tc - len(table)
                cuts[2] = tc
                tc = len(table)
                if len(table)>0:
                    msk = table.field('ZENITH_ANGLE')<105
                    table = table[msk]
                    tc = tc - len(table)
                    cuts[3] = tc
                    tc = len(table)
                    if len(table)>0:
                        if cls!=-1:
                            msk = table.field('CONVERSION_TYPE')==cls
                            table = table[msk]
                            tc = tc - len(table)
                            cuts[4] = tc
                            tc = len(table)
                        if len(table)>0:
                            msk = self.dirmask(srcs[0],rad,table)
                            for j in range(len(srcs)-1):
                                msk = msk | self.dirmask(srcs[j+1],rad,table)
                            table = table[msk]
                            tc = tc - len(table)
                            cuts[5] = tc
        tc = len(table)
        print 'TOTAL: %d    TIMEC: %d    ENERGY: %d    THETA: %d    ZENITH: %d    ECLASS: %d    POS: %d    EXAM: %d'%(total,cuts[0],cuts[1],cuts[2],cuts[3],cuts[4],cuts[5],tc)
        return table
    

    #direction mask - masks area around sources to speed execution
    def dirmask(self,sd,deg,tb):
        ra,dec=sd.ra(),sd.dec()
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
                            #ramin < 0 go to 360+ramin               racut                                  deccut
            mask = ( ((tb.field('RA')>360+ramin)&(ramin<0))  |  ((tb.field('RA')>ramin)&(tb.field('RA')<ramax))  )\
            & ( (tb.field('DEC')>decmin)&(tb.field('DEC')<decmax) )
        return mask

    #sets scaled angular deviations
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

########################################### TEST METHODS  #################################################

def test():
    plr = os.environ['POINTLIKEROOT']
    fdir = plr+'/python/uw/utilities/boresighttest/'
    al = StackLoader(lis='cgv',files=['test'],datadr=fdir,ft2dr=fdir,srcdr=fdir)
    al.loadphotons(1,1000,2e6,0,999999999,-1)
    al.solverot()
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
    al.makeplot()
    return ret
    
###################################################################################################