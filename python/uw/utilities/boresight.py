import numpy as N
import pylab as p
import skymaps as s
import pointlike as pl
import pyfits as pf
import glob as glob
import scipy.integrate as si
import os as os
from uw.utilities.minuit import Minuit

debug=False
monthlist = ['aug2008','sep2008','oct2008','nov2008','dec2008','jan2009','feb2009','mar2009','apr2009','may2009','jun2009','jul2009','aug2009','sep2009','oct2009','nov2009','dec2009','jan2010','feb2010','mar2010']
try:
    CALDBdir=os.environ['CALDB']
except:
    CALDBdir = r'd:\fermi/CALDB/v1r1/CALDB/data/glast/lat/'
datadir = r'd:\common\mar0\data/'#r'd:\fermi\data\flight/'
ft2dir = datadir#r'd:\fermi\data\flight/'
srcdir = r'd:\common\mar0\sourcelists/'
files = []
ft2s = []
for month in monthlist:
    files.append(datadir+'%s-ft1.fits'%month)
    ft2s.append(datadir+'%s-ft2.fits'%month)

############################################ START HEPROTATION CLASS ##########################################################

# Realization of the rotation matrix in CLHEP (takes either three angles, three vectors, or another matrix)
class HepRotation(object):

    def __init__(self,c,axes=True):
        if len(c)==1:
            self.matrix = c[0]
        if len(c)==3 and axes:
            vecx,vecy,vecz=c[0],c[1],c[2]
            self.matrix = N.matrix([[vecx.x(),vecy.x(),vecz.x()],[vecx.y(),vecy.y(),vecz.y()],[vecx.z(),vecy.z(),vecz.z()]])
        if len(c)==3 and not axes:
            x,y,z=c[0],c[1],c[2]
            #non-Euler
            self.matrix = N.matrix([[1,0,0],[0,N.cos(x),N.sin(x)],[0,-N.sin(x),N.cos(x)]])
            self.matrix = self.matrix * N.matrix([[N.cos(y),0,-N.sin(y)],[0,1,0],[N.sin(y),0,N.cos(y)]])
            self.matrix = self.matrix * N.matrix([[N.cos(z),N.sin(z),0],[-N.sin(z),N.cos(z),0],[0,0,1]])
            #Euler
            """self.matrix = N.matrix([[N.cos(x),N.sin(x),0],[-N.sin(x),N.cos(x),0],[0,0,1]])
            self.matrix = self.matrix * N.matrix([[1,0,0],[0,N.cos(y),N.sin(y)],[0,-N.sin(y),N.cos(y)]])
            self.matrix = self.matrix * N.matrix([[N.cos(z),N.sin(z),0],[-N.sin(z),N.cos(z),0],[0,0,1]])"""

    #print self
    def echo(self):
        for i in range(3):
                print self.matrix[i,0],self.matrix[i,1],self.matrix[i,2]

    #multiplication with either matrix or vector
    def m(self,other):
        if other.matrix.shape[1]==3:
            return HepRotation([self.matrix * other.matrix])
        else:
            return Hep3Vector([self.matrix * other.vec])

    #make inverse
    def inverse(self):
        return HepRotation([N.linalg.inv(self.matrix)])

################################################## END HEPROTATION CLASS ################################################


################################################## START HEP3VECTOR CLASS ###############################################

# Realization of the Vector class in CLHEP with SkyDir support

class Hep3Vector(object):

    def __init__(self,c):
        if len(c)==0:
            self.vec=c
        if len(c)==1:
            self.vec = c[0]
        if len(c)==2:
            phi,theta = c[0],c[1]
            self.vec = N.matrix([[N.cos(phi)*N.cos(theta)],[N.sin(phi)*N.cos(theta)],[N.sin(theta)]])
        if len(c)==3:
            x,y,z=c[0],c[1],c[2]
            self.vec = N.matrix([[x],[y],[z]])/N.sqrt(x*x+y*y+z*z)
        self.matrix =self.vec

    #convert from a skydir
    def __call__(self,sd):
        phi,theta = sd.ra()*N.pi/180,sd.dec()*N.pi/180
        self.vec = N.matrix([[N.cos(phi)*N.cos(theta)],[N.sin(phi)*N.cos(theta)],[N.sin(theta)]])

    def phi(self):
        x,y=self.x(),self.y()
        if x>0 and y>0:
            return N.arctan(y/x)
        if x<0 and y>0:
            return N.pi-N.arctan(y/x)
        if x<0 and y<0:
            return N.pi+N.arctan(y/x)
        if x>0 and y<0:
            return 2*N.pi-N.arctan(y/x)
        if x==0 and y>0:
            return N.pi/2
        if x==0 and y<0:
            return 3.*N.pi/2
        if y==0:
            return 0

    def dir(self):
        return s.SkyDir(self.ra(),self.dec())
    
    def ra(self):
        return self.phi()*180/N.pi
    
    def dec(self):
        return self.theta()*180/N.pi
    
    def theta(self):
        return N.arcsin(self.z())

    def x(self):
        return self.vec[0][0].item()
    
    def y(self):
        return self.vec[1][0].item()
    
    def z(self):
        return self.vec[2][0].item()

    def cross(self,other):
        return Hep3Vector([N.matrix([[self.y()*other.z()-self.z()*other.y()],[self.z()*other.x()-self.x()*other.z()],[self.x()*other.y()-self.y()*other.x()]])])

    def dot(self,other):
        return self.x()*other.x()+self.y()*other.y()+self.z()*other.z()

    def echo(self):
        for i in range(3):
            print self.vec[i,0]

################################################# END HEP3VECTOR CLASS #######################################################


################################################# START PHOTON CLASS #########################################################

# "Smart photon": associated with a source, and can be rotated in the GLAST frame

class Photon(object):

    def __init__(self,ra,dec,en,time,ec,hepx,hepz,srcdir):
        self.vec=Hep3Vector([ra*N.pi/180,dec*N.pi/180])
        self.energy=en
        self.time=time
        self.event_class=ec
        self.hepx=hepx
        self.hepz=hepz
        self.hepy=hepz.cross(hepx)
        self.rot=HepRotation([hepx,self.hepy,hepz])
        self.srcdir=srcdir

    # return rotated vector
    def rotate(self,rot):
        return self.rot.m(rot.m(self.rot.inverse().m(self.vec)))

    #return angular separation from source (after rotation in GLAST frame)
    def diff(self,rot):
        h1=Hep3Vector([0,0])
        h1(self.srcdir)
        return N.arccos(h1.dot(self.rotate(rot)))
    
    #return angular separation from source
    def srcdiff(self):
        h1=Hep3Vector([0,0])
        h1(self.srcdir)
        return N.arccos(h1.dot(self.vec))

###################################################  END PHOTON CLASS #######################################################


###################################################  START ALIGNMENT CLASS ##################################################    

class Alignment(object):


    #Constructor - sets up sources and pointing history file
    #lis = text file with all the source names, ra's, and dec's
    #tev = University of Washington flag for TEV machines
    #rot = optional rotation information input beforehand (in radians)
    #files = list of filenames in the form ['name1'] corresponding to 'name1-ft1.fits' and 'name1-ft2.fits'
    #CALDBdr = CALDBdir pointlike dir (should be '......./glast/lat/')
    #ft2dr = ft2 file directory
    #datadr = ft1 file directory
    #srcdr = source list directory (ascii file with 'name    ra     dec', first line is skipped)
    def __init__(self,lis='strong',tev=False,rot=[0,0,0],files=files,CALDBdr=CALDBdir,datadr=datadir,ft2dr=ft2dir,srcdr=srcdir,test=False):
        
        self.firstlight = False
        self.rot = HepRotation(rot,False)
        self.rot.echo()
        self.tev = tev
        self.CALDBdir = CALDBdr
        self.datadir = datadr
        self.srcdir = srcdr

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
        s.IParams.init('P6_v1')
        self.atb = []
        self.srcs = []
        sf = file(self.srcdir+lis+'.txt')
        header = sf.readline()
        for lines in sf:
            line = lines.split()
            ra = float(line[1])
            dec = float(line[2])
            self.srcs.append(s.SkyDir(ra,dec))

    #loads photons from FT1 file near sources  
    #rad = ROI in degrees
    #emin = minimum energy in MeV
    #emax = maximum energy in MeV
    #start = start time (MET)
    #stop = end time (MET)
    #cls = conversion type: 0=front,1=back,-1=all
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
        zens = []
        if debug:
            bpd = s.BinnedPhotonData(4)
        for j,tb in enumerate(self.atb):
            print 'Examining %d events'%len(tb)
            if len(tb)>0:
                print '    *Loading pointing history'
                phist = pl.PointingHistory(self.ft2s[j])
                print '    *Loading events'
                for k in range(len(tb)):
                    event = tb[k]
                    #print event.field('RA'),event.field('DEC')
                    sd = s.SkyDir(float(event.field('RA')),float(event.field('DEC')))
                    rsrc = self.srcs[0]
                    diff = 1e40
                    for src in self.srcs:
                        tdiff = sd.difference(src)*180/N.pi
                        if tdiff<diff:
                            rsrc=src
                            diff=tdiff
                    if diff<rad:
                        time = event.field('TIME')
                        pi = phist(time)
                        xax = pi.xAxis()
                        zax = pi.zAxis()
                        zen = pi.zenith()
                        zentheta = zax.difference(zen)*180/N.pi
                        #if abs(zentheta-35)<1 or abs(zentheta-50)<1 or abs(zentheta-39)<1:
                        zens.append(zentheta)
                        xv = Hep3Vector([])
                        zv = Hep3Vector([])
                        xv(xax)
                        zv(zax)
                        photon = Photon(sd.ra(),sd.dec(),event.field('ENERGY'),time,event.field('EVENT_CLASS'),xv,zv,rsrc)
                        self.ebar = self.ebar+event.field('ENERGY')
                        self.photons.append(photon)
                        #self.aeff.append([N.cos(event.field('THETA')*N.pi/180.),event.field('ENERGY')])
                        if debug:
                            bpd.addPhoton(pl.Photon(sd,float(event.field('ENERGY')),float(time),int(event.field('EVENT_CLASS'))))
                del phist
        #p.ioff()
        #p.hist(zens,bins=100,log=True)
        #p.grid()
        #p.ylabel('# in bin')
        #p.xlabel('Zenith Angle')
        #p.suptitle('Rocking angle of events in PSF/boresight analysis')
        #p.ion()
        #p.show()
        #for i in self.diffs:
        #    print i
        self.aeff = N.array(self.aeff).transpose()
        self.ebar = self.ebar/len(self.photons)
        if debug:
            bpd.write(r'd:\common\mar0\images\alignmentbpd.fits')

    #Find the uniform background component
    def solveback(self):
        alphas = [0.5]
        self.alphamin = Minuit(lambda x: self.loglike(x[0]),[alphas[0]])
        self.alphamin.minimize()
        self.il = self.alphamin.fval
        self.alpha = self.alphamin.params[0]
        #self.m = self.alphamin.params[1]
        self.sigalph = self.alphamin.errors()[0][0]
        print self.alpha,self.sigalph

    # Finds boresight alignment solution
    def solverot(self):
        steps=[20,20,30]
        params=[0,0,0]
        self.minuit = Minuit(lambda x: self.like(x[0],x[1],x[2],self.alpha),params,steps=steps)
        self.minuit.minimize()
        fl = self.minuit.fval
        TS = -2*(fl-self.il)
        print 'Significance of rotation was TS = %d'%TS
        errs = []
        self.errs = self.minuit.errors()
        for i in range(3):
            print self.minuit.params[i],N.sqrt(self.errs[i][i].item())
            errs.append(N.sqrt(self.errs[i][i].item()))
        self.params=self.minuit.params
        self.errs = errs

    #solve for psf parameters
    def solvepsf(self):
        self.getds()
        energy = self.ebar
        params=[s.IParams.sigma(float(energy),int(self.cls)),s.IParams.gamma(float(energy),int(self.cls))]   
        ta = [0.5]
        afmin = Minuit(lambda x: self.likesp(params[0],params[1],x),ta)
        afmin.minimize()
        alpha = afmin.params[0]
        self.psfmin=Minuit(lambda x: self.likesp(x[0],x[1],alpha),params)
        self.psfmin.minimize()
        self.sigerrs = self.psfmin.errors()
        self.sigerr = N.sqrt(self.sigerrs[0][0].item())
        self.gamerr = N.sqrt(self.sigerrs[1][1].item())
        self.cov = self.sigerrs[0][1].item()
        self.sigma = self.psfmin.params[0]
        self.gamma = self.psfmin.params[1]
        sigma = self.sigma
        gamma = self.gamma
        phots = len(self.photons)
        sigmamc = s.IParams.sigma(self.ebar,self.cls)
        gammamc = s.IParams.gamma(self.ebar,self.cls)
        p.ioff()
        #p.hold(True)
        dmax = self.rad*N.pi/180
        p.figure(figsize=(10,10))
        model = N.array([phots*self.integratemodel(N.exp(y),N.exp(y+2*self.hist[2]),alpha,sigma,gamma)/self.integratemodel(0,dmax,alpha,sigma,gamma) for y in self.hist[1]])
        mc = N.array([phots*self.integratemodel(N.exp(y),N.exp(y+2*self.hist[2]),alpha,sigmamc,gammamc)/self.integratemodel(0,dmax,alpha,sigmamc,gammamc) for y in self.hist[1]])
        back = N.array([phots*self.integrateback(N.exp(y),N.exp(y+2*self.hist[2]),alpha,sigma)/self.integrateback(0,dmax,alpha,sigma)*(1-alpha) for y in self.hist[1]])
        p.hist(N.log(self.ds),bins=200)
        p1=p.plot(self.hist[1]+self.hist[2],model,'g-',lw=4)
        p2=p.plot(self.hist[1]+self.hist[2],mc,'y-',lw=4)
        p3=p.plot(self.hist[1]+self.hist[2],back,'r-',lw=4)
        p.legend((p1,p2,p3),('opt','mc','back'))
        p.grid()
        p.ion()
        p.show()
        
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
                msk = table.field('THETA')<66.4
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
                cdec = N.cos((dec-deg)*N.pi/180.)
                ramin = ra-deg/cdec
                ramax = ra+deg/cdec
                decmin = dec-deg
                decmax = dec+deg
            else:
                cdec = N.cos((dec+deg)*N.pi/180.)
                ramin = ra-deg/cdec
                ramax = ra+deg/cdec
                decmin = dec-deg
                decmax = dec+deg
                            #ramin < 0 go to 360+ramin               racut                                  deccut
            mask = ( ((tb.field('RA')>360+ramin)&(ramin<0))  |  ((tb.field('RA')>ramin)&(tb.field('RA')<ramax))  )\
            & ( (tb.field('DEC')>decmin)&(tb.field('DEC')<decmax) )
        return mask


    #uniform background for calculating boresight alignment
    def like(self,x1,y1,z1,alpha):
        acc = 0
        x,y,z=x1*N.pi/180./3600.,y1*N.pi/180./3600.,z1*N.pi/180./3600.
        rot = HepRotation([x,y,z],False).m(self.rot)
        #rot.echo()
        for photon in self.photons:
            sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
            g = s.IParams.gamma(float(photon.energy),int(photon.event_class))
            umax = (self.rad*N.pi/sig/180)
            umax = umax*umax/2
            u = photon.diff(rot)/sig
            u = u*u/2
            fint = self.ipsf(umax,g)
            f = self.psf(u,g)
            logl = -N.log(alpha*f/fint+(1-alpha)/umax)
            acc = acc +logl
        #print '%d    %d   %d    %1.4f    %1.4f'%(x1,y1,z1,alpha,acc)
        return acc

    #uniform background for calculating background fraction
    def loglike(self,alpha):
        acc =0
        if alpha<0 or alpha>1:
            return 1e40
        for photon in self.photons:
            sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
            g = s.IParams.gamma(float(photon.energy),int(photon.event_class))
            umax = (self.rad*N.pi/sig/180)
            umax = umax*umax/2
            u = photon.srcdiff()/sig
            u = u*u/2
            fint = self.ipsf(umax,g)
            f = self.psf(u,g)
            logl = -N.log(alpha*f/fint+(1-alpha)/umax)
            acc = acc +logl
        #print '%1.4f    %1.4f'%(alpha,acc)
        return acc


    #linear model of background for background fraction
    def loglikep(self,a0,m):
        acc =0
        if a0<0 or a0>1:
            return 1e40
        for photon in self.photons:
            sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
            g = s.IParams.gamma(float(photon.energy),int(photon.event_class))
            umax = (self.rad*N.pi/sig/180)
            umax = umax*umax/2
            u = photon.srcdiff()/sig
            u = u*u/2
            fint = self.ipsf(umax,g)
            f = self.psf(u,g)
            alpha = self.alp(photon.energy,m,a0)
            logl = -N.log(alpha*f/fint+(1-alpha)/umax)
            acc = acc +logl
        #print '%1.4f    %1.4f    %1.4f'%(a0,m,acc)
        return acc

    #linear model for boresight alignment
    def likep(self,x1,y1,z1,a0,m):
        acc = 0
        x,y,z=x1*N.pi/180./3600.,y1*N.pi/180./3600.,z1*N.pi/180./3600.
        rot = HepRotation([x,y,z],False).m(self.rot)
        #rot.echo()
        for photon in self.photons:
            sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
            g = s.IParams.gamma(float(photon.energy),int(photon.event_class))
            umax = (self.rad*N.pi/sig/180)
            umax = umax*umax/2
            u = photon.diff(rot)/sig
            u = u*u/2
            fint = self.ipsf(umax,g)
            f = self.psf(u,g)
            alpha = self.alp(photon.energy,m,a0)
            logl = -N.log(alpha*f/fint+(1-alpha)/umax)
            acc = acc +logl
        #print '%d    %d   %d    %1.4f    %1.4f'%(x1,y1,z1,alpha,acc)
        return acc

    def likesp(self,sig,g,alpha):
        acc =0
        if sig<0 or g<1:
            return 1e40
        for j,num in enumerate(self.hist[0]):
            umax = (self.rad/sig*N.pi/180)
            umax = umax*umax/2
            u = N.exp(self.hist[1][j]+self.hist[2])/sig
            u = u*u/2
            fint = self.ipsf(umax,g)
            f = self.psf(u,g)
            logl = -N.log(alpha*f/fint+(1-alpha)/umax)+2*N.log(sig)
            acc = acc + num*logl
        #print '%1.4f    %1.4f    %1.4f'%(sig*57,g,acc)
        return acc

    #psf definition
    def psf(self,u,g):
        #g=100 #gaussian
        return (1-1/g)*(1+u/g)**(-g)

    #psf integral
    def ipsf(self,umax,g):
        #g=100 #gaussian
        return 1-(1+umax/g)**(-g+1)

    #returns scaled angular deviations
    def getus(self,x1,y1,z1):
        self.us=[]
        x,y,z=x1*N.pi/180./3600.,y1*N.pi/180./3600.,z1*N.pi/180./3600.
        rot = HepRotation([x,y,z],False)
        for photon in self.photons:
            sig = s.IParams.sigma(float(photon.energy),int(photon.event_class))
            umax = (self.rad*N.pi/sig/180)
            umax = umax*umax/2
            u = photon.diff(rot)/sig
            u = u*u/2
            #print photon.energy,sig,photon.srcdiff(),umax
            self.us.append(u)

    #returns angular deviations
    def getds(self):
        self.ds=[]
        for photon in self.photons:
            u = photon.diff(self.rot)
            #print photon.energy,sig,photon.srcdiff(),umax
            self.ds.append(u)
        self.ds = N.array(self.ds)

        self.hist = N.histogram(N.log(self.ds),bins=200,new=True)
        self.hist = [self.hist[0],self.hist[1],(self.hist[1][1]-self.hist[1][0])/2.]
        #print len(self.hist[0]),len(self.hist[1])

    #returns a vector of energies
    def getens(self):
        energy=[]
        for photon in self.photons:
            energy.append(photon.energy)
        return N.array(energy)

    #linear background model
    def alp(self,e,m,a0):
        alph=m*(e-self.emin)+a0
        if alph>1:
            return 1
        if alph<0:
            return 0
        return alph

    def integratemodel(self,dmin,dmax,alpha,sigma,gamma):
        umax = (self.rad/sigma*N.pi/180)
        umax = umax*umax/2
        fint = self.ipsf(umax,gamma)
        integral = si.quadrature(lambda x: alpha*self.psf(x,gamma)/fint+(1-alpha)/umax,dmin*dmin/(2*sigma*sigma),dmax*dmax/(2*sigma*sigma))
        return integral[0]

    def integrateback(self,dmin,dmax,alpha,sigma):
        umax = (self.rad/sigma*N.pi/180)
        umax = umax*umax/2
        integral = si.quadrature(lambda x: (1-alpha)/umax,dmin*dmin/(2*sigma*sigma),dmax*dmax/(2*sigma*sigma))
        return integral[0]

################################################### END ALIGNMENT CLASS ###########################################
  
class Opt(object):

    def __init__(self,tstart,tend,rad=1,emin=1000,emax=2e6,tev=False,lis='strong',r=[0,0,0],cls=-1):
        self.tstart=tstart
        self.tend=tend
        self.rad=rad
        self.emin=emin
        self.emax=emax
        self.tev=tev
        self.lis=lis
        self.r=r
        self.cls=cls

    def run(self,increment):
        self.al = Alignment(self.lis,self.tev,self.r)
        runs = (self.tend-self.tstart)/increment
        data = []
        for i in range(runs):
            st = self.tstart+i*increment
            sp = min(self.tend,self.tstart+(i+1)*increment)
            self.al.loadphotons(self.rad,self.emin,self.emax,st,sp,self.cls)
            self.al.solveback()
            self.al.solverot()
            data.append([st,sp,self.al.params[0],self.al.errs[0],self.al.params[1],self.al.errs[1],self.al.params[2],self.al.errs[2]])
        return N.array(data)

    def runone(self):
        return self.run(self.tend-self.tstart)[0]

    def iterate(self,n=3):
        rv = [0,0,0]
        for i in range(n):
            rvect = runone()
            rv = [rv[0]+rvect[2]*N.pi/180./3600.,rv[1]+rvect[4]*N.pi/180./3600.,rv[2]+rvect[6]*N.pi/180./3600.]
        print rv[0]*3600.*180./N.pi,rv[1]*3600.*180./N.pi,rv[2]*3600.*180./N.pi

    def runonepsf(self,cls):
        self.al = Alignment(self.lis,self.tev,self.r)
        self.al.loadphotons(self.rad,self.emin,self.emax,self.tstart,self.tend,cls)
        #al.solveback()
        self.al.solvepsf()
        return [self.al.sigma,self.al.sigerr,self.al.gamma,self.al.gamerr,self.al.cov]


def test():
    plr = os.environ['POINTLIKEROOT']
    fdir = plr+'/python/uw/utilities/boresighttest/'
    al = Alignment(lis='cgv',files=['test'],datadr=fdir,ft2dr=fdir,srcdr=fdir)
    al.loadphotons(1,1000,2e6,0,999999999,-1)
    al.solveback()
    al.solverot()
    ret = 0
    if (al.params[0]+77.6158515)>1:
        ret = ret + 1
    if (al.params[1]+73.5283828322)>1:
        ret = ret + 2
    if (al.params[2]-67.346498744)>1:
        ret = ret + 4
    if (al.errs[0]-60.9043937028)>1:
        ret = ret+8
    if (al.errs[1]-59.0699490493)>1:
        ret = ret+16
    if (al.errs[2]-91.2188846332)>1:
        ret = ret +32
    return ret
    
############################################# OLD CODE ###################################################
"""        if not True:
            self.alphap = so.fmin_powell(lambda x: self.loglike(x),0.5,ftol=1e-5,disp=1,full_output=1)
            self.bestp = so.fmin_powell(lambda x: self.like(x[0],x[1],x[2],self.alphap[0]),params,ftol=1e-5,disp=1,full_output=1)
            print self.alphap
            print self.bestp
            self.params = self.bestp[0]
            self.errs = steps
            mx = self.bestp[1].item()
            rg = N.arange(-50,60,10)
            tl = len(rg)
            self.xy = N.array([N.sqrt(self.like(x,y,self.params[2],self.alphap[0])-mx) for x in rg for y in rg]).reshape(tl,tl)
            self.xz = N.array([N.sqrt(self.like(x,self.params[1],y,self.alphap[0])-mx) for x in rg for y in rg]).reshape(tl,tl)
            self.yz = N.array([N.sqrt(self.like(self.params[0],x,y,self.alphap[0])-mx) for x in rg for y in rg]).reshape(tl,tl)
            p.figure(figsize=(10,10))
            p.subplot(2,2,1)
            p.contour(rg,rg,self.xy,[1,2,3,4,5])
            p.subplot(2,2,2)
            p.contour(rg,rg,self.xz,[1,2,3,4,5])
            p.subplot(2,2,3)
            p.contour(rg,rg,self.yz,[1,2,3,4,5])
            self.getus(0,0,0)
            p.ioff()
            p.figure(figsize=(6,6))
            bins = int((1-self.alpha)*len(self.photons))/20
            his = p.hist(self.us,bins=bins,fc='none')
            p.grid()
            p.xlabel(r'$\theta^2/2\sigma^2$')
            p.ylabel(r'# in bin')
            p.plot([0,max(his[1])],[20,20],'r-')
            p.ion()
            p.show()"""