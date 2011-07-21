import skymaps as sk
from uw.utilities.minuit import Minuit
import uw.stacklike.stacklike as b
import numpy as np
import scipy.optimize as so
import glob as glob
import uw.like.pypsf as up
from uw.like.pycaldb import CALDBManager
import uw.stacklike.angularmodels as ua
import pyfits as pf
import pylab as py
import numpy as N
import copy as cp
import matplotlib
import os as os
import time as t
from uw.stacklike.CLHEP import Photon

rd = 180/np.pi   #DEG2RAD
trad = rd
ft1dir = r'/phys/groups/tev/scratch1/users/Fermi/mar0/data/6.3/'          #ft1 fits files directory
stacklist = 'cgv'                    #text file containing 'name ra dec' for each stacked source
cdb = r'/phys/groups/tev/scratch1/users/Fermi/CALDB/v1r1/CALDB/data/glast/lat'     #CALDB directory
sdir='/phys/users/mar0/sourcelists/'

## Main method
# @param ec Conversion Type (0-front 1-back)
def main(ec):
    f = Fitter(ec,ft1dir,stacklist)


################################## BEGIN FITTER CLASS ####################################

## Fitter class
#  attempts to fit a smooth PSF as a function of energy
#  by fitting the R68 and R95 containment
#  
#  Fit function:
#  multiple scattering: R68 ~ E^-0.8
#  instrument pitch: R68 ~ constant
#  single scattering tails: R95 ~ E^0.1 above 10 GeV
#
#  R68/R95^2 = c0^2*(E/100)^(-2*b)+c1^2+c2^2*(E/100)^(-2*b1)
class Fitter(object):

    ## Constructor
    #  runs optimization
    #  energy bins (GeV): [1.0-1.7],[1.7-3.2],[3.2-5.6],[5.6-10],[10-17],[17-32],[32-100]
    def __init__(self,ec,ftdir,slist):
        self.lastlike=0
        self.lastbin=[]
        self.iterations=0

        #set up the initial condition based on P6v3
        if ec==0:
            pars = [4.0,0.11,-0.84,0.25,-0.05,11.,0.3,-0.84,0.6,-0.05]
        else:
            pars = [6.7,0.26,-0.79,0.6,-0.05,21.,0.76,-0.79,1.,-0.05]

        
        steps = abs(np.array(pars))*0.2
        lims = [[0.01,20.],[0.00001,1.],[-2.,-0.001],[0.01,40.],[0.00001,1.],[-2.,-0.001]]

        emin=[]
        emax=[]
        irf = 'P6_v3_diff'
        st1=0.25
        st2=0.25
        st3=0.5
        sc1 = np.arange(3.,4.,st1)
        sc2 = np.arange(4.,4.5,st2)
        sc3 = np.arange(4.5,5.,st3)

        for i in sc1:
            emin.append(10**(i))
            emax.append(10**(i+st1))
        for i in sc2:
            emin.append(10**(i))
            emax.append(10**(i+st2))
        for i in sc3:
            emin.append(10**(i))
            emax.append(10**(i+st3))
        
        altb=[]

        files = glob.glob(ftdir+'*-ft1.fits')

        ## setup the stacker for each energy bin
        for it,en in enumerate(emin):
            al = b.StackLoader(lis=slist,quiet=False,useft2s=False,irf=irf,ctmin=0.4,srcdir=sdir,ft1dir=ft1dir)
            al.files=files
            al.loadphotons(0.,4.,en,emax[it],0,999999999,ec)
            events = len(al.photons)
            al.bindata()
            al.solveback()
            altb.append(al)
            #al.makeplot(name='/phys/users/mar0/figures/%s%d%d'%(slist,al.ebar,al.cls),scale='linear',bin=len(hist[0]))

        ## optimize curves
        almin = so.fmin_powell(lambda x:self.likelihood(altb,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]),pars,disp=1,full_output=1)
        
        ## print out fit parameters
        print almin[0]


    ## likelihood function
    #  make a king function from the 68 and 95% containment calculated from the smooth
    #  function for each energy bin
    #
    #  fit a combination of point-source and isotropic background
    #  and return likelihood
    
    def likelihood(self,altb,c068,c168,b68,c268,b168,c095,c195,b95,c295,b195):
        print c068,c168,b68,c268,b168,c095,c195,b95,c295,b195

        #catch garbage values
        if np.isnan(c068) or np.isnan(c095):
            return 1e40
        acc=0.
        print '######################################'
        self.iterations=self.iterations+1
        print 'iteration: %d'%self.iterations
        
        ## calculate the likelihood of the parameterization for each energy bin
        for it,tb in enumerate(altb):
            eb = tb.ebar
            r68 = self.scale(eb,c068,c168,b68,c268,b168)/rd
            r95 = self.scale(eb,c095,c195,b95,c295,b195)/rd

            #R68 should always be less than R95!
            if r95<=r68:
                like = -1.e40
            else:
                if c068<(10.*c268) or c095<(10.*c295):
                    like=-1.e40
                else:
                    like = tb.TS(ct=[r68,r95])/2.
            if self.iterations==1:
                self.lastbin.append(like)
            acc = acc-like
            print 'Energy: %1.0f  R68: %1.2f  R95: %1.2f  like: %1.0f  diff= %1.0f'%(eb,r68*rd,r95*rd,like,self.lastbin[it]-like)
            self.lastbin[it]=like
        print '--------------------------------------'
        print '68  c0: %1.2f  c1: %1.2f  b: %1.2f  c2: %1.2f  b1: %1.2f'%(c068,c168,b68,c268,b168)
        print '95  c0: %1.2f  c1: %1.2f  b: %1.2f  c2: %1.2f  b1: %1.2f'%(c095,c195,b95,c295,b195)
        print ' Total: %1.2f      Diff: %1.2f'%(acc,acc-self.lastlike)
        print '######################################'
        self.lastlike=acc
        return acc

    ## function that describes the smoothed R68 and R95
    def scale(self,e,c0,c1,b,c2,b1):
        return np.sqrt((c0*((e/100.)**b))**2+c1**2+(c2*((e/100.)**b1))**2)

############################## END FITTER CLASS ###################################

## calculates the chisq for a 3 parameter scale function and the sigma parameter
#  helps set up the scaling function parameters
#  @param el energy list
#  @param er1 sigma error list
#  @param c0 PSF sigma @ 100 MeV
#  @param c1 PSF sigma @ 100 GeV
#  @param b power law of multiple scattering
def chisq(sigmal,el,erl,c0,c1,b):
    if b>0:
        return 1e40
    return sum([(sigmal[x]-scale2(el[x],c0,c1,b))**2/(erl[x]**2) for x in range(len(el))])
   

## calculates the chisq for a both front and parameter scale function and the sigma parameter
#  helps set up the scaling function paramters
#  @param el energy list
#  @param er1 sigma error list
#  @param c0 PSF sigma @ 100 MeV
#  @param c1 instrument pitch
#  @param b power law of multiple scattering
#  @param 
def chisq2(sigmaf,sigmab,el,erf,erb,c0f,c1f,c0b,c1b,b):
    if b>0:
        return 1e40
    acc = 0
    acc = acc + chisq(sigmaf,el,erf,c0f,c1f,b)
    acc = acc + chisq(sigmab,el,erb,c0b,c1b,b)
    return acc

def scale(e,cl):
    if cl==0:
        p0,p1=58e-3,377e-6
    else:
        p0,p1=96e-3,1300e-6

    return N.sqrt((p0*(e/100)**-0.8)**2+p1**2)

## 5 parameter scale function
#  @param e energy in MeV
#  @param c0 sigma at 100 MeV
#  @param c1 sigma from instrument pitch
#  @param b power law of multiple scattering
#  @param c2 sigma from single scattering tail
#  @param b1 power law of single scattering tail
def scale3(e,c0,c1,b,c2,b1):
    return N.sqrt((c0*((e/100.)**b))**2+c1**2+(c2*((e/100.)**b1))**2)

## 3 parameter scale function
#  @param e energy in MeV
#  @param c0 sigma at 100 MeV
#  @param c1 sigma from instrument pitch
#  @param b power law of multiple scattering
def scale2(e,c0,c1,b):
    return N.sqrt((c0*((e/100.)**b))**2+c1**2)

## makes a CALDB PSF output from the uncommented parameter sets 
#  for 5 parameter containment function
#  
#  calculates new scaling function to make PSF values between bins
#  continuous
#
#  r68f - parameterization for R68 front
#  r95f - parameterization for R95 front
#  r68b - parameterization for R68 back
#  r95b - parameterization for R95 back
#
#  some common ones are listed below
def getirfparams(irf='P6_v8_diff',out='P6_v10_diff'):
    bpsf = ua.PSF(lims=[0,10],model_par=[0.001,2.25])
    fpsf = ua.PSF(lims=[0,10],model_par=[0.001,2.25])

    ######## P6_v11 3 pars ###############
    #r68fp = [4.022,0.1086,-0.842]
    #r95fp = [10.912,0.296,-0.8405]
    #r68bp = [6.7,0.2632,-0.796]
    #r95bp = [20.997,0.7558,-0.786]

    ######## P6_v11 5 pars ###############
    #r68fp = [4.022,0.0289,-0.842,0.1628,-0.0906]
    #r95fp = [10.912,0.0841,-0.8405,0.3695,-0.0725]
    #r68bp = [6.7,0.0961,-0.796,0.4235,-0.1016]
    #r95bp = [20.997,0.4271,-0.786,1.021,-0.0663]

    ######## P6_v12 5 pars ###############
    r68fp = [3.639,-0.006,-0.8354,0.1943,-0.0786]
    r95fp = [8.297,0.0288,-0.7702,0.4022,-0.1096]
    r68bp = [7.411,0.045,-0.861,0.4557,-0.0927]
    r95bp = [17.956,0.3554,-0.7484,0.9265,-0.0624]

    ######## P7_v3 5 pars ###############
    #r68fp = [3.701,0.0289,-0.8256,0.1628,-0.0906]
    #r95fp = [9.635,0.0841,-0.8181,0.3695,-0.0725]
    #r68bp = [7.589,0.0961,-0.8597,0.4235,-0.1016]
    #r95bp = [18.59,0.4271,-0.7649,1.021,-0.0663]

    ######## P7_v3 ###############
    #r68fp = [3.919,0.1241,-0.848]
    #r95fp = [10.954,0.329,-0.845]
    #r68bp = [6.77,0.2853,-0.803]
    #r95bp = [19.983,0.8327,-0.784]

    ######## P7_v4 uc 5 pars ###############
    #r68fp = [4.004,0.0000,-0.8106,0.2035,-0.0731]
    #r95fp = [13.01,0.0029,-0.8549,0.5428,-0.0405]
    #r68bp = [7.664,0.1361,-0.8396,0.4373,-0.0792]
    #r95bp = [26.22,0.6692,-0.8148,1.025,-0.0596]

    ######## P7_v4 src 5 pars ###############
    #r68fp = [4.155,0.0161,-0.8160,0.2156,-0.0739]
    #r95fp = [13.01,0.0407,-0.8389,0.5099,0.0165]
    #r68bp = [7.697,0.0793,-0.8401,0.4721,-0.0856]
    #r95bp = [24.82,0.3825,-0.8018,0.9080,0.0075]

    ######## P7_v4 clean 5 pars ###############
    #r68fp = [4.057,0.0099,-0.8145,0.2090,-0.0793]
    #r95fp = [12.17,0.0613,-0.8285,0.4340,0.0034]
    #r68bp = [7.575,0.1139,-0.8332,0.4637,-0.0784]
    #r95bp = [25.52,0.3097,-0.7823,0.8766,0.0470]

    ######## P6_v12 5 pars ###############
    r68fp = [3.525,-0.002,-0.8261,0.1997,-0.0825]
    r95fp = [7.846,0.0263,-0.7468,0.3948,-0.1152]
    r68bp = [5.905,0.037,-0.8749,0.4633,-0.0903]
    r95bp = [12.423,0.2504,-0.81,0.9353,-0.0510]

    ######## P6_v12 5 pars cgv ###############
    #r68fp = [3.477,-0.0000,-0.8271,0.2207,-0.0525]
    #r95fp = [9.662,0.0126,-0.8254,0.5893,-0.0962]
    #r68bp = [6.658,0.0334,-0.8516,0.4959,-0.0642]
    #r95bp = [16.172,1.003,-0.8369,1.3854,-0.0490]

    pname = cdb+r'\bcf\psf\psf_%s_'%irf
    print pname
    cwd = os.getcwd()
    cmd1 = r'copy %sfront.fits %s\\psf_%s_front.fits'%(pname,cwd,out)
    cmd2 = r'copy %sback.fits %s\\psf_%s_back.fits'%(pname,cwd,out)
    print cmd1
    print cmd2
    os.system(cmd1)
    os.system(cmd2)
    ff = pf.open('psf_%s_front.fits'%out,mode='update')
    bf = pf.open('psf_%s_back.fits'%out,mode='update')

    cth = ff[1].data.field('CTHETA_LO')[0]
    en = ff[1].data.field('ENERG_LO')[0]
    ct = len(ff[1].data.field('CTHETA_LO')[0])
    e = len(ff[1].data.field('ENERG_LO')[0])

    ens = []
    fss =[]
    bss =[]

    for i in N.arange(0,e,1):
        energy = N.sqrt(ff[1].data.field('ENERG_LO')[0][i]*ff[1].data.field('ENERG_HI')[0][i])
        ens.append(energy)
        """r68f = scale2(energy,r68fp[0],r68fp[1],r68fp[2])
        r95f = scale2(energy,r95fp[0],r95fp[1],r95fp[2])
        r68b = scale2(energy,r68bp[0],r68bp[1],r68bp[2])
        r95b = scale2(energy,r95bp[0],r95bp[1],r95bp[2])"""

        r68f = scale3(energy,r68fp[0],r68fp[1],r68fp[2],r68fp[3],r68fp[4])
        r95f = scale3(energy,r95fp[0],r95fp[1],r95fp[2],r95fp[3],r95fp[4])
        r68b = scale3(energy,r68bp[0],r68bp[1],r68bp[2],r68bp[3],r68bp[4])
        r95b = scale3(energy,r95bp[0],r95bp[1],r95bp[2],r95bp[3],r95bp[4])


        fpsf.fromcontain(1.,r95f/r68f,0.68,0.95)
        bpsf.fromcontain(1.,r95b/r68b,0.68,0.95)

        print r68f,r95f
        print r68b,r95b

        fsigma = fpsf.model_par[0]*r68f
        bsigma = bpsf.model_par[0]*r68b
        fgamma = fpsf.model_par[1]
        bgamma = bpsf.model_par[1]
        fss.append(fsigma*trad)
        bss.append(bsigma*trad)

        scf = scale(energy,0)
        scb = scale(energy,1)

        for j in range(ct):
            try:
                print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,ff[1].data.field('SIGMA')[0][e*j+i],fsigma*trad/scf)
                ff[1].data.field('SIGMA')[0][e*j+i]=fsigma*trad/scf
                print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,bf[1].data.field('SIGMA')[0][e*j+i],bsigma*trad/scb)
                bf[1].data.field('SIGMA')[0][e*j+i]=bsigma*trad/scb
            except:
                print 'New style psf'
            try:
                print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,ff[1].data.field('SCORE')[0][e*j+i],fsigma*trad/scf)
                ff[1].data.field('SCORE')[0][e*j+i]=fsigma*trad/scf
                print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,bf[1].data.field('SCORE')[0][e*j+i],bsigma*trad/scb)
                bf[1].data.field('SCORE')[0][e*j+i]=bsigma*trad/scb
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,ff[1].data.field('STAIL')[0][e*j+i],fsigma*trad/scf)
                ff[1].data.field('STAIL')[0][e*j+i]=fsigma*trad/scf
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,bf[1].data.field('STAIL')[0][e*j+i],bsigma*trad/scb)
                bf[1].data.field('STAIL')[0][e*j+i]=bsigma*trad/scb
            except:
                print 'Old style psf'
            print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,ff[1].data.field('NTAIL')[0][e*j+i],0)
            ff[1].data.field('NTAIL')[0][e*j+i]=0
            print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,bf[1].data.field('NTAIL')[0][e*j+i],0)
            bf[1].data.field('NTAIL')[0][e*j+i]=0
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,ff[1].data.field('GCORE')[0][e*j+i],fgamma)
            ff[1].data.field('GCORE')[0][e*j+i]=fgamma
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,bf[1].data.field('GCORE')[0][e*j+i],bgamma)
            bf[1].data.field('GCORE')[0][e*j+i]=bgamma
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,ff[1].data.field('GTAIL')[0][e*j+i],fgamma)
            ff[1].data.field('GTAIL')[0][e*j+i]=fgamma
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,bf[1].data.field('GTAIL')[0][e*j+i],bgamma)
            bf[1].data.field('GTAIL')[0][e*j+i]=bgamma
        #t.sleep(2)

    fss = N.array(fss)
    fse = 0.10*fss
    bss = N.array(bss)
    bse = 0.10*bss

    pars = [58e-3,377e-6,96e-3,1300e-6,-0.8]
    best = Minuit(lambda x: chisq2(fss,bss,ens,fse,bse,x[0],x[1],x[2],x[3],x[4]),params=pars,up=1,printMode=0,strategy=2,tolerance=1e-12)
    best.minimize()

    for it,energy in enumerate(ens):
        print energy,fss[it],scale2(energy,best.params[0],best.params[1],best.params[4]),bss[it],scale2(energy,best.params[2],best.params[3],best.params[4])

    for it,par in enumerate(best.params):
        ff[2].data.field('PSFSCALE')[0][it]=par
        bf[2].data.field('PSFSCALE')[0][it]=par

    for i in N.arange(0,e,1):
        energy = N.sqrt(ff[1].data.field('ENERG_LO')[0][i]*ff[1].data.field('ENERG_HI')[0][i])
        """r68f = scale2(energy,r68fp[0],r68fp[1],r68fp[2])
        r95f = scale2(energy,r95fp[0],r95fp[1],r95fp[2])
        r68b = scale2(energy,r68bp[0],r68bp[1],r68bp[2])
        r95b = scale2(energy,r95bp[0],r95bp[1],r95bp[2])"""

        r68f = scale3(energy,r68fp[0],r68fp[1],r68fp[2],r68fp[3],r68fp[4])
        r95f = scale3(energy,r95fp[0],r95fp[1],r95fp[2],r95fp[3],r95fp[4])
        r68b = scale3(energy,r68bp[0],r68bp[1],r68bp[2],r68bp[3],r68bp[4])
        r95b = scale3(energy,r95bp[0],r95bp[1],r95bp[2],r95bp[3],r95bp[4])

        fpsf.fromcontain(1.,r95f/r68f,0.68,0.95)
        bpsf.fromcontain(1.,r95b/r68b,0.68,0.95)

        print r68f,r95f
        print r68b,r95b

        fsigma = fpsf.model_par[0]*r68f
        bsigma = bpsf.model_par[0]*r68b
        fgamma = fpsf.model_par[1]
        bgamma = bpsf.model_par[1]

        scf = scale2(energy,best.params[0],best.params[1],best.params[4])
        scb = scale2(energy,best.params[2],best.params[3],best.params[4])

        for j in range(ct):
            try:
                print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,ff[1].data.field('SIGMA')[0][e*j+i],fsigma*trad/scf)
                ff[1].data.field('SIGMA')[0][e*j+i]=fsigma*trad/scf
                print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,bf[1].data.field('SIGMA')[0][e*j+i],bsigma*trad/scb)
                bf[1].data.field('SIGMA')[0][e*j+i]=bsigma*trad/scb
            except:
                print 'New style psf'
            try:
                print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,ff[1].data.field('SCORE')[0][e*j+i],fsigma*trad/scf)
                ff[1].data.field('SCORE')[0][e*j+i]=fsigma*trad/scf
                print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,bf[1].data.field('SCORE')[0][e*j+i],bsigma*trad/scb)
                bf[1].data.field('SCORE')[0][e*j+i]=bsigma*trad/scb
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,ff[1].data.field('STAIL')[0][e*j+i],fsigma*trad/scf)
                ff[1].data.field('STAIL')[0][e*j+i]=fsigma*trad/scf
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,bf[1].data.field('STAIL')[0][e*j+i],bsigma*trad/scb)
                bf[1].data.field('STAIL')[0][e*j+i]=bsigma*trad/scb
            except:
                print 'Old style psf'
            print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,ff[1].data.field('NTAIL')[0][e*j+i],0)
            ff[1].data.field('NTAIL')[0][e*j+i]=0
            print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,bf[1].data.field('NTAIL')[0][e*j+i],0)
            bf[1].data.field('NTAIL')[0][e*j+i]=0
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,ff[1].data.field('GCORE')[0][e*j+i],fgamma)
            ff[1].data.field('GCORE')[0][e*j+i]=fgamma
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,bf[1].data.field('GCORE')[0][e*j+i],bgamma)
            bf[1].data.field('GCORE')[0][e*j+i]=bgamma
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,ff[1].data.field('GTAIL')[0][e*j+i],fgamma)
            ff[1].data.field('GTAIL')[0][e*j+i]=fgamma
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,bf[1].data.field('GTAIL')[0][e*j+i],bgamma)
            bf[1].data.field('GTAIL')[0][e*j+i]=bgamma


    ff.flush()
    bf.flush()

def makeplots(irfs=['P6_v11_diff'],fact=1.,num=32.,cts=[34,68,95]):
    
    nms = [str(x)+'% '+y for x in cts for y in irfs]
    clr = ['b','k','g','r']
    #sty = ['-','--','-.',':']
    sty = ['solid','dashed','dashdot','dotted']

    energy = N.arange(1.75,5.75,1./num)
    energy=10**energy
    energy_r = np.array([energy[-i-1] for i in range(len(energy))])
    py.ioff()
    py.hold(True)
    py.figure(figsize=(10.5,4.75))
    py.subplot(1,2,1)
    py.loglog()
    #print energy
    ps = []
    for it1,ct in enumerate(cts):
        for it2,irf in enumerate(irfs):
            psf = up.CALDBPsf(CALDBManager(irf=irf))
            rct = [np.sqrt(psf.inverse_integral(x/fact,0,ct,ctbin=0))*np.sqrt(psf.inverse_integral(x*fact,0,ct,ctbin=0)) for x in energy]
            rct60 = [np.sqrt(psf.inverse_integral(x/fact,0,ct,ctbin=5))*np.sqrt(psf.inverse_integral(x*fact,0,ct,ctbin=5)) for x in energy_r]
            p1 = py.fill(np.hstack([energy,energy_r]),np.hstack([rct,rct60]),color=clr[it2])#,linestyle=sty[it2])#,alpha=0.4)
            p1[0].set_alpha(0.4)
            #print rct
            ps.append(p1)
    py.grid()
    py.ylim(1e-2,1e2)
    prop = matplotlib.font_manager.FontProperties(size=9) 
    py.legend(ps,nms,bbox_to_anchor=(1.0, 1.0),prop=prop)
    py.xlabel('Energy (MeV)',fontsize=9)
    py.ylabel('PSF containment(deg)',fontsize=9)
    py.title('Front')
    py.subplot(1,2,2)
    py.loglog()
    ps = []
    for it1,ct in enumerate(cts):
        for it2,irf in enumerate(irfs):
            psf = up.CALDBPsf(CALDBManager(irf=irf))
            rct = [np.sqrt(psf.inverse_integral(x/fact,1,ct,ctbin=0))*np.sqrt(psf.inverse_integral(x*fact,1,ct,ctbin=0)) for x in energy]
            rct60 = [np.sqrt(psf.inverse_integral(x/fact,1,ct,ctbin=5))*np.sqrt(psf.inverse_integral(x*fact,1,ct,ctbin=5)) for x in energy_r]
            p1 = py.fill(np.hstack([energy,energy_r]),np.hstack([rct,rct60]),color=clr[it2])#,linestyle=sty[it2])#,alpha=0.4)
            p1[0].set_alpha(0.4)
            ps.append(p1)
    py.grid()
    py.ylim(1e-2,1e2)
    prop = matplotlib.font_manager.FontProperties(size=9) 
    py.legend(ps,nms,bbox_to_anchor=(1.0, 1.0),prop=prop)
    py.xlabel('Energy (MeV)',fontsize=9)
    py.ylabel('PSF containment(deg)',fontsize=9)
    py.title('Back')
    title = ''
    it = 0
    for irf in irfs:
        title=title+irf
        if it<(len(irfs)-1):
            title = title+'__'
        it = it +1
    title2=(title.replace('__',', ')).replace('_',' ')
    py.suptitle(title2)
    py.savefig(title.strip('_')+'.png')

def newmethod(fname,irf):
    out = 'P7SOURCE_V11'
    ff1 = open(fname+'0.txt')
    bf1 = open(fname+'1.txt')

    offset = 3
    copyabove = 13

    febar = []
    fsig = []
    fgam = []
    fms = []
    fmg = []
    bebar = []
    bsig = []
    bgam = []
    bms = []
    bmg = []
    for lines in ff1:
        line = lines.split()
        febar.append(float(line[0]))
        fsig.append(float(line[1])*rd)
        fgam.append(float(line[2]))
        fms.append(float(line[3]))
        fmg.append(float(line[4]))
    for lines in bf1:
        line = lines.split()
        bebar.append(float(line[0]))
        bsig.append(float(line[1])*rd)
        bgam.append(float(line[2]))
        bms.append(float(line[3]))
        bmg.append(float(line[4]))
    ff1.close()
    bf1.close()
    fgam = np.array(fgam)
    fsig = np.array(fsig)/trad
    fsig2 = fsig/np.sqrt(1-1./fgam)
    bgam = np.array(bgam)
    bsig = np.array(bsig)/trad
    bsig2 = bsig/np.sqrt(1-1./bgam)
    pars = [1.8/trad,0.065/trad,1.8/trad,0.065/trad,-0.8]
    minu = Minuit(lambda x: chisq2(fsig2,bsig2,febar,np.sqrt(fsig2),np.sqrt(bsig2),x[0],x[1],x[2],x[3],x[4]),pars)
    minu.minimize()
    pars = minu.params
    print pars

    pname = cdb+r'/bcf/psf/psf_%s_'%irf
    print pname
    cwd = os.getcwd()
    cmd1 = r'cp %sfront.fits %s/psf_%s_front.fits'%(pname,cwd,out)
    cmd2 = r'cp %sback.fits %s/psf_%s_back.fits'%(pname,cwd,out)
    print cmd1
    print cmd2
    os.system(cmd1)
    os.system(cmd2)
    ff = pf.open('psf_%s_front.fits'%out,mode='update')
    bf = pf.open('psf_%s_back.fits'%out,mode='update')

    cth = ff[1].data.field('CTHETA_LO')[0]
    cth2 = ff[1].data.field('CTHETA_HI')[0]
    en = ff[1].data.field('ENERG_LO')[0]
    ct = len(ff[1].data.field('CTHETA_LO')[0])
    e = len(ff[1].data.field('ENERG_LO')[0])

    psfsc = ff[2].data
    for it,par in enumerate(pars):
        psfsc.field('PSFSCALE')[0][it]=par
    ff[2].data = psfsc
    bf[2].data = psfsc

    ftb = ff[1].data
    btb = bf[1].data

    for i in range(e-offset):
        for j in range(ct):
            ctbar = 0.5*(cth[j]+cth2[j])
            idx = e*j+i+offset
            ebar = np.sqrt(ff[1].data.field('ENERG_LO')[0][i+offset]*ff[1].data.field('ENERG_HI')[0][i+offset])
            scf = scale2(ebar,pars[0],pars[1],pars[4])
            scb = scale2(ebar,pars[2],pars[3],pars[4])
            
            if (i+offset)>copyabove:
                it = len(fsig)-1
            else:
                it = i

            fsfact=((1.-ctbar)*fms[it]+1.)
            fgfact=((1.-ctbar)*fmg[it]+1.)
            bsfact=((1.-ctbar)*bms[it]+1.)
            bgfact=((1.-ctbar)*bmg[it]+1.)
            ftb.field('SCORE')[0][idx]=fsig[it]/scf*fsfact
            #ftb.field('SIGMA')[0][idx]=fsig[it]/scf*fsfact
            ftb.field('STAIL')[0][idx]=fsig[it]/scf*fsfact
            ftb.field('GCORE')[0][idx]=fgam[it]*fgfact
            ftb.field('GTAIL')[0][idx]=fgam[it]*fgfact
            ftb.field('NTAIL')[0][idx]=0

            btb.field('SCORE')[0][idx]=bsig[it]/scb*bsfact
            btb.field('STAIL')[0][idx]=bsig[it]/scb*bsfact
            #btb.field('SIGMA')[0][idx]=bsig[it]/scb*bsfact
            btb.field('GCORE')[0][idx]=bgam[it]*bgfact
            btb.field('GTAIL')[0][idx]=bgam[it]*bgfact
            btb.field('NTAIL')[0][idx]=0

    ff[1].data = ftb
    bf[1].data = btb
    ff.flush()
    bf.flush()

def mergeresults(resdir):
    flist = np.sort(glob.glob('%s/*[0-9]_0.txt'%resdir))
    blist = np.sort(glob.glob('%s/*[0-9]_1.txt'%resdir))
    fout = open('%s/summary_0.txt'%resdir,'w')
    bout = open('%s/summary_1.txt'%resdir,'w')
    ebars = []
    sigs = []
    gams = []
    ms = []
    mg = []
    tbars = []
    for fl in flist:
        data = fl.split('/')[-1].split('_')
        ebar = np.sqrt(int(data[0])*int(data[1]))
        ebars.append(ebar)
        ff = open(fl)
        line = ff.readline().split()
        sigs.append(float(line[0]))
        gams.append(float(line[2]))
        ms.append(float(line[4]))
        mg.append(float(line[6]))
        tbars.append(float(line[8]))
        ff.close()
    ebars = np.array(ebars)
    sigs = np.array(sigs)
    gams = np.array(gams)
    ms = np.array(ms)
    mg = np.array(mg)
    tbars = np.array(tbars)
    mask = np.argsort(ebars)
    ebars = ebars[mask];sigs=sigs[mask];gams=gams[mask];ms=ms[mask];mg=mg[mask];tbars=tbars[mask];
    for i in range(len(ebars)):
        print >>fout,'%d    %1.6f    %1.6f    %1.6f    %1.6f    %1.6f'%(ebars[i],sigs[i],gams[i],ms[i],mg[i],tbars[i])
    fout.close()
    ebars = []
    sigs = []
    gams = []
    ms = []
    mg = []
    tbars = []
    for fl in blist:
        data = fl.split('/')[-1].split('_')
        ebar = np.sqrt(int(data[0])*int(data[1]))
        ebars.append(ebar)
        ff = open(fl)
        line = ff.readline().split()
        sigs.append(float(line[0]))
        gams.append(float(line[2]))
        ms.append(float(line[4]))
        mg.append(float(line[6]))
        tbars.append(float(line[8]))
        ff.close()
    ebars = np.array(ebars)
    sigs = np.array(sigs)
    gams = np.array(gams)
    ms = np.array(ms)
    mg = np.array(mg)
    tbars = np.array(tbars)
    mask = np.argsort(ebars)
    ebars = ebars[mask];sigs=sigs[mask];gams=gams[mask];ms=ms[mask];mg=mg[mask];tbars=tbars[mask];
    for i in range(len(ebars)):
        print >>bout,'%d    %1.6f    %1.6f    %1.6f    %1.6f    %1.6f'%(ebars[i],sigs[i],gams[i],ms[i],mg[i],tbars[i])
    fout.close()