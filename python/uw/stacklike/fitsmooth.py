import skymaps as sk
from uw.utilities.minuit import Minuit
import uw.stacklike.stacklike as b
import numpy as np
import scipy.optimize as so
import glob as glob
import uw.like.pypsf as up
from uw.like.pycaldb import CALDBManager
from uw.like.pypsf import CALDBPsf
import uw.stacklike.angularmodels as ua
#import uw.utilities.assigntasks as uta
from uw.stacklike.stcaldb import scale2
import astropy.io.fits as pf
import pylab as py
import numpy as N
import copy as cp
import matplotlib
import os as os
import time as t
import string
import scipy.interpolate as intp
from uw.stacklike.CLHEP import Photon
import uw.stacklike.binned as ub
import cPickle
from uw.utilities.minuit import Minuit
from IPython.parallel import Client
from scipy.interpolate import UnivariateSpline
import uw.stacklike.dataset as ud


rd = 180/np.pi   #DEG2RAD
inf = 1e40
trad = rd
ft1dir = r'/phys/groups/tev/scratch1/users/Fermi/mar0/data/7.6src/'          #ft1 fits files directory
stacklist = 'agn-psf-study-bright'                    #text file containing 'name ra dec' for each stacked source
cdb = os.environ['CALDB']+'/data/glast/lat'     #CALDB directory
sdir='/phys/users/mar0/sourcelists/'
cachedir = '/phys/groups/tev/scratch1/users/Fermi/mar0/cache/'
np.seterr(all='ignore')

def doublepsfint(dm,pars):
    s1,g1,n1,s2,g2,n2 = pars
    return n1*intg(dm,s1,g1)+n2*intg(dm,s2,g2)

def intg(dm,s1,g1):
    return 1.-(1+(dm/s1)**2/(2*g1))**(-g1+1)

def rfrac(frac,s1,g1):
    return s1*np.sqrt(2*g1*((1.-frac)**(1./(-g1+1.))-1))

def double2single(pars):
    #s1,g1,n1,s2,g2,n2=pars
    norm = doublepsfint(inf,pars)

    r68 = so.fmin_powell(lambda x: (doublepsfint(x[0],pars)/norm-0.68)**2 if x[0]>0 else inf,[1.54*(s1*s2)**0.5],disp=0,full_output=1)
    r95 = so.fmin_powell(lambda x: (doublepsfint(x[0],pars)/norm-0.95)**2 if x[0]>0 else inf,[2.37*(s1*s2)**0.5],disp=0,full_output=1)
    print r68[0],r95[0],r95[0]/r68[0]
    psf = ua.PSF([0,4],[0.1,2.])
    psf.fromcontain([r68[0],r95[0]],[0.68,0.95])
    return psf.model_par

def rfrac2(pars,frac):
    norm = doublepsfint(inf,pars)
    best = so.fmin_powell(lambda x: (doublepsfint(x[0],pars)/norm-frac)**2 if x[0]>0 else inf,[1.54*(pars[0]*pars[2])**0.5],disp=0,full_output=1)
    return best[0].item()

## Main method
# @param ec Conversion Type (0-front 1-back)
def main(ec,slist,name,verbose=0):
    f = Fitter(ec,slist,name,single=False,verbose=verbose)
    return f

######
class ContainIntp(object):

    def __init__(self):
        self.makefunc()

    def makefunc(self):
        gam = np.linspace(1.2,5.,200.)[::-1]
        ratio = np.array([rfrac(0.95,1.,gm)/rfrac(0.68,1.,gm) for gm in gam])
        self.f = UnivariateSpline(ratio,gam,s=0)

    def getpars(self,r68,r95):
        gam = self.f(r95/r68)
        sig = r68/rfrac(0.68,1.,gam)
        return sig,gam

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
    def __init__(self,ec,slist,name,single=True,verbose=0,ctbin='0210',suff=''):
        self.single=single
        self.lastlike=0
        self.suff=suff
        self.ctbin=ctbin
        self.iterations=0
        self.verbose=verbose
        self.ec = ec
        #self.ftdir=ftdir
        self.name=name
        self.slist = slist
        self.intp = ContainIntp()
        #set up the initial condition based on P6v3  
        self.loaddata()
        #self.solve()

    def loaddata(self):
        import cPickle
        #print len(pars),len(steps),len(lims)
        emin=[]
        emax=[]
        irf = 'P6_v11_diff'
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
        
        self.altb=[]

        #files = np.sort(glob.glob(self.ftdir+'*.fit*'))
        template = 'fitter_%s_%s_ec%d_'%(self.name,self.slist,self.ec)
        ## setup the stacker for each energy bin
        gacc = 0
        """for it,en in enumerate(emin):
            temp2 = 'emin_%d_emax_%d.pickle'%(en,emax[it])
            if not os.path.exists(cachedir+template+temp2):
                al = b.StackLoader(name=self.slist,quiet=self.verbose<2)
                al.loadphotons(0.,4.,en,emax[it],0,999999999,self.ec)
                events = len(al.photons)
                #del al.cm.minuit
                al.spickle(cachedir+template+temp2)
            else:
                al = cPickle.load(open(cachedir+template+temp2))
            al.bindata()
            al.solveback()
            #al.solvepsf()
            self.altb.append(al)"""
        template = ud.basedir+'cache/likelihoods*%s%d*%s.pickle'%(self.ctbin,self.ec,self.suff)
        print template
        files = np.sort(glob.glob(template))
        for fil in files:
            print fil
            cl = cPickle.load(open(fil))
            self.altb.append(cl)

    def profilelikelihoods(self,grid=21):
        fdir='/phys/groups/tev/scratch1/users/Fermi/mar0/figures/'
        py.figure(figsize=(8,8))
        for tb in self.altb:
            tb.solvepsf()
            py.clf()
            l0 = tb.cm.minuit.fval
            sran = np.linspace(max(tb.sigma-5*tb.sigmae,0),tb.sigma+5*tb.sigmae,grid)
            gran = np.linspace(max(tb.gamma-5*tb.gammae,0),tb.gamma+5*tb.gammae,grid)
            surf = np.array([tb.solveback(model_par=[sra,gra])-l0 for gra in gran for sra in sran]).reshape(grid,-1)
            py.contourf(sran/tb.sigma,gran/tb.gamma,surf,levels=[0.5,2,4.5,8,12.5,18])
            py.xlim(0,2)
            py.ylim(0,2)
            py.colorbar()
            py.savefig(fdir+'%d%d.png'%(tb.ebar,tb.cls))


    def solve(self,**kwargs):
        self.__dict__.update(kwargs)
        for tb in self.altb:
            tb.quiet = self.verbose<2
        if self.single:
            if self.ec==0:
                pars = [3.88,0.06,0.81,10.62,0.16,0.85]
            else:
                pars = [5.74,0.11,0.76,17.67,0.47,0.87]
            lims = np.array([[0.01,20.],[0.00001,1.],[0.001,2.],[0.01,40.],[0.00001,1.],[0.001,2.]])
        else:
            if self.ec==0:
                pars = [4.0,-0.84,0.25,-0.05,0,11.,-0.84,0.6,0.05,0]
            else:
                pars = [6.7,-0.79,0.6,-0.05,0,21.,-0.79,1.,0.05,0]
            lims = np.array([[0.01,20.],[0.00001,1.],[0.001,2.],[0.00001,1.],[0.001,2.],[0.01,40.],[0.00001,1.],[0.001,2.],[0.00001,1.],[0.001,2.]])

        
        #pars = np.log10(np.array(pars))
        pars = np.array(pars)
        steps = abs(np.array(pars))*0.0001
        lims = np.log10(lims)


        #al.makeplot(name='/phys/users/mar0/figures/%s%d%d'%(slist,al.ebar,al.cls),scale='linear',bin=len(hist[0]))

        fixed = [True for x in pars]
        fixed[0]=False
        ## optimize curves
        #minuit = Minuit(lambda x:self.likelihood(x),pars,fixed=fixed,printMode=2)#limits=lims,steps=steps,
        #minuit.minimize()
        self.best=0
        self.lastbin=np.zeros(len(self.altb))
        sample = abs(self.likelihood(pars))
        self.best=sample
        ftol = 0.5/sample
        if self.verbose>0:
            print 'Using function tolerance: %1.2g (%1.1f)'%(ftol,sample)
        almin = so.fmin_powell(lambda x:self.likelihood(x),pars,ftol=ftol,disp=1,full_output=1,callback=self.update(self.lastlike))
        
        ## print out fit parameters
        #print 10**almin[0]
        #self.pars = 10**almin[0]
        print almin[0]
        self.pars = almin[0]
        #print 10**minuit.params
        #self.pars = 10**minuit.params
        self.fval = almin[1]
        return self.pars

        
    def update(self,last):
        self.best=last

    ## likelihood function
    #  make a king function from the 68 and 95% containment calculated from the smooth
    #  function for each energy bin
    #
    #  fit a combination of point-source and isotropic background
    #  and return likelihood
    
    def likelihood(self,pars):#c068,c168,b68,c268,b168,c095,c195,b95,c295,b195):
        #pars = 10**np.array(pars)
        #print pars

        #catch garbage values
        if np.any(np.isnan(pars)):# or np.isnan(c095):
            return 1e40
        acc=0.

        if self.verbose>1:
            print '######################################'

            print 'iteration: %d'%self.iterations
        ## calculate the likelihood of the parameterization for each energy bin
        for it,tb in enumerate(self.altb):
            eb = tb.ebar
            if self.single:
                c068,c168,b68,c095,c195,b95 = pars
                r68 = scale2(eb,c068,c168,-b68)/rd
                r95 = scale2(eb,c095,c195,-b95)/rd
            else:
                c068,b68,c268,b168,a68,c095,b95,c295,b195,a95 = pars
                r68 = self.scale(eb,pars[:len(pars)/2])/rd
                r95 = self.scale(eb,pars[len(pars)/2:])/rd

            #R68 should always be less than R95!
            if self.single:
                check = (r95/r68<1.84 or r95/r68>100 or c068<0 or c095<0)
            else:
                check = (r95/r68<1.84 or r95/r68>100 or c068<0 or c095<0 or c268<0 or c295<0)

            if check:
                continue
            else:
                if not self.single:
                    if c068<(10.*c268) or c095<(10.*c295):
                        continue
                    else:
                        sig,gam = self.intp.getpars(r68*rd,r95*rd)
                        like = tb.psflikelihood([sig,gam])#like = tb.rcontainlike(ct=[r68,r95])
                else:
                    sig,gam = self.intp.getpars(r68*rd,r95*rd)#like = tb.rcontainlike(ct=[r68,r95])
                    like = tb.psflikelihood([sig,gam])
            if self.iterations==0:
                self.lastbin[it]=like
            acc = acc+like
            if self.verbose>1:
                print 'Energy: %1.0f  sigma= %1.3f  gamma= %1.3f  R68: %1.2f  R95: %1.2f  like: %1.0f  diff= %1.0f'%(eb,sig,gam,r68*rd,r95*rd,like,self.lastbin[it]-like)
            self.lastbin[it]=like
        if self.verbose>1:
            print '--------------------------------------'
            if self.single:
                print '68  c0: %1.2f  c1: %1.2f  b: %1.2f'%(c068,c168,b68)
                print '95  c0: %1.2f  c1: %1.2f  b: %1.2f'%(c095,c195,b95)
            else:
                print '68  c0: %1.2f  b: %1.2f  c2: %1.2f  b1: %1.2f'%(c068,b68,c268,b168)
                print '95  c0: %1.2f  b: %1.2f  c2: %1.2f  b1: %1.2f'%(c095,b95,c295,b195)
            print ' Total: %1.2f      Diff: %1.2f'%(acc,acc-self.lastlike)
            print '######################################'
        pstring = ''
        for par in pars:
            pstring += '%1.3e   '%(par)
        pstring += '%1.2f'%(acc-self.best)
        if self.verbose>0:
            print pstring
        self.lastlike=acc
        self.iterations=self.iterations+1
        return acc

    def qprofile(self,it):
        ftol = 0.01/0.25
        almin = so.fmin_powell(lambda x: (self.fval-self.likelihood([x[0] if it==it2 else self.pars[it2] for it2 in range(len(self.pars))])+0.5)**2,[self.pars[it]],disp=1,full_output=1,ftol=ftol)
        return almin

    def solveerrs(self):
        self.r68e,self.r95e = [],[]
        for tb in self.altb:
            tb.solvepsf()
            self.r68e.append(tb.r68e)
            self.r95e.append(tb.r95e)

    def getuncs(self):
        self.uncs=[]
        for it,par in enumerate(self.pars):
            sig1 = self.qprofile(it)
            self.uncs.append(np.sqrt((sig1[0].item()-par)**2))

    ## function that describes the smoothed R68 and R95
    def scale(self,e,pars):
        c0,b,c2,b1,a0=pars
        return np.sqrt((c0*((e/100.)**(b+a0*np.log(e))))**2+(c2*((e/100.)**b1))**2)

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
   

def scale2unc(en,pars,uncs):
    t0,t1,b0 = pars
    st0,st1,sb0 = uncs

    fv = scale2(en,t0,t1,-b0)
    ex = fv**2-t1**2
    err1 = (np.sqrt(ex)/fv*st0)**2
    err2 = (t1/fv*st1)**2
    err3 = (ex*np.log(en/100.)/fv*sb0)**2
    func = err1+err2+err3
    #print np.sqrt(err1),np.sqrt(err2),np.sqrt(err3),np.sqrt(ex),fv
    return np.sqrt(func)

## calculates the chisq for a both front and parameter scale function and the sigma parameter
#  helps set up the scaling function paramters
#  @param el energy list
#  @param er1 sigma error list
#  @param c0 PSF sigma @ 100 MeV
#  @param c1 instrument pitch
#  @param b power law of multiple scattering
#  @param 
def chisq2(sigmaf,sigmab,el,erf,erb,c0f,c1f,c0b,c1b,b):
    if b>0 or c0f<0 or c1f<0 or c0b<0 or c1b<0:
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
def scale3(e,pars):
    c0,b,c2,b1 = pars
    return N.sqrt((c0*((e/100.)**(b)))**2+(c2*((e/100.)**b1))**2)



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

    ######## P7_v8 5 pars cgv ###############
    r68fp = [ 4.20809294,0., -0.80217512,   0.26678933,   -0.1567494]
    r95fp = [15.33020716,0., -0.81261126,   0.14620895,  0.24271978]
    r68bp = [7.95071816, 0.,-0.84758715,   0.43720761,   -0.11798399]
    r95bp = [29.6898745, 0.,-0.84844305,   0.2916015,   0.26842457]

    #0.8-1.0
    r68fp = [ 4.04011621,0,   -0.82593281,   0.30040763,   -0.17397118]
    r95fp = [14.08348845,0,   -0.82971535,   0.13511403,  0.2576134]
    r68bp = [7.40104724,0,   -0.87068742,   0.40068761,   -0.11358863]
    r95bp = [25.74968446,0,   -0.85169841,   0.16804804,  0.29550437]
    
    #0.6-0.8
    #r68fp=[  4.56939118,0,   -0.81566421,   0.24209374,   -0.14208137]
    #r95fp= [    17.28561549,0,   -0.84011464,   0.10370819,  0.31374605]
    #r68bp=    [  8.25820957,0,   -0.84130906,   0.40417903,   -0.10810866]
    #r95bp=    [28.76750319,0,   -0.81233442,   0.47846827,  0.17074784]

    #0.2-0.6
    #r68fp=[  5.45433437,0,   -0.791509374,   0.191437601, -0.0935396294]
    #r95fp=[25.6334213,0,   -0.884321342, 0.766460729,   -0.0251667491]
    #r68bp=[  8.78455277,0,   -0.80194256,   0.43848317,   -0.10880232]
    #r95bp= [  28.36437816, 0,  -0.84534901,   0.94483145,  0.11911295]
    
    """#SOURCE02-06
    allpars=[[[5.23857857e+00,   -8.20795258e-01,   1.75548040e-01, -7.95343457e-02],
    [   2.08210327e+01, -8.52300971e-01,   7.75651267e-01,   -9.10045251e-03],
    [8.24621427e+00,  -8.27639242e-01,   3.64923474e-01,    -8.65548078e-02],
    [2.28680747e+01,-7.07672505e-01,   1.32342024e+00,   7.46508133e-02]],

    #06-08
    [[  4.60192931e+00,  -8.57376328e-01,   1.68603947e-01,-9.22038520e-02],
    [1.76176859e+01,-8.59780731e-01,   3.16241369e-01,   6.98472663e-02],
    [  7.82112178e+00,  -8.84973756e-01,   2.66668147e-01, -8.46488723e-02],
    [   2.33117324e+01,-6.94111810e-01,   4.19253591e-01,   8.03534836e-02]],

    #08-10
    [[  4.11266993e+00,  -8.67261035e-01,   2.50921530e-01,-1.45804873e-01],
    [   1.45770725e+01, -8.24783573e-01,   1.17220829e-01,   2.90414166e-01],
    [  7.51906842e+00,  -9.41678371e-01,   2.75720534e-01,-3.24545580e-02], 
    [   2.49499271e+01, -7.55555753e-01,   7.99968308e-02,   4.96254239e-01]]]"""

    #CLEAN02-06
    allpars = [[[4.97427873e+00,  -8.16872319e-01,   2.03156702e-01,-7.62511798e-02],
    [   1.90086552e+01,   -8.89804859e-01,   5.92363712e-01,   8.57664341e-02],
    [  7.48915166e+00,  -8.16907271e-01,   3.68621132e-01, -7.36714911e-02],
    [2.31540210e+01, -7.56993104e-01,   1.46246091e+00,   7.03964889e-02]],

    #06-08
    [[  4.34627407e+00,  -8.61956353e-01,   1.72604875e-01, -1.05576383e-01],
    [1.50449030e+01, -8.22714121e-01,   2.58606419e-01,   1.01315119e-01],
    [  7.71643952e+00,  -8.83027878e-01,   2.27643235e-01,-5.86944443e-02],
    [2.33047514e+01, -7.16012650e-01,   3.87204673e-01,   8.34804948e-03]],

    #08-10
    [[  3.97339407e+00,  -8.64029541e-01,   2.55489154e-01,-1.50298532e-01],
    [1.35283285e+01,-8.22397000e-01,   1.25321889e-01,   2.76959006e-01],
    [  7.27692814e+00,  -9.15601608e-01,   2.71996594e-01,-2.54543408e-02],
    [2.36346689e+01,-8.05587385e-01,   8.41914171e-02,   4.71889582e-01]]]
    

    #02-10
    """allpars = [[[  4.21686675e+00,  -8.44843190e-01,   2.34435814e-01,-1.56357794e-01],
    [   1.53165076e+01, -8.15209164e-01,   1.46550021e-01,   2.09081999e-01],
    [  7.37424681e+00,  -8.79920587e-01,   3.09835237e-01,-9.58531901e-02],
    [   2.28097448e+01,-7.17841791e-01,   3.83354674e-01,   1.25992575e-01]],
    [[  4.21686675e+00,  -8.44843190e-01,   2.34435814e-01,-1.56357794e-01],
    [   1.53165076e+01, -8.15209164e-01,   1.46550021e-01,   2.09081999e-01],
    [  7.37424681e+00,  -8.79920587e-01,   3.09835237e-01,-9.58531901e-02],
    [   2.28097448e+01,-7.17841791e-01,   3.83354674e-01,   1.25992575e-01]],
    [[  4.21686675e+00,  -8.44843190e-01,   2.34435814e-01,-1.56357794e-01],
    [   1.53165076e+01, -8.15209164e-01,   1.46550021e-01,   2.09081999e-01],
    [  7.37424681e+00,  -8.79920587e-01,   3.09835237e-01,-9.58531901e-02],
    [   2.28097448e+01,-7.17841791e-01,   3.83354674e-01,   1.25992575e-01]]]"""

    xvec = [0.52,0.72,0.92]
    xi2 = [xi*xi for xi in xvec]
    xi = sum(xvec)
    ii = len(xi2)
    detA = sum(xi2)*ii-(xi)**2

    pname = cdb+r'/bcf/psf/psf_%s_'%irf
    print pname
    cwd = os.getcwd()
    cmd1 = r'cp %sfront.fits %s/psf_%s_front.fits'%(pname,cwd,out)
    cmd2 = r'cp %sback.fits %s/psf_%s_back.fits'%(pname,cwd,out)
    print cmd1
    print cmd2
    os.system(cmd1)
    os.system(cmd2)
    ff = pf.open('%s/psf_%s_front.fits'%(cwd,out),mode='update')
    bf = pf.open('%s/psf_%s_back.fits'%(cwd,out),mode='update')

    cth = ff[1].data.field('CTHETA_LO')[0]
    ctbar = (ff[1].data.field('CTHETA_LO')[0]+ff[1].data.field('CTHETA_HI')[0])*0.5
    en = ff[1].data.field('ENERG_LO')[0]
    ct = len(ff[1].data.field('CTHETA_LO')[0])
    e = len(ff[1].data.field('ENERG_LO')[0])

    ens = []
    fss =[]
    bss =[]

    intp = ContainIntp()

    for i in N.arange(0,e,1):
        energy = N.sqrt(ff[1].data.field('ENERG_LO')[0][i]*ff[1].data.field('ENERG_HI')[0][i])
        ens.append(energy)
        """r68f = scale2(energy,r68fp[0],r68fp[1],r68fp[2])
        r95f = scale2(energy,r95fp[0],r95fp[1],r95fp[2])
        r68b = scale2(energy,r68bp[0],r68bp[1],r68bp[2])
        r95b = scale2(energy,r95bp[0],r95bp[1],r95bp[2])"""
        yi68f= []
        yi95f=[]
        yi68b=[]
        yi95b=[]
        for tpars in allpars:
            yi68f.append(scale3(energy,tpars[0]))
            yi95f.append(scale3(energy,tpars[1]))
            yi68b.append(scale3(energy,tpars[2]))
            yi95b.append(scale3(energy,tpars[3]))
        #t.sleep(1.)
        print energy,yi68f,yi95f
        dotp = lambda x: sum([xvec[it]*x[it] for it in range(len(xvec))])
        mcalc = lambda x: (dotp(x)*ii-xi*sum(x))/detA
        bcalc = lambda x: (sum(xi2)*sum(x)-xi*dotp(x))/detA
        def lin(m,b,x):
            return m*x+b
        m68f = mcalc(yi68f)
        b68f = bcalc(yi68f)
        m95f = mcalc(yi95f)
        b95f = bcalc(yi95f)
        m68b = mcalc(yi68b)
        b68b = bcalc(yi68b)
        m95b = mcalc(yi95b)
        b95b = bcalc(yi95b)
        #print m68f,b68f
        #fpsf.fromcontain(1.,r95f/r68f,0.68,0.95)
        #bpsf.fromcontain(1.,r95b/r68b,0.68,0.95)

        r68f = lin(m68f,b68f,0.82)
        r95f = lin(m95f,b95f,0.82)
        r68b = lin(m68b,b68b,0.82)
        r95b = lin(m95b,b95b,0.82)
        #print r68f,r95f
        #print r68b,r95b

        fsigma,fgamma = intp.getpars(r68f,r95f)#fpsf.model_par[0]*r68f
        bsigma,bgamma = intp.getpars(r68b,r95b)#bpsf.model_par[0]*r68b
        fsigma/=rd
        bsigma/=rd
        #fgamma = fpsf.model_par[1]
        #bgamma = bpsf.model_par[1]
        fss.append(fsigma)
        bss.append(bsigma)

        scf = scale(energy,0)
        scb = scale(energy,1)

        t68 = []
        t95 = []
        for j in range(ct):
        
            
            r68f = lin(m68f,b68f,ctbar[j])
            r95f = lin(m95f,b95f,ctbar[j])
            r68b = lin(m68b,b68b,ctbar[j])
            r95b = lin(m95b,b95b,ctbar[j])
            t68.append(r68f);t95.append(r95f)
            #print r68f,r95f
            #print r68b,r95b
            #t.sleep(1.)
            fsigma,fgamma = intp.getpars(r68f,r95f)#fpsf.model_par[0]*r68f
            bsigma,bgamma = intp.getpars(r68b,r95b)#bpsf.model_par[0]*r68b
            fsigma/=rd
            bsigma/=rd
            try:
                """print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,ff[1].data.field('SIGMA')[0][j][i],fsigma*trad/scf)
                print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,bf[1].data.field('SIGMA')[0][j][i],bsigma*trad/scb)"""
                ff[1].data.field('SIGMA')[0][j][i]=fsigma*trad/scf
                bf[1].data.field('SIGMA')[0][j][i]=bsigma*trad/scb
            except:
                pass#print 'New style psf'
            try:
                """print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,ff[1].data.field('SCORE')[0][j][i],fsigma*trad/scf)
                print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,bf[1].data.field('SCORE')[0][j][i],bsigma*trad/scb)
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,ff[1].data.field('STAIL')[0][j][i],fsigma*trad/scf)
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,bf[1].data.field('STAIL')[0][j][i],bsigma*trad/scb)"""
                ff[1].data.field('SCORE')[0][j][i]=fsigma*trad/scf
                bf[1].data.field('SCORE')[0][j][i]=bsigma*trad/scb
                ff[1].data.field('STAIL')[0][j][i]=fsigma*trad/scf
                bf[1].data.field('STAIL')[0][j][i]=bsigma*trad/scb
            except:
                pass#print 'Old style psf'
            """print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,ff[1].data.field('NTAIL')[0][j][i],0)
            print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,bf[1].data.field('NTAIL')[0][j][i],0)
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,ff[1].data.field('GCORE')[0][j][i],fgamma)
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,bf[1].data.field('GCORE')[0][j][i],bgamma)
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,ff[1].data.field('GTAIL')[0][j][i],fgamma)
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,bf[1].data.field('GTAIL')[0][j][i],bgamma)"""
            ff[1].data.field('NTAIL')[0][j][i]=0
            bf[1].data.field('NTAIL')[0][j][i]=0
            ff[1].data.field('GCORE')[0][j][i]=fgamma
            bf[1].data.field('GCORE')[0][j][i]=bgamma
            ff[1].data.field('GTAIL')[0][j][i]=fgamma
            bf[1].data.field('GTAIL')[0][j][i]=bgamma
        #t.sleep(2)
        print t68
        print t95
    fss = N.array(fss)
    fse = fss**1.5
    bss = N.array(bss)
    bse = bss**1.5


    print fss
    print bss
    pars = [58e-3,377e-6,96e-3,1300e-6,-0.8]
    #best = Minuit(lambda x: chisq2(fss,bss,ens,fse,bse,x[0],x[1],x[2],x[3],x[4]),params=pars,up=1,printMode=0,strategy=2,tolerance=1e-12)
    #best.minimize()
    best = so.fmin_powell(lambda x: chisq2(fss,bss,ens,fse,bse,x[0],x[1],x[2],x[3],x[4]),pars,ftol=1e-12,full_output=1,disp=0)
    bestparams = best[0]
    print bestparams

    for it,energy in enumerate(ens):
        print energy,fss[it]*rd,scale2(energy,bestparams[0],bestparams[1],bestparams[4])*rd,bss[it]*rd,scale2(energy,bestparams[2],bestparams[3],bestparams[4])*rd

    #t.sleep(60.)
    for it,par in enumerate(bestparams):
        ff[2].data.field('PSFSCALE')[0][it]=par
        bf[2].data.field('PSFSCALE')[0][it]=par

    for i in N.arange(0,e,1):
        energy = N.sqrt(ff[1].data.field('ENERG_LO')[0][i]*ff[1].data.field('ENERG_HI')[0][i])
        """r68f = scale2(energy,r68fp[0],r68fp[1],r68fp[2])
        r95f = scale2(energy,r95fp[0],r95fp[1],r95fp[2])
        r68b = scale2(energy,r68bp[0],r68bp[1],r68bp[2])
        r95b = scale2(energy,r95bp[0],r95bp[1],r95bp[2])"""

        yi68f= []
        yi95f=[]
        yi68b=[]
        yi95b=[]
        for tpars in allpars:
            yi68f.append(scale3(energy,tpars[0]))
            yi95f.append(scale3(energy,tpars[1]))
            yi68b.append(scale3(energy,tpars[2]))
            yi95b.append(scale3(energy,tpars[3]))
        dotp = lambda x: sum([xvec[it]*x[it] for it in range(len(xvec))])
        mcalc = lambda x: (dotp(x)*ii-xi*sum(x))/detA
        bcalc = lambda x: (sum(xi2)*sum(x)-xi*dotp(x))/detA
        def lin(m,b,x):
            return m*x+b
        m68f = mcalc(yi68f)
        b68f = bcalc(yi68f)
        m95f = mcalc(yi95f)
        b95f = bcalc(yi95f)
        m68b = mcalc(yi68b)
        b68b = bcalc(yi68b)
        m95b = mcalc(yi95b)
        b95b = bcalc(yi95b)
        print m68f,b68f
        #fpsf.fromcontain(1.,r95f/r68f,0.68,0.95)
        #bpsf.fromcontain(1.,r95b/r68b,0.68,0.95)

        r68f = lin(m68f,b68f,0.82)
        r95f = lin(m95f,b95f,0.82)
        r68b = lin(m68b,b68b,0.82)
        r95b = lin(m95b,b95b,0.82)
        print r68f,r95f
        print r68b,r95b

        fsigma,fgamma = intp.getpars(r68f,r95f)#fpsf.model_par[0]*r68f
        bsigma,bgamma = intp.getpars(r68b,r95b)#bpsf.model_par[0]*r68b
        fsigma/=rd
        bsigma/=rd
        #fgamma = fpsf.model_par[1]
        #bgamma = bpsf.model_par[1]
        #fss.append(fsigma)
        #bss.append(bsigma)

        scf = scale2(energy,bestparams[0],bestparams[1],bestparams[4])
        scb = scale2(energy,bestparams[2],bestparams[3],bestparams[4])

        for j in range(ct):

            r68f = lin(m68f,b68f,ctbar[j])
            r95f = lin(m95f,b95f,ctbar[j])
            r68b = lin(m68b,b68b,ctbar[j])
            r95b = lin(m95b,b95b,ctbar[j])
            print r68f,r95f
            print r68b,r95b
            #t.sleep(1.)
            fsigma,fgamma = intp.getpars(r68f,r95f)#fpsf.model_par[0]*r68f
            bsigma,bgamma = intp.getpars(r68b,r95b)#bpsf.model_par[0]*r68b
            fsigma/=rd
            bsigma/=rd
            try:
                print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,ff[1].data.field('SIGMA')[0][j][i],fsigma*trad/scf)
                ff[1].data.field('SIGMA')[0][j][i]=fsigma/scf
                print 'E = %1.0f Changing sigma from %1.3f -> %1.3f'%(energy,bf[1].data.field('SIGMA')[0][j][i],bsigma*trad/scb)
                bf[1].data.field('SIGMA')[0][j][i]=bsigma/scb
            except:
                print 'New style psf'
            try:
                print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,ff[1].data.field('SCORE')[0][j][i],fsigma*trad/scf)
                ff[1].data.field('SCORE')[0][j][i]=fsigma/scf
                print 'E = %1.0f Changing score from %1.3f -> %1.3f'%(energy,bf[1].data.field('SCORE')[0][j][i],bsigma*trad/scb)
                bf[1].data.field('SCORE')[0][j][i]=bsigma/scb
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,ff[1].data.field('STAIL')[0][j][i],fsigma*trad/scf)
                ff[1].data.field('STAIL')[0][j][i]=fsigma/scf
                print 'E = %1.0f Changing stail from %1.3f -> %1.3f'%(energy,bf[1].data.field('STAIL')[0][j][i],bsigma*trad/scb)
                bf[1].data.field('STAIL')[0][j][i]=bsigma/scb
            except:
                print 'Old style psf'
            print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,ff[1].data.field('NTAIL')[0][j][i],0)
            ff[1].data.field('NTAIL')[0][j][i]=0
            print 'E = %1.0f Changing NTAIL from %1.2f -> %1.2f'%(energy,bf[1].data.field('NTAIL')[0][j][i],0)
            bf[1].data.field('NTAIL')[0][j][i]=0
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,ff[1].data.field('GCORE')[0][j][i],fgamma)
            ff[1].data.field('GCORE')[0][j][i]=fgamma
            print 'E = %1.0f Changing gcore from %1.2f -> %1.2f'%(energy,bf[1].data.field('GCORE')[0][j][i],bgamma)
            bf[1].data.field('GCORE')[0][j][i]=bgamma
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,ff[1].data.field('GTAIL')[0][j][i],fgamma)
            ff[1].data.field('GTAIL')[0][j][i]=fgamma
            print 'E = %1.0f Changing gtail from %1.2f -> %1.2f'%(energy,bf[1].data.field('GTAIL')[0][j][i],bgamma)
            bf[1].data.field('GTAIL')[0][j][i]=bgamma

    ff.flush()
    bf.flush()
    cmd1 = r'cp %s/psf_%s_front.fits %s'%(cwd,out,cdb+r'/bcf/psf/')
    cmd2 = r'cp %s/psf_%s_back.fits %s'%(cwd,out,cdb+r'/bcf/psf/')
    print os.system(cmd1)
    print os.system(cmd2)

## makeplots makes plots of the PSF for different IRFS
#  @param irfs list of IRF definitions
#  @param fact smoothing factor (averages the value at x*fact and x/fact)
#  @param num bins per decade
#  @param cts containment fraction plotted 
def makeplots(irfs=['P6_v11_diff'],fact=1.,num=32.,cts=[34,68,95]):
    
    nms = [str(x)+'% '+y for x in cts for y in irfs]
    clr = ['b','k','g','r']
    sty2 = ['-','--','-.',':']
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
            #p1 = py.plot(energy,[np.sqrt(psf.inverse_integral(x/fact,0,ct))*np.sqrt(psf.inverse_integral(x*fact,0,ct)) for x in energy],'%s%s'%(clr[it2],sty2[it1]))
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
            #p1 = py.plot(energy,[np.sqrt(psf.inverse_integral(x/fact,1,ct))*np.sqrt(psf.inverse_integral(x*fact,1,ct)) for x in energy],'%s%s'%(clr[it2],sty2[it1]))
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

################################################    THETA DEPENDENT METHOD    ####################################################

## likelihoodtable - calculates the best parameters for a given CALDB table element based on neighbors
#  @param enr1 energy index in the table [0-17]
#  @param ctr1 cosine theta index in table [0-7]
#  @param ctype conversion type 0-front, 1-back
#  @param weight defines the smoothing kernel size, \propto sigma^-1
#  @param irf reference CALDB file, used to extract scaling function
#  @param out output CALDB file
#  @param jg generate likelihood pickles, but don't analyze

def likelihoodtable(enr1,ctr1,ctype,weight=0.5,irf ='P7SOURCE_V6',mcirf='P7SOURCE_V4MC',out='P7SOURCE_V11',jg=False):
    import cPickle
    days =4*365      #number of days of data to use
    bins = 12      #number of angular bins
    pname = os.environ['CALDB']+r'/data/glast/lat/bcf/psf/psf_%s_'%irf
    pulsdir = r'/phys/groups/tev/scratch1/users/Fermi/mar0/data/7.3srcrp/pulsar/'
    agndir = r'/phys/groups/tev/scratch1/users/Fermi/mar0/data/7.3srcrp/'
    suffix = 'clean'
    print pname

    ######  get reference CALDB details  ###########
    ff = pf.open('%sfront.fits'%pname)
    bf = pf.open('%sback.fits'%pname)
    cthlo = ff[1].data.field('CTHETA_LO')[0]         #cosine theta bin edges
    cthhi = ff[1].data.field('CTHETA_HI')[0]
    enlo = ff[1].data.field('ENERG_LO')[0]           #energy bin edges
    enhi = ff[1].data.field('ENERG_HI')[0]
    ctbins = len(ff[1].data.field('CTHETA_LO')[0])   #number of cosine theta and energy bins
    ebins = len(ff[1].data.field('ENERG_LO')[0])

    #######  define the binning for the analysis, doesn't have to match the reference CALDB file       #######################
    a_ebins = np.array([100,133,178,237,316,421,562,750,1000,1778,3162,5623,10000,17783,31623,100000])   #energy bin edges
    a_ebars = np.sqrt(a_ebins[:-1]*a_ebins[1:])                                                          #center of energy bin
    a_ctbins = np.array([0.2,0.6,0.8,1.0])                                                               #cosine theta bin edges
    a_ctbar = (a_ctbins[1:]+a_ctbins[:-1])/2.                                                            #center of cosine theta bins
    #a_eidx = np.array([ for aeb in a_ebars])                                #convert energies to table indices
    #a_cidx = 10*a_ctbar-2.5                                                                              #convert cosin thetas to table indices

    def mkeidx(aeb):
        return 4*np.log10(aeb/min(enlo))-0.5

    def mkcidx(ct):
        return 10*ct-2.5

    ########### Extract scaling function from reference file  #################
    psfsc = ff[2].data[0][0]
    psf = CALDBPsf(CALDBManager(irf=irf))

    ########### Setup Likelihoods for each energy and cosine theta bin for analysis #############
    likelihoods=[]
    for emin,emax in zip(a_ebins[:-1],a_ebins[1:]):
        enlike = []
        ebar = np.sqrt(emin*emax)
        if True:
            tctbins = a_ctbins
        else:
            tctbins = np.array([a_ctbins[0],a_ctbins[-1]])
        for ctmin,ctmax in zip(tctbins[:-1],tctbins[1:]):

            maxr = psf.inverse_integral(ebar,ctype,95)*2.                    #Go out to tails of angular distribution for the ROI cut
            agnlis = ['psfstudy'+suffix] if ebar>1778 else []              #Use AGN information above 1.7 GeV
            psrs = ub.pulsars if ebar<10000 else []             #Use Pulsar information below 10 GeV
            def formatint(pint,num):
                return (num-len(str(pint)))*'0'+str(pint)

            #set up pickle
            pickles = '/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/likelihoods%s%s%s%s%d%s%s.pickle'%(formatint(emin,6),formatint(emax,6),formatint(int(ctmin*10),2),formatint(int(ctmax*10),2),ctype,irf,suffix)

            #load pickle if it exists
            if os.path.exists(pickles):
                pfile = open(pickles)
                print 'loaded %s'%pickles
                cl = cPickle.load(pfile)

            #make a new one
            else:
                bins = 16 if emin<5623 else 12
                cl = ub.CombinedLike(irf=irf,mode=-1,pulsars=psrs[0:2],agnlist=agnlis,verbose=False,ctmin=ctmin,ctmax=ctmax,agndir=agndir,pulsdir=pulsdir)
                cl.loadphotons(0,maxr,emin,emax,239557417,239517417+days*86400,ctype)
                cl.bindata(bins)
                cl.fit(qandd=True)

                #Minuit objects cannot be pickled!
                del cl.minuit
                pfile = open(pickles,'w')
                cPickle.dump(cl,pfile)
                pfile.close()
                cl.makeplot('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/likelihoods%s%s%s%s%d%s%s.png'%(formatint(emin,6),formatint(emax,6),formatint(int(ctmin*10),2),formatint(int(ctmax*10),2),ctype,cl.irf,suffix))
            del psrs
            enlike.append(cp.copy(cl))
        likelihoods.append(enlike)

    ## lin - 2D linear helper function
    #  @param pars [b,m1,m2], intercept and slope1/slope2
    #  @param x1 dimension 1 value
    #  @param x2 dimension 2 value
    def lin(pars,x1,x2):
        return pars[0]+pars[1]*x1+pars[2]*x2

    ## weightedlike - calculates the likelihood of a particular linear parameterization (helper)
    #  @param likel 2D array of likelihoods [energy bin][costheta bin]
    #  @param ti reference table row
    #  @param tj reference table column
    #  @param pars parameters for the king function, passed from fmin
    #  @param scl inverse kernel size
    #  @param spars sigma scaling function parameters, extracted from CALDB
    #  @param ct conversion type
    def weightedlike(likel,ti,tj,pars,scl,spars,ct,ref=0):
        acc =0
        #energy iteration
        for i in range(len(likel)):
            #cosine theta iteration
            for j in range(len(likel[i])):
                eidx = mkeidx(likel[i][j].ebar)
                cidx = mkcidx(likel[i][j].thetabar)
                scb = scale2(a_ebars[i],spars[ct*2],spars[ct*2+1],spars[4])*rd  #calulate scale width at the mean bin energy in degrees
                dist = (((ti-eidx)*scl)**2+((tj-cidx)*scl)**2)        #distance between current table position and reference table position
                if True:
                    fact = np.exp(-dist/2.)/np.sqrt(2*np.pi)/scl                                         #weighting factor for likelihood is just the gaussian
                    tss=lin(pars[0:3],eidx,cidx)                          #SCORE linear approx
                    tsg=lin(pars[3:6],eidx,cidx)                          #GCORE linear approx
                    tss2 = lin(pars[6:9],eidx,cidx)#tss                                                      #set STAIL = SCORE for the time
                    tsg2 = lin(pars[9:12],eidx,cidx)                       #GTAIL linear approx
                    tsf = lin(pars[12:15],eidx,cidx)                       #NCORE/NTAIL linear approx, will be determined by min(GCORE,GTAIL)

                    #make sure parameters are bounded, otherwise they don't contribute to the likelihood
                    if tss<0 or tsg<1 or tsg2<1 or tss2<0 or tsf>1 or tsf<0 or tsg>5. or tsg2>5.:
                        continue
                    
                    #evaluate the likelihood of the parameterization in this bin, weight it, and sum
                    tlike = likel[i][j].psflikelihood([scb*tss,tsg,tsf,scb*tss2,tsg2])#[scb*tss,tsg])
                    acc += tlike*fact
                    #print '%d %d %1.0f %1.1f %1.3f %1.3f %1.1f %1.2f %1.1f %1.1f'%(eidx,cidx,a_ebars[i],dist,tss,tsg,fact,tlike,tlike*fact,acc)
                    #t.sleep(0.1)
        print string.join(['%1.4f'%(prs) for prs in pars],' ') + ' %1.3f'%(np.sqrt(abs(ref-acc))*np.sign(ref-acc))
        #print '#################################################'
        return acc

    #if not the generation step, run analysis
    if not jg:
        idx = ebins*ctr1+enr1             #determine unique index for entry in table

        initial = [0.9,0,0,2.1,0,0,1.1,0,0,1.9,0,0,0.5,0,0]
        l0 = weightedlike(likelihoods,enr1,ctr1,initial,weight,psfsc,ctype)
        ftol = abs(0.05/l0)
        #calcuate best parameters
        print ftol
        if True:
            minuit = so.fmin_powell(lambda x: weightedlike(likelihoods,enr1,ctr1,x,weight,psfsc,ctype,l0),initial,full_output=1,disp=1,ftol=ftol)
            prs = minuit[0]                   
        else:
            minuit = Minuit(lambda x: weightedlike(likelihoods,enr1,ctr1,x,weight,psfsc,ctype,l0),initial,mode=1)
            minuit.minimize()
            prs = minuit.params
        bests = lin(prs[0:3],enr1,ctr1)
        bestg = lin(prs[3:6],enr1,ctr1)
        bests2 = lin(prs[6:9],enr1,ctr1)
        bestg2 = lin(prs[9:12],enr1,ctr1)
        bestf = lin(prs[12:15],enr1,ctr1)
        outfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/%d_%d.txt'%(idx,ctype),'w')

        #output the best fit parameters for this table entry
        print >>outfile,bests,bestg,bests2,bestg2,bestf
        print bests,bestg,bests2,bestg2,bestf
        return idx,bests,bestg,bests2,bestg2,bestf

## makelikeirfs - read in best parameters and create the corresponding CALDB files
#  @param irf reference CALDB file
#  @param out output CALDB file
def makelikeirfs(irf='P7SOURCE_V6',out='P7SOURCE_V11'):
    
    #get the parameter files
    os.chdir('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/')
    frontfiles = ['%d_0.txt'%x for x in range(144)]
    backfiles = ['%d_1.txt'%x for x in range(144)]
    pname = os.environ['CALDB']+r'/data/glast/lat/bcf/psf/psf_%s_'%irf
    cwd = os.getcwd()

    #copy in reference file to modify
    cmd1 = r'cp %sfront.fits %s/psf_%s_front.fits'%(pname,cwd,out)
    cmd2 = r'cp %sback.fits %s/psf_%s_back.fits'%(pname,cwd,out)
    os.system(cmd1)
    os.system(cmd2)
    ff = pf.open('psf_%s_front.fits'%out,mode='update')
    bf = pf.open('psf_%s_back.fits'%out,mode='update')
    ftb=ff[1].data
    btb=bf[1].data
    dim1 = len(ftb.field('SCORE')[0][0])
    dim2 = len(ftb.field('SCORE')[0])

    #update the front CALDB file entries
    for ffil in frontfiles:
        tf = open(ffil)
        idx = int(ffil.split('_')[0])

        enidx=idx%dim1
        ctidx=idx/dim1
        line = tf.readline()
        """sigma,gamma=line.split(' ')
        sigma,gamma=float(sigma),float(gamma)
        #print idx,sigma,gamma
        ftb.field('SCORE')[0][idx]=sigma
        ftb.field('STAIL')[0][idx]=sigma
        ftb.field('GCORE')[0][idx]=gamma
        ftb.field('GTAIL')[0][idx]=gamma
        ftb.field('NTAIL')[0][idx]=0"""
        lines = line.split(' ')
        lines = [float(line) for line in lines]
        #print idx,sigma,gamma


        #figure out if we're looking at core or tail
        if lines[1]>lines[3]:
            score,gcore,stail,gtail,ncore=lines
            ntail= 1-ncore
        else:
            stail,gtail,score,gcore,ntail=lines
            ncore= 1-ntail
        ncore = max(min(1,ncore),1e-4)
        ntail = max((1-ncore)/ncore*(score/stail)**2,0)
        print enidx,ctidx,score,gcore,stail,gtail,ncore
        ftb.field('SCORE')[0][ctidx][enidx]=score
        ftb.field('STAIL')[0][ctidx][enidx]=stail
        ftb.field('GCORE')[0][ctidx][enidx]=gcore
        ftb.field('GTAIL')[0][ctidx][enidx]=gtail
        ftb.field('NCORE')[0][ctidx][enidx]=ncore
        ftb.field('NTAIL')[0][ctidx][enidx]=ntail

    #update the back CALDB file entries
    for ffil in backfiles:
        tf = open(ffil)
        idx = int(ffil.split('_')[0])
        enidx=idx%dim1
        ctidx=idx/dim1
        line = tf.readline()
        """sigma,gamma=line.split(' ')
        sigma,gamma=float(sigma),float(gamma)
        #print idx,sigma,gamma
        btb.field('SCORE')[0][idx]=sigma
        btb.field('STAIL')[0][idx]=sigma
        btb.field('GCORE')[0][idx]=gamma
        btb.field('GTAIL')[0][idx]=gamma
        btb.field('NTAIL')[0][idx]=0"""
        lines = line.split(' ')
        lines = [float(line) for line in lines]
        #print idx,sigma,gamma
        if lines[1]>lines[3]:
            score,gcore,stail,gtail,ncore=lines
            ntail= 1-ncore
        else:
            stail,gtail,score,gcore,ntail=lines
            ncore= 1-ntail
        ncore = max(min(1,ncore),1e-4)
        ntail = max((1-ncore)/ncore*(score/stail)**2,0)
        print enidx,ctidx,score,gcore,stail,gtail,ncore
        btb.field('SCORE')[0][ctidx][enidx]=score
        btb.field('STAIL')[0][ctidx][enidx]=stail
        btb.field('GCORE')[0][ctidx][enidx]=gcore
        btb.field('GTAIL')[0][ctidx][enidx]=gtail
        btb.field('NCORE')[0][ctidx][enidx]=ncore
        btb.field('NTAIL')[0][ctidx][enidx]=ntail

    #output the tables and copy into irfs
    ff[1].data=ftb
    bf[1].data=btb
    ff.flush()
    bf.flush()
    pname = os.environ['CALDB']+r'/data/glast/lat/bcf/psf/psf_%s_'%out
    cmd1 = r'cp %s/psf_%s_front.fits %sfront.fits '%(cwd,out,pname)
    cmd2 = r'cp %s/psf_%s_back.fits %sback.fits'%(cwd,out,pname)
    os.system(cmd1)
    os.system(cmd2)
    #plot the new PSF
    #makeplots(['P7SOURCE_V4MC',irf,out],1.1,4,cts=[68.,95.])

def makeallpsftasks():
    mcff = pf.open(os.environ['CALDB']+r'/data/glast/lat/bcf/psf/psf_P7SOURCE_V4MC_front.fits')
    ctbins = len(mcff[1].data.field('CTHETA_LO')[0])
    ebins = len(mcff[1].data.field('ENERG_LO')[0])
    tasks = ['uf.likelihoodtable(%d,%d,%d,0.15)'%(enr,ctr,y) for enr in range(ebins) for ctr in range(ctbins) for y in range(2)]
    return tasks

## run the full likelihood analysis on different cores and output results
def psftable():

    #names of machines to run ipengines on
    #machines = 'tev01 tev02 tev03 tev04 tev05 tev06 tev07 tev08 tev09 tev10 tev11'.split()

    #make sure the pickle objects are made
    likelihoodtable(0,0,0,jg=True)
    likelihoodtable(0,0,1,jg=True)

    #setup the engines to run a table entry on each
    setup_string = 'import uw.stacklike.fitsmooth as uf;reload(uf);from uw.stacklike.fitsmooth import *'
    mcff = pf.open(os.environ['CALDB']+r'/data/glast/lat/bcf/psf/psf_P7SOURCE_V4MC_front.fits')

    from IPython.parallel import Client
    clfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/allhalos.txt','w')
    clfile.close()
    
    tasks = makeallpsftasks()
    
    tpr = len(tasks)
    rc = Client()
    dview = rc[:]

    dview.execute(setup_string)
    t.sleep(20)
    dview.execute('tasks = uf.makeallpsftasks()')
    t.sleep(20)
    dview.scatter('nums',range(len(tasks)))
    t.sleep(20)
    cmd = 'x=[eval(tasks[y]) for y in nums]'
    print cmd
    result = dview.execute(cmd)
    for engine in result.metadata:
        print engine.status
    logfile = open('/phys/groups/tev/scratch1/users/Fermi/mar0/python/mec.log','w')

    #makelikeirfs(irf='P6_v7_diff',out='P6_v12_diff')

def plot_params():
    ffiles = glob.glob('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/*_0.txt')
    bfiles = glob.glob('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/*_1.txt')
    frn = range(len(ffiles))

    ffiles = np.array(['/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/%d_0.txt'%ran for ran in frn])
    bfiles = np.array(['/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/%d_1.txt'%ran for ran in frn])
    enr = np.arange(0,18,1)
    ctr = np.arange(0,8,1)
    fpars = []
    bpars = []
    for fil in ffiles:
        fpars.append(np.loadtxt(fil))
    for fil in bfiles:
        bpars.append(np.loadtxt(fil))
    fpars = np.array(fpars).transpose()
    bpars = np.array(bpars).transpose()
    mi=[0.5,1.,0.5,1.,0.,0.]
    ma=[1.5,3.,1.5,3.,1.,1.]

    cts = [0.34,0.68,0.90,0.95,0.99]

    py.figure(figsize =(36,9))
    for it,ct in enumerate(cts):
        py.subplot(2,3,it+1)
        par = np.array([fcontain(ct,fpars[0][it2],fpars[1][it2],fpars[2][it2],fpars[3][it2],fpars[4][it2]) for it2 in range(len(fpars[0]))])
        py.pcolor(enr,ctr,par.reshape(-1,18),cmap=py.cm.hot)#,vmin=0.,vmax=10)
        py.colorbar()
        py.xlabel(r'$\rm{Energy \/Bin}$')
        py.ylabel(r'$cos(\theta)\/\rm{bin}$')
        py.title(ct)
        #py.xlim(0,17)
        #py.ylim(0,7)
    py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/pars_f.png')
    py.clf()
    for it,ct in enumerate(cts):
        py.subplot(2,3,it+1)
        par = np.array([fcontain(ct,bpars[0][it2],bpars[1][it2],bpars[2][it2],bpars[3][it2],bpars[4][it2]) for it2 in range(len(bpars[0]))])
        py.pcolor(enr,ctr,par.reshape(-1,18),cmap=py.cm.hot)#,vmin=0.,vmax=10)
        py.colorbar()
        py.xlabel(r'$\rm{Energy \/Bin}$')
        py.ylabel(r'$cos(\theta)\/\rm{bin}$')
        py.title(ct)
        #py.xlim(0,17)
        #py.ylim(0,7)
    py.savefig('/phys/groups/tev/scratch1/users/Fermi/mar0/figures/psftable/pars_b.png')

def fcontain(frac,s1,g1,s2,g2,n1):
    n1 = min(n1,1.0)
    n2 = 1-n1
    fint = doublepsfint(inf,s1,g1,n1,s2,g2,n2)
    xr = np.logspace(-2,1.5,1001)
    #print xr
    yr = np.array([doublepsfint(x,s1,g1,n1,s2,g2,n2) for x in xr])/fint
    #print yr
    fnc = intp.UnivariateSpline(np.log(yr),np.log(xr),s=0)
    ret = fnc(np.log(frac))
    #print frac,np.exp(ret)
    #for it in range(len(yr[:-1])):
    #    print 0.5*(yr[it]+yr[it+1]),0.5*(xr[it]+xr[it+1]),np.exp(fnc(np.log(0.5*(yr[it]+yr[it+1]))))[0]
    #print fnc1
    return np.exp(ret)

def setup_engines():
    import subprocess
    machines=['tev01','tev02','tev04','tev05','tev06','tev08','tev09','tev10','tev11']#'tev03','tev07',
    uses = []

    for it,machine in enumerate(machines):
        if os.system('ping %s -c 2 -i 0.5'%machine)==0:
            uses.append(it)
    subp = []
    #cmd = 'ipcontroller --ip=\'*\''
    #subp.append(subprocess.Popen(cmd , shell=True, stdout=subprocess.PIPE).communicate())
    for it in uses:
        for engines in range(16):
            cmd = 'ssh -n -f %s "sh -c \\" nohup ipengine > /dev/null 2>&1 &\\""'%machines[it]
            subp.append(subprocess.Popen(cmd , shell=True, stdout=subprocess.PIPE).communicate())
            t.sleep(0.5)

def kill_engines():
    machines=['tev01','tev02','tev04','tev05','tev06','tev08','tev09','tev10','tev11']#'tev03','tev07',
    uses = []
    for it,machine in enumerate(machines):
        if os.system('ping %s -c 2 -i 0.5'%machine)==0:
            uses.append(it)
    for it in uses:
        os.system('ssh %s killall ipengine -9'%machine)
        os.system('ssh %s killall ipcontroller -9'%machine)


if __name__=='__main__':
    psftable(sys.argv[0])