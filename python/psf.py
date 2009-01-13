import numpy as N
from math import log

class PSF:

   def __init__(self,**kwargs):

      self.init()
      self.__dict__.update(kwargs)       

      if self.use_mc:
         self.fp=N.array([0.014,-1.400,1.056,1.434])#allgammav15r0
         self.bp=N.array([0.033,-1.734,6.270,0.599])
      else:
         self.fp = N.array([0.045,-1.303,.767,2.360]) #day-126
         self.bp = N.array([-0.123,-1.302,1.094,15.109])
         #self.fp = N.array([0.063,-1.673,1.850,-0.039]) #p6cuts
         #self.bp = N.array([-0.112,-1.561,3.586,12.753])
         
      self.gm = GammaManager(use_mc=self.use_mc)

      if self.use_psf_init:

         psf_init = PSFInit()
         self.fp = psf_init.fparams
         self.bp = psf_init.bparams
         self.gm.f_elist = psf_init.energy
         self.gm.b_elist = psf_init.energy
         self.gm.f_gamma = psf_init.fgam
         self.gm.b_gamma = psf_init.bgam

   def init(self):
      self.use_mc = False
      self.use_psf_init = False
      self.umax = 50

   def sigma(self,e,event_class=0):
      a,b,c,d = self.fp if event_class == 0 else self.bp
      return ( a**2 + c*(e/100.)**b + d*(e/100.)**(2*b) )**0.5

   def gamma(self,e,event_class=0):
      return self.gm.gamma(e,event_class)

   def conf_region(self,e,percentage=68,event_class = 0,use_u = False):
      "Return the _percentage_ confidence region radius for the provided energy.  Unit is degrees."""
      gamma = self.gamma(e,event_class)
      u = gamma*((1-percentage/100.)**(1./(1-gamma))-1)
      if use_u: return u
      return self.u2deg(e,u,event_class)

   def u2deg(self,e,u,event_class=0):
      return self.sigma(e,event_class)*(2*u)**0.5

   def deg2u(self,e,deg,event_class=0):
      return 0.5*(deg/self.sigma(e,event_class))**2
      
   #def set_sigmas(self):
   #   from skymaps import IParams
   #   IParams.set_fp(self.fp)
   #   IParams.set_bp(self.bp)

   def set_psf_params(self,psl):
      for s in psl:
         e = (s.band().emin()*s.band().emax())**0.5
         s.setgamma(self.gamma(e,s.band().event_class()))
         s.setsigma(self.sigma(e,s.band().event_class())/180*N.pi)
         s.recalc()

   def fmax(self,e,u,event_class=0):
      g = self.gamma(e,event_class)
      return 1-(1+u/g)**(1-g)

   def __call__(self,e,delta,event_class=0):
      delta = delta * 180 / N.pi
      u = (delta / self.sigma(e,event_class))**2/2
      gam = self.gamma(e,event_class)
      fmax = 1-(1+u/gam)**(1-gam)
      return 1./fmax*(1-1/gam)*(1+u/gam)**-gam
      

class GammaManager(object):
   def __init__(self,use_mc=False):

      try:
         #f = open('boogie woogie')#just get an exception!
         path = 'd:/common/sourcelikefilter/first_light/psf/126_day'
         if use_mc:
            files = [path+'/m_front.txt',path+'/m_back.txt']
         else:
            files = [path+'/front.txt',path+'/back.txt']
         for i in xrange(len(files)): self.__parse_file__(files[i],i)

      except:         
         print 'Unable to open PSF files -- reverting to hardcoded defaults.'
         #Defaults if can't find files
         self.f_elist = self.b_elist = 100*10**N.arange(0,2.01,1./10) #10 per decade in this analysis
         self.f_gamma = [2.42,2.47,2.84,2.07,2.08,2.06,2.03,1.84,1.90,1.81,
                         1.85,1.88,1.82,1.97,1.94,1.86,1.89,1.90,1.89,1.95,1.84]
         self.b_gamma = [3.46,2.24,2.77,2.38,1.86,1.86,1.93,1.82,1.81,1.79,
                         2.16,1.73,1.83,1.79,1.94,1.76,1.73,1.74,1.68] + [1.68]*2 #Just replicate last value

   def gamma(self,e,event_class=0):
      ens = self.f_elist if event_class == 0 else self.b_elist
      g_list = self.f_gamma if event_class == 0 else self.b_gamma
      if e < ens[0]: return g_list[0]
      if e > ens[-1]: return g_list[-1]
      for i in xrange(len(ens)-1):
         if ens[i]<=e and ens[i+1]>e: break
      e1,e2,f1,f2 = ens[i],ens[i+1],g_list[i],g_list[i+1]
      return f1*(e/e1)**(log(f2/f1)/log(e2/e1)) 

   def __parse_file__(self,filename,event_class):

      f = open(filename)
      g = f.readlines()
      f.close()
      g = N.array([x.strip().split() for x in g]).astype(float)
      ens = g[:,0]
      gams = g[:,3]

      letter = 'f' if event_class==0 else 'b'
      exec('self.%s_elist = ens'%letter)
      exec('self.%s_gamma = gams'%letter)

if __name__=='__main__':
   #Testing script
   psf = PSF()
   ens = [100*10**(.2*x) for x in xrange(10)]
   regs = [psf.conf_region(x) for x in ens]
   print 'Events converted in thin-foil TKR, 68 percent'
   for i in xrange(len(ens)):
      print 'Energy: %d MeV  Containment: %.2f deg'%(ens[i],regs[i])
   regs = [psf.conf_region(x,event_class=1) for x in ens]
   print '\nEvents converted in thick-foil TKR, 68 percet'
   for i in xrange(len(ens)):
      print 'Energy: %d MeV  Containment: %.2f deg'%(ens[i],regs[i])


"""
#### defaults for setting up psf parameters #####
Author: M. Roth
$Header: /nfs/slac/g/glast/ground/cvs/skymaps/python/psf_defaults.py,v 1.2 2008/11/10 22:37:24 burnett Exp $
"""
import pyfits as pf
import numpy as N
import pylab as p
import scipy.optimize

mf = 1.5
mb = 1.2

"""---------------------------Helper functions-------------------------"""
#psf scale function from handoff response
def scale(e,s,p):
	pa=[0,0]

	#front
	if p==0:
		pa[0] = 0.058
		pa[1] = 0.000377
	#back
	if p==1:
		pa[0] = 0.096
		pa[1] = 0.0013
	rs=[]
	for i in N.arange(0,len(s),1):
		rs.append(s[i]*N.sqrt(pa[0]*pa[0]*(e[i]/100)**(-1.6)+pa[1]*pa[1])*180/N.pi)
	return rs

#calculate chi-squared
def chi(a,b,c,d,sarr,Earr):
	logl=0
	for i in N.arange(0,len(sarr)):
		logl+=(sarr[i]-func(a,b,c,d,Earr[i]))**2/(sarr[i]*0.07)**2
	return logl

#fit function
def func(a,b,c,d,le):
	ssq = a**2+c*(le/100)**(b)+d*(le/100)**(2*b)
	if ssq<0 or c<0:
		return 1e80
	return N.sqrt(ssq)

class PSFInit(object):
   """----------------------------------------------------------------------------------------"""


   """-------------------------------FITS table lookup----------------------------------------"""

   def __init__(self):
      path=r'F:/glast/caldb/v0r7p1/CALDB/data/glast/lat/bcf/psf/'

      #names of front and back CALDB files
      ffile=path+'psf_P6_v1_diff_front.fits'
      bfile=path+'psf_P6_v1_diff_back.fits'

      #check for environment variable
      import os
      if 'CALDB' in os.environ and ffile=='' and bfile=='':
         ffile = os.path.join(os.environ['CALDB'],'v0r7p1','CALDB','data', 'glast', 'lat', 'bsf', 'psf','psf_P6_v1_diff_front.fits')
         bfile = os.path.join(os.environ['CALDB'],'v0r7p1','CALDB','data', 'glast', 'lat', 'bsf', 'bsf','psf_P6_v1_diff_back.fits')

      #open fits files and point to tables
      frontfile = pf.open(ffile)
      backfile  = pf.open(bfile)
      fpsftable = frontfile[1].data
      bpsftable = backfile[1].data

      #x corresponds to energies, y to cos theta
      xit = N.arange(0,18,1)
      yit = N.arange(0,8,1)
      sf=[];sb=[]
      gf=[];gb=[]

      #go through energies
      for i in xit:
         w=[0,0]
         s=[0,0]
         g=[0,0]
         for j in yit:

            #weight contribution for sigma and gamma based on effective area
            #slopes in cos theta are defined at the top of the file as 'mf' and 'mb' for front/back
            costh=(j-7)/10.-0.5
            wf=mf*costh+1
            if wf<0:
               wf=0
            wb=mb*costh+1
            if wb<0:
               wb=0
            w[0]=w[0]+wf
            w[1]=w[1]+wb
            s[0]=s[0]+wf*fpsftable.field('SIGMA')[0][18*j+i]
            s[1]=s[1]+wb*bpsftable.field('SIGMA')[0][18*j+i]

            #weight gamma by core photon fraction
            ncore = fpsftable.field('NCORE')[0][18*j+i]
            if ncore>1:
               ncore=1
            if ncore<0:
               ncore=0
            g[0]=g[0]+wf*ncore*fpsftable.field('GCORE')[0][18*j+i]
            g[0]=g[0]+wf*(1-ncore)*fpsftable.field('GTAIL')[0][18*j+i]
            ncore = bpsftable.field('NCORE')[0][18*j+i]
            if ncore>1:
               ncore=1
            if ncore<0:
               ncore=0
            g[1]=g[1]+wb*ncore*bpsftable.field('GCORE')[0][18*j+i]
            g[1]=g[1]+wb*(1-ncore)*bpsftable.field('GTAIL')[0][18*j+i]
         sf.append(s[0]/w[0])
         sb.append(s[1]/w[1])
         gf.append(g[0]/w[0])
         gb.append(g[1]/w[1])

      #setup energy bins as 10/decade starting with 17 MeV (same as CALDB)
      energy = 10**(1.25+0.25*xit)

      #scale sigma values using PSF scale function from handoff-response
      sf = scale(energy,sf,0)
      sb = scale(energy,sb,1)

      #fit sigma for IParams class
      self.fparams = N.array(scipy.optimize.fmin_powell(lambda x: chi(x[0],x[1],x[2],x[3],sf,energy),(0.01,-1.43,1.46,2.43), ftol = 0.0001,full_output=1, disp=0)[0])
      self.bparams = N.array(scipy.optimize.fmin_powell(lambda x: chi(x[0],x[1],x[2],x[3],sb,energy),(0.01,-1.43,1.46,2.43), ftol = 0.0001,full_output=1, disp=0)[0])

      #set gammas
      self.fgam = gf
      self.bgam = gb

      self.energy = energy