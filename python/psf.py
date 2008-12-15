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

   def init(self):
      self.use_mc = False
      self.umax = 50

   def sigma(self,e,event_class=0):
      a,b,c,d = self.fp if event_class == 0 else self.bp
      return ( a**2 + c*(e/100)**b + d*(e/100)**(2*b) )**0.5

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