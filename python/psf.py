import numpy as N
from math import log

class PSF:

   def __init__(self,sigma_parameters=None, gamma_parameters=None):

      if sigma_parameters is None:
         self.fp=N.array([0.0166,-1.55,1.69,1.34])#allgammav15r21
         self.bp=N.array([0.0373,-1.72,6.13,4.91])#allgammav15r21
         #self.fp=N.array([0.060,-1.422,1.532,0.966])
         #self.bp=N.array([0.142,-1.65,3.881,21.671])#MET=243873727 (80 days)

      #ignore user-supplied gammas for now
      self.gm = GammaManager()

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
      
   def set_sigmas(self):
      import pointlike as pl
      pl.IParams.set_fp(self.fp)
      pl.IParams.set_bp(self.bp)

   def fmax(self,e,u,event_class=0):
      g = self.gamma(e,event_class)
      return 1-(1+u/g)**(1-g)
      

class GammaManager(object):
   def __init__(self):

      try:
         #f = open('boogie woogie')#just get an exception!
         path = 'd:/common/sourcelikefilter/first_light/psf/day1-80'
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