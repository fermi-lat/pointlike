import numpy as N
import math as M
from scipy.integrate import quad,Inf


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class Dispersion(object):
   """Manage the energy dispersion pdf."""
   def __init__(self,input_min=18,input_max=700000,output_min=30,output_max=520000):
      self.input_min=input_min #Artificial Monte Carlo bounds
      self.input_max=input_max
      self.output_min=output_min
      self.output_max=output_max
   def __call__(self,e1,e2): #P(e1|e2)
      if (e2<self.input_min or e2>self.input_max): return 0.
      elif (e1<self.output_min or e1>self.output_max): return 0
      var=(e2*0.08+0.6*e2**0.5)**2
      return (2*M.pi*var)**(-0.5)*M.exp(-0.5*(e1-e2)**2/var)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
class Model(object):
   def error_string(self):
      p=N.array(self.p)
      errors=N.diag(self.cov_matrix)**0.5
      ratios=errors/p
      m=max([len(n) for n in self.param_names])
      l=[]
      for i in xrange(len(self.param_names)):
         n=self.param_names[i][:m]
         t_n=n+(m-len(n))*' '
         l+=[t_n+': (1 +/- %.3f) %.3g'%(ratios[i],p[i])]
         #t_p='%.3g'%p[i]
         #t_e='%.3g'%errors[i]
         #l+=[t_n+(12-len(t_n))*' '+': '+t_p+' +/- '+t_e]
         #l+=[t_n+(12-len(t_n))*' '
      return '\n'.join(l)
   def i_flux(self,emin=100,emax=Inf,e_weight=0,cgs=False):
      """Return integrated flux."""
      func = self if e_weight == 0 else lambda e: self(e)*e**e_weight
      units = 1 if not cgs else (e_weight*1.60218e-6)
      return units*quad(func,emin,emax)[0]
   def copy(self):
      param_string=','.join( ( str(p) for p in self.p ) )
      return eval(self.name+'(parameters=('+param_string+'))')



#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PowerLaw(Model):
   def __init__(self,parameters=(1e-9,1.6),e0=100.):
      self.p=parameters
      self.e0=e0
      self.name='PowerLaw'
      self.param_names=['Norm','Index']
   def __call__(self,e):
      n0,gamma=self.p
      return n0*(self.e0/e)**gamma

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PowerLawScaled(Model):
   def __init__(self,parameters=(1.,2),e0=100.):
      self.p=parameters
      self.e0=e0
      self.name='PowerLawScaled'
      self.param_names=['Norm','Index']
   def __call__(self,e):
      n0,gamma=self.p
      return n0/1e9*(self.e0/e)**gamma

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PowerLawFlux(Model):
   """Use flux (ph cm^-2 s^-1) integrated above e0 to normalize."""
   def __init__(self,parameters=(1e-7,1.5),e0=100.):
      self.p=parameters
      self.e0=e0
      self.name='PowerLawFlux'
      self.param_names=['Int. Flux','Index']
   def __call__(self,e):
      flux,gamma=self.p
      return flux*(gamma-1)*(self.e0/e)**(gamma-1)/e

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class BrokenPowerLaw(Model):
   def __init__(self,parameters=(1e-9,1.5,2.5,1000),e0=1):
      self.p=parameters
      self.name='BrokenPowerLaw'
      self.param_names=['Norm','Index 1','Index 2', 'E_break']
      self.e0=e0
   def __call__(self,e):
      e = N.array([e]).ravel()
      n0,gamma1,gamma2,e_break=self.p
      return n0*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class BrokenPowerLawF(Model):
   def __init__(self,parameters=(3.9e-11,1.7,3.2),e_break=1620.):
      self.p=parameters
      self.name='BrokenPowerLawF'
      self.param_names=['Norm','Index 1','Index 2']
      self.e_break = e_break
   def __call__(self,e):
      e = N.array([e]).ravel()
      n0,gamma1,gamma2=self.p
      e_break = self.e_break
      return n0*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ExpCutoff(Model):
   """Is this the right form?"""
   def __init__(self,parameters=(10**-8,2,150000),e_break=1,e0=100.):
      self.p=parameters
      self.e0=e0
      self.e_break=e_break
      self.name='ExpCutoff'
      self.param_names=['Norm','Index','p1']
   def __call__(self,e):
      n0,gamma,p1=self.p
      e_break=self.e_break
      #if e_break<=0 or -e/e_break>700: return 0
      if p1<0: return 0
      return n0*(self.e0/e)**gamma*N.where(e>e_break,N.exp(-(e-e_break)/p1),1.)
      
      #if e > e_break:
      #   return n0*(self.e0/e)**gamma*M.exp(-(e-e_break)/p1)
      #else:
      #   return n0*(self.e0/e)**gamma
      #return n0*(self.e0/e)**gamma*M.exp(-e/e_break)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PLSuperExpCutoff(Model):
   def __init__(self,parameters=(10**-8,2,2000,.5),e0=100):
      self.p=parameters
      self.e0=e0
      self.name='PLSuperExpCutoff'
      self.param_names=['Norm','Index 1','E_cutoff', 'Index 2']
   def __call__(self,e):
      n0,gamma1,e_cutoff,gamma2=self.p
      return n0*(self.e0/e)**gamma1*N.exp(-(e/e_cutoff)**gamma2)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class FixedPLSuperExpCutoff(Model):
   def __init__(self,parameters=(10**-8,2,2000),e0=100):
      self.p=parameters
      self.e0=e0
      self.name='FixedPLSuperExpCutoff'
      self.param_names=['Norm','Index 1','E_cutoff']
   def __call__(self,e):
      n0,gamma1,e_cutoff=self.p
      return n0*(self.e0/e)**gamma1*N.exp(-(e/e_cutoff))

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

if __name__=='__main__':
   disp=Dispersion()
   print disp(300,330)