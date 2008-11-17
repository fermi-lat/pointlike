"""A set of classes to implement spectral models.

"""

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

class DefaultModelValues(object):

   names = ['PowerLaw', 'PowerLawFlux', 'BrokenPowerLaw', 'BrokenPowerLawF','ExpCutoff',
            'PLSuperExpCutoff']

   @staticmethod
   def setup(the_model,**kwargs):
      the_model.e0 = 1000.
      the_model.flux_scale = 1e7
      the_model.good_fit = False
      classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__
      if classname == 'PowerLaw':
         the_model.p = [1e-4,2]
         the_model.param_names = ['Norm','Index']
      elif classname == 'PowerLawFlux':
         the_model.p = [1e-3,2]
         the_model.param_names = ['Int. Flux','Index']
         the_model.emin = 100.
         the_model.emax = 1e6
      elif classname == 'BrokenPowerLaw':
         the_model.p = [1e-4,1.5,2.5,1000]
         the_model.param_names = ['Norm','Index 1','Index 2', 'E_break']
      elif classname == 'BrokenPowerLawF':
         the_model.p = [1e-4,1.5,2.5,1000]
         the_model.param_names = ['Norm','Index 1','Index 2']
         the_model.e_break = 10000
      elif classname == 'ExpCutoff':
         the_model.p = [1e-4,2,5e3]
         the_model.param_names = ['Norm','Index','Cutoff']
      elif classname == 'PLSuperExpCutoff':
         the_model.p = [1e-4,2,2e3,2]
         the_model.param_names = ['Norm','Index','Cutoff', 'b']
      elif classname == 'MixedModel':
         #model names are in kwargs
         if 'models' in kwargs:
            models = kwargs['models']
            kwargs.pop('models')
            the_model.models = models
         else:
            the_model.p = None
            the_model.param_names = ['Undefined Model']
            models = []
         the_model.param_names,the_model.p,the_model.spec_models=[],[],[]
         for model in models:
            exec('this_model = %s(**kwargs)'%model)
            DefaultModelValues.setup(this_model)
            the_model.spec_models += [this_model]
            the_model.param_names += this_model.param_names
            the_model.p += this_model.p
         the_model.n = [len(x.p) for x in the_model.spec_models]
         the_model.pretty_name = '+'.join(models)
    
      else:
         the_model.p = None
         the_model.param_names = ['Undefined Model']
      

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class Model(object):
   def __init__(self,**kwargs):
      DefaultModelValues.setup(self,**kwargs)
      self.__dict__.update(**kwargs)
   def __str__(self):
      p=N.array(self.p)
      m=max([len(n) for n in self.param_names])
      l=[]
      try:
         errors=N.diag(self.cov_matrix)**0.5
         ratios=errors/p         
         try:
            s = self.systematics
            sys_flag = True
            up_sys_errs = self.systematics[:,1]
            low_sys_errs = self.systematics[:,0]
            print up_sys_errs
            print low_sys_errs
            up_sys_errs = N.nan_to_num((up_sys_errs**2 - ratios**2)**0.5)
            print up_sys_errs
            low_sys_errs = N.nan_to_num((low_sys_errs**2 - ratios**2)**0.5)
            print low_sys_errs
         except AttributeError:
            sys_flag = False         
         for i in xrange(len(self.param_names)):
            n=self.param_names[i][:m]
            t_n=n+(m-len(n))*' '
            if not sys_flag:
               l+=[t_n+': (1 +/- %.3f) %.3g'%(ratios[i],p[i])]
            else:
               l+=[t_n+': (1 + %.3f - %.3f)(1 +/- %.3f) %.3g'%(up_sys_errs[i],low_sys_errs[i],ratios[i],p[i])]
         return '\n'.join(l)
      except:
         for i in xrange(len(self.param_names)):
            n=self.param_names[i][:m]
            t_n=n+(m-len(n))*' '
            l+=[t_n+': %.3g'%(p[i])]
         return '\n'.join(l)

   def i_flux(self,emin=100,emax=Inf,e_weight=0,cgs=False,error=False):
      """Return integrated flux."""
      if not self.good_fit:
         print 'Model fit is not good, cannot calculate integrated flux.'
         return
      func = self if e_weight == 0 else lambda e: self(e)*e**e_weight
      units = 1.60218e-6**(e_weight) if cgs else 1. #extra factor from integral!
      flux =  units*quad(func,emin,emax)[0]
      if error:
         variances = N.diag(self.cov_matrix)
         args = (emin,emax,e_weight,cgs,False)
         err = (variances*self.__flux_derivs__(*args)**2).sum()**0.5
      if not cgs: flux*=self.flux_scale
      if error: return (flux,err)
      return flux

   def copy(self):
      param_string=','.join( ( str(p) for p in self.p ) )
      a = eval(self.name+'(**self.__dict__)')
      a.p = [x for x in self.p]
      try:
         a.cov_matrix = self.cov_matrix.__copy__()
      except:
         pass
      return a

   def __flux_derivs__(self,*args):
      delta = .02
      hi,lo = self.copy(),self.copy()
      derivs = []
      for i in xrange(len(self.p)):
         hi.p[i]*=(1+delta/2.)
         lo.p[i]*=(1-delta/2.)
         derivs += [(hi.i_flux(*args) - lo.i_flux(*args))/(delta*self.p[i])]
         hi.p[i]/=(1+delta/2.)
         lo.p[i]/=(1-delta/2.)

      return N.asarray(derivs)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ModelManager(object):

   pass

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PowerLaw(Model):
   def __call__(self,e):
      n0,gamma=self.p
      return (n0/self.flux_scale)*(self.e0/e)**gamma

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PowerLawFlux(Model):
   """Use flux (ph cm^-2 s^-1) integrated from emin to emax."""
   def __call__(self,e):
      flux,gamma=self.p
      #return (flux/self.flux_scale)*(gamma-1)*(self.e0/e)**(gamma-1)/e
      return (flux/self.flux_scale)*(1-gamma)*e**(-gamma)/(self.emax**(1-gamma)-self.emin**(1-gamma))

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class BrokenPowerLaw(Model):
   def __call__(self,e):
      e = N.array([e]).ravel()
      n0,gamma1,gamma2,e_break=self.p
      return (n0/self.flux_scale)*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class BrokenPowerLawF(Model):
   def __call__(self,e):
      e = N.array([e]).ravel()
      n0,gamma1,gamma2=self.p
      e_break = self.e_break
      return (n0/self.flux_scale)*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class ExpCutoff(Model):
   def __call__(self,e):
      n0,gamma,cutoff=self.p
      if cutoff < 0: return 0
      return (n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-e/cutoff)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class PLSuperExpCutoff(Model):
   def __call__(self,e):
      n0,gamma,cutoff,b=self.p
      if cutoff < 0: return 0
      return (n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-(e/cutoff)**b)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class MixedModel(Model):
   def __call__(self,e):
      counter = 0
      for i in xrange(len(self.n)):
         self.spec_models[i].p = self.p[counter:counter+self.n[i]]
         counter += self.n[i]
      return N.array([model(e) for model in self.spec_models]).sum(axis=0)


#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#

class Zero(Model):
   def __call__(self,e):
      return N.asarray(0*e)

#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------#
if __name__=='__main__':
   disp=Dispersion()
   print disp(300,330)