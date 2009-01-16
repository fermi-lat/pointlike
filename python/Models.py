"""A set of classes to implement spectral models.

   $Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/Models.py,v 1.13 2008/12/15 22:12:26 burnett Exp $

   author: Matthew Kerr

"""

import numpy as N
import math as M

#===============================================================================================#

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

#===============================================================================================#

class DefaultModelValues(object):
   """Static methods and class members to assign default values to the spectral models."""

   simple_models = {
      'PowerLaw'        : {'p':[-11.,2],'param_names':['Norm','Index']},
      'PowerLawFlux'    : {'p':[-7.,2],'param_names':['Int. Flux','Index'],'emin':100,'emax':1e6},
      'BrokenPowerLaw'  : {'p':[-11.,1.5,2.5,1000],'param_names':['Norm','Index 1','Index 2', 'E_break']},
      'BrokenPowerLawF' : {'p':[-11.,1.5,2.5],'param_names':['Norm','Index 1','Index 2'],'e_break':1000},
      'LogParabola'     : {'p':[-11,2.,0.,1000.],'param_names':['Norm','Index','beta','E_break']},
      'ExpCutoff'       : {'p':[-11.,2,5e3],'param_names':['Norm','Index','Cutoff']},
      'PLSuperExpCutoff': {'p':[-11.,2,2e3,2.],'param_names':['Norm','Index','Cutoff', 'b']},
      'PLSuperExpCutoffF': {'p':[-11.,2,2e3],'param_names':['Norm','Index','Cutoff'],'b':2},
      'Constant'        : {'p':[0.],'param_names':['Scale']}
      }

   names = simple_models.keys()+['MixedModel']

   @staticmethod
   def setup(the_model,**kwargs):
      """Pass a model instance to give it default values.  The keyword arguments are used
         only for MixedModel, in which case they contain the 'simple_models' keyword argument to
         describe which simple models the MixedModel comprises."""
      
      DefaultModelValues.start(the_model)
      classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__
      
      if classname in DefaultModelValues.simple_models:
         for key,val in DefaultModelValues.simple_models[classname].items():
            exec('the_model.%s = val'%key)
      
      if classname == 'MixedModel':
         val = kwargs['simple_models'] if 'simple_models' in kwargs else ['PowerLaw']
         if type(val) == type(dict()): #a dictionary with kwargs for each simple model has been passed
            the_model.models = val.keys()
            default_dicts = [val[model] for model in the_model.models]
         else: #a list has been passed; use default values
            the_model.models = val
            default_dicts = [DefaultModelValues.simple_models[model] for model in the_model.models]
         the_model.param_names,the_model.p,the_model.spec_models=[],[],[]
         for i,model in enumerate(the_model.models):
               exec('this_model = %s(**default_dicts[i])'%model)
               the_model.spec_models += [this_model]
               the_model.param_names += this_model.param_names
               the_model.p += list(this_model.p)
         the_model.n = [len(x.p) for x in the_model.spec_models]
         the_model.pretty_name = '+'.join(the_model.models)

      DefaultModelValues.finish(the_model)

   @staticmethod
   def start(the_model):
      """Common values independent of the model type."""
      the_model.e0 = 1000.
      the_model.flux_scale = 1.
      the_model.good_fit = False
      the_model.p = None
      the_model.param_names = ['Undefined Model']

   @staticmethod
   def finish(the_model):
      """Common values that can be written once the model type has been sussed out."""
      the_model.cov_matrix = N.zeros([len(the_model.p),len(the_model.p)]) #default covariance matrix
      the_model.free = N.asarray([True] * len(the_model.p))
      the_model.p = N.asarray(the_model.p)
      

#===============================================================================================#

class Model(object):
   """Spectral model giving dN/dE for a point source.  Default units are ph/cm^2/s/MeV."""
   def __init__(self,**kwargs):
      """
Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  e0          [1000] value in MeV at which normalization is given
  flux_scale  [1e7] multiplier for actual value to make units more convenient
  p           [p1,p2,...] default values of spectral parameters; see docstring individual model classes
  models      [] for MixedModel, a list of simple model names composing MixedModel
  free        [True, True,...] a boolean list the same length as p giving the free (True) and fixed (False) parameters
  =========   =======================================================
     """
      DefaultModelValues.setup(self,**kwargs)
      if 'p' in kwargs: kwargs['p'][0] = N.log10(kwargs['p'][0]) #convert norm to log
      self.__dict__.update(**kwargs)

   def get_parameters(self):
      """Return FREE parameters; used for spectral fitting."""
      return self.p[self.free]

   def set_parameters(self,new_vals):
      """Set FREE parameters; new_vals should have length equal to number of free parameters."""
      self.p[self.free] = new_vals

   def freeze(self,parameter,freeze=True):
      """Freeze one of the spectral parameters from fitting.
      
         parameter: a parameter name or index.
         freeze   : if True, freeze parameter; if False, free it
         """
      if type(parameter) == type(''):
         for n,name in enumerate(self.param_names):
            if parameter == name: parameter = n; break
      self.free[parameter] = not freeze

   def set_cov_matrix(self,new_cov_matrix):
      #if len(new_cov_matrix) == 1: new_cov_matrix = new_cov_matrix.flatten()
      self.cov_matrix[N.outer(self.free,self.free)] = N.ravel(new_cov_matrix)

   def statistical(self):
      """Return the parameter values and fractional statistical errors.
         If no error estimates are present, return 0 for the fractional error."""
      p=N.array(self.p).copy()
      p[0] = 10**p[0] #convert from log measurement
      try: #see if error estimates are present
         errors=N.diag(self.cov_matrix)**0.5
         ratios=errors/p
         ratios[0] = (10**errors[0]-10**-errors[0])/2 #kludge to take care of logarithm
         return p,ratios
      except:
         return p,N.zeros_like(p)

   def systematic(self,ratios):
      try: #see if systematic errors are present
         s = self.systematics
         up_sys_errs = self.systematics[:,1]
         low_sys_errs = self.systematics[:,0]
         up_sys_errs = N.nan_to_num((up_sys_errs**2 - ratios**2)**0.5)
         low_sys_errs = N.nan_to_num((low_sys_errs**2 - ratios**2)**0.5)
         return up_sys_errs,low_sys_errs
      except AttributeError:
         return N.zeros_like(ratios),N.zeros_like(ratios)

   def __str__(self):
      """Return a pretty print version of parameter values and errors."""
      p,ratios = self.statistical()
      m=max([len(n) for n in self.param_names])
      l=[]
      if ratios.sum() > 0: #if statistical errors are present
         up_sys_errs,low_sys_errs = self.systematic(ratios)    
         for i in xrange(len(self.param_names)):
            n=self.param_names[i][:m]
            t_n=n+(m-len(n))*' '
            frozen = '' if self.free[i] else '(FROZEN)'
            if up_sys_errs.sum()==0: #if no systematic errors present
               l+=[t_n+': (1 +/- %.3f) %.3g %s'%(ratios[i],p[i],frozen)]
            else:
               l+=[t_n+': (1 + %.3f - %.3f)(1 +/- %.3f) %.3g %s'%(up_sys_errs[i],low_sys_errs[i],ratios[i],p[i],frozen)]
         return '\n'.join(l)
      else: #if no errors are present
         for i in xrange(len(self.param_names)):
            n=self.param_names[i][:m]
            t_n=n+(m-len(n))*' '
            l+=[t_n+': %.3g'%(p[i])]
         return '\n'.join(l)

   def i_flux(self,emin=100,emax=1e6,e_weight=0,cgs=False,error=False):
      """Return the integral flux.         

Optional keyword arguments:

  =========   =======================================================
  Keyword     Description
  =========   =======================================================
  emin        [100] lower bound in MeV
  emax        [1e6] upper bound in MeV
  e_weight    [0] energy power by which to scale dN/dE
  cgs         [False] if True, energy in ergs
  error       [False] if True, return value is a tuple with flux and estimated error
  =========   =======================================================
      """
      try:
         from scipy.integrate import quad
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
      except:
         print 'Encountered a numerical error when attempting to calculate integral flux.'

   def copy(self):
      self.p[0] = 10**self.p[0] #account for logarithm -- will be restored by next line
      a = eval(self.name+'(**self.__dict__)')
      a.p = [x for x in self.p]
      try: a.cov_matrix = self.cov_matrix.__copy__()
      except: pass
      return a

   def fast_iflux(self,emin=100,emax=1e6):
      """Return a quick calculation for photon flux for models where it is analytically available."""
      return self.i_flux(emin=emin,emax=emax)

   def __flux_derivs__(self,*args):
      """Use finite differences to estimate the gradient of the integral flux wrt the spectral parameters."""
      delta = .02
      hi,lo = self.copy(),self.copy()
      derivs = []
      for i in xrange(len(self.p)):
         my_delta = delta if self.p[i] >= 0 else -delta
         hi.p[i]*=(1+my_delta/2.)
         lo.p[i]*=(1-my_delta/2.)
         derivs += [(hi.i_flux(*args) - lo.i_flux(*args))/(delta*self.p[i])]
         hi.p[i]/=(1+my_delta/2.)
         lo.p[i]/=(1-my_delta/2.)

      return N.asarray(derivs)

#===============================================================================================#

class PowerLaw(Model):
   """Implement a power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  gamma      (absolute value of) spectral index
      """
   def __call__(self,e):
      n0,gamma=self.p
      return (10**n0/self.flux_scale)*(self.e0/e)**gamma

   def fast_iflux(self,emin=100,emax=1e6):
      n0,gamma = self.p
      return 10**n0/(1-gamma)*self.e0**gamma*(emax**(1-gamma)-emin**(1-gamma))

#===============================================================================================#

class PowerLawFlux(Model):
   """Implement a power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  flux       integrated flux from emin to emax MeV
  gamma      (absolute value of) spectral index
      """
   def __call__(self,e):
      flux,gamma=self.p
      #return (flux/self.flux_scale)*(gamma-1)*(self.e0/e)**(gamma-1)/e
      return (10**flux/self.flux_scale)*(1-gamma)*e**(-gamma)/(self.emax**(1-gamma)-self.emin**(1-gamma))

   def fast_iflux(self,emin=100,emax=1e6):
      n0,gamma = self.p
      return 10**n0*(emax**(1-gamma) - emin**(1-gamma)) / (self.emax**(1-gamma) - self.emin**(1-gamma))

#===============================================================================================#

class BrokenPowerLaw(Model):
   """Implement a broken power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  gamma1     (absolute value of) spectral index for e < e_break
  gamma2     (absolute value of) spectral index for e > e_break
  e_break    break energy (free parameter)
      """
   def __call__(self,e):
      e = N.array([e]).ravel()
      n0,gamma1,gamma2,e_break=self.p
      return (10**n0/self.flux_scale)*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

#===============================================================================================#

class BrokenPowerLawF(Model):
   """Implement a broken power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  gamma1     (absolute value of) spectral index for e < e_break
  gamma2     (absolute value of) spectral index for e > e_break
  e_break    break energy (nota bene: fixed!)
      """
   def __call__(self,e):
      e = N.array([e]).ravel()
      n0,gamma1,gamma2=self.p
      e_break = self.e_break
      return (10**n0/self.flux_scale)*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

#===============================================================================================#

class LogParabola(Model):
   """Implement a log parabola (for blazars.)  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  alpha      (absolute value of) constant spectral index
  beta       co-efficient for energy-dependent spectral index
  e_break    break energy
      """
   def __call__(self,e):
      n0,alpha,beta,e_break=self.p
      return (10**n0/self.flux_scale)*(e_break/e)**(alpha - beta*N.log(e_break/e))

#===============================================================================================#

#===============================================================================================#

class ExpCutoff(Model):
   """Implement a power law with exponential cutoff.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  gamma      (absolute value of) spectral index
  cutoff     e-folding cutoff energy (MeV)
      """
   def __call__(self,e):
      n0,gamma,cutoff=self.p
      if cutoff < 0: return 0
      return (10**n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-e/cutoff)

#===============================================================================================#

class PLSuperExpCutoff(Model):
   """Implement a power law with hyper-exponential cutoff.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  gamma      (absolute value of) spectral index
  cutoff     e-folding cutoff energy (MeV)
  b          additional power in the exponent
      """
   def __call__(self,e):
      n0,gamma,cutoff,b=self.p
      if cutoff < 0: return 0
      return (10**n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-(e/cutoff)**b)

#===============================================================================================#

class PLSuperExpCutoffF(Model):
   """Implement a power law with fixed hyper-exponential cutoff.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  gamma      (absolute value of) spectral index
  cutoff     e-folding cutoff energy (MeV)
  b          additional power in the exponent
      """
   def __call__(self,e):
      n0,gamma,cutoff=self.p
      if cutoff < 0: return 0
      return (10**n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-(e/cutoff)**self.b)

#===============================================================================================#

class MixedModel(Model):
   """Implement a composite model.  The value is the sum of the simple models.
      See constructor docstring for further keyword arguments.
      NOTA BENE: specify the simple models via the keyword arguments 'models'
      """
   def __call__(self,e):
      counter = 0
      for i in xrange(len(self.n)):
         self.spec_models[i].p = self.p[counter:counter+self.n[i]]
         counter += self.n[i]
      return N.array([model(e) for model in self.spec_models]).sum(axis=0)


#===============================================================================================#

class Constant(Model):
   def __call__(self,e):
      return N.ones_like(e)*10**self.p[0]
   
   def fast_iflux(self,emin=100,emax=1e6):
      return (emax-emin)*10**self.p[0]

#===============================================================================================#
if __name__=='__main__':
   disp=Dispersion()
   print disp(300,330)