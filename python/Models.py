"""A set of classes to implement spectral models.

   $Header$

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
      'PowerLaw'        : {'p':[1e-4,2],'param_names':['Norm','Index']},
      'PowerLawFlux'    : {'p':[1.,2],'param_names':['Int. Flux','Index'],'emin':100,'emax':1e6},
      'BrokenPowerLaw'  : {'p':[1e-4,1.5,2.5,1000],'param_names':['Norm','Index 1','Index 2', 'E_break']},
      'BrokenPowerLawF' : {'p':[1e-4,1.5,2.5,1000],'param_names':['Norm','Index 1','Index 2'],'e_break':1000},
      'ExpCutoff'       : {'p':[1e-4,2,5e3],'param_names':['Norm','Index','Cutoff']},
      'PLSuperExpCutoff': {'p':[1e-4,2,2e3,2.],'param_names':['Norm','Index','Cutoff', 'b']} }

   names = simple_models.keys()+['MixedModel']

   @staticmethod
   def setup(the_model,**kwargs):
      """Pass a model instance to give it default values.  The keyword arguments are used
         only for MixedModel, in which case they contain the 'models' keyword argument to
         describe which simple models the MixedModel comprises."""
      
      DefaultModelValues.common_values(the_model)
      classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__
      
      if classname in DefaultModelValues.simple_models:
         for key,val in DefaultModelValues.simple_models[classname].items():
            exec('the_model.%s = val'%key)
      
      if classname == 'MixedModel':         
         the_model.models = kwargs['models'] if 'models' in kwargs else ['PowerLaw'] #defult
         the_model.param_names,the_model.p,the_model.spec_models=[],[],[]
         for model in the_model.models:
            exec('this_model = %s(**kwargs)'%model) #instantiate simple model
            the_model.spec_models += [this_model]
            the_model.param_names += this_model.param_names
            the_model.p += this_model.p
         the_model.n = [len(x.p) for x in the_model.spec_models]
         the_model.pretty_name = '+'.join(the_model.models)
      the_model.cov_matrix = N.zeros([len(the_model.p),len(the_model.p)]) #default covariance matrix
      the_model.free = N.asarray([True] * len(the_model.p))
      the_model.p = N.asarray(the_model.p)

   @staticmethod
   def common_values(the_model):
      the_model.e0 = 1000.
      the_model.flux_scale = 1e7
      the_model.good_fit = False
      the_model.p = None
      the_model.param_names = ['Undefined Model']
      

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

   def __str__(self):
      """Return a pretty print version of parameter values and errors."""
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
            up_sys_errs = N.nan_to_num((up_sys_errs**2 - ratios**2)**0.5)
            low_sys_errs = N.nan_to_num((low_sys_errs**2 - ratios**2)**0.5)
         except AttributeError:
            sys_flag = False         
         for i in xrange(len(self.param_names)):
            n=self.param_names[i][:m]
            t_n=n+(m-len(n))*' '
            frozen = '' if self.free[i] else '(FROZEN)'
            if not sys_flag:
               l+=[t_n+': (1 +/- %.3f) %.3g %s'%(ratios[i],p[i],frozen)]
            else:
               l+=[t_n+': (1 + %.3f - %.3f)(1 +/- %.3f) %.3g %s'%(up_sys_errs[i],low_sys_errs[i],ratios[i],p[i],frozen)]
         return '\n'.join(l)
      except:
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
      param_string=','.join( ( str(p) for p in self.p ) )
      a = eval(self.name+'(**self.__dict__)')
      a.p = [x for x in self.p]
      try: a.cov_matrix = self.cov_matrix.__copy__()
      except: pass
      return a

   def __flux_derivs__(self,*args):
      """Use finite differences to estimate the gradient of the integral flux wrt the spectral parameters."""
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

#===============================================================================================#

class PowerLaw(Model):
   """Implement a power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0         differential flux at e0 MeV
  gamma      (absolute value of) spectral index
      """
   def __call__(self,e):
      n0,gamma=self.p
      return (n0/self.flux_scale)*(self.e0/e)**gamma

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
      return (flux/self.flux_scale)*(1-gamma)*e**(-gamma)/(self.emax**(1-gamma)-self.emin**(1-gamma))

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
      return (n0/self.flux_scale)*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

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
      return (n0/self.flux_scale)*N.append( (e_break/e[e<e_break])**gamma1 , (e_break/e[e>=e_break])**gamma2 )

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
      return (n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-e/cutoff)

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
      return (n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-(e/cutoff)**b)

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

class Zero(Model):
   def __call__(self,e):
      return N.asarray(0*e)

#===============================================================================================#
if __name__=='__main__':
   disp=Dispersion()
   print disp(300,330)