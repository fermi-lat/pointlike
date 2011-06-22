"""A set of classes to implement spectral models.

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/Models.py,v 1.46 2011/06/18 01:51:03 lande Exp $

    author: Matthew Kerr, Joshua Lande

"""
import copy
import numpy as N
import numpy as np
import math as M
from scipy.integrate import quad
from scipy.interpolate import interp1d

#===============================================================================================#

class DefaultModelValues(object):
    """Static methods and class members to assign default values to the spectral models."""

    simple_models = {
        'PowerLaw'            : {'_p':[1e-11, 2.0],             'param_names':['Norm','Index'],'index_offset':0},
        'PowerLawFlux'        : {'_p':[1e-7 , 2.0],             'param_names':['Int_Flux','Index'],'emin':100,'emax':1e6},
        'BrokenPowerLaw'      : {'_p':[1e-11, 2.0, 2.0 ,1e3],   'param_names':['Norm','Index_1','Index_2', 'E_break']},
        'BrokenPowerLawFlux'  : {'_p':[1e-7, 2.0, 2.0 ,1e3],    'param_names':['Int_Flux','Index_1','Index_2', 'E_break'],'emin':100,'emax':1e6},
        'BrokenPowerLawCutoff': {'_p':[1e-11,2,2,1e3,3e3],      'param_names':['Norm','Index_1','Index_2','E_break','Cutoff']},
        'SmoothBrokenPowerLaw': {'_p':[1e-11,2.0,2.0, 1e3],     'param_names':['Norm','Index_1','Index_2','E_break'],'beta':0.1},
        'DoublePowerLaw'      : {'_p':[5e-12, 2.0, 2.0, 1],     'param_names':['Norm','Index_1','Index_2','Ratio']},
        'DoublePowerLawCutoff': {'_p':[5e-12,2,2,1e3,1],        'param_names':['Norm','Index_1','Index_2','Cutoff','Ratio']},
        'LogParabola'         : {'_p':[1e-11, 2.0, 1e-5,2e3],   'param_names':['Norm','Index','beta','E_break']},
        'ExpCutoff'           : {'_p':[1e-11, 2.0, 2e3],        'param_names':['Norm','Index','Cutoff']},
        'ExpCutoffPlusPL'     : {'_p':[1e-11,2.0,2e3,1e-12,1.5],'param_names':['Norm1','Index1','Cutoff1','Norm2','Index2']},
        'AllCutoff'           : {'_p':[1e-11, 1e3],             'param_names':['Norm','Cutoff']},
        'PLSuperExpCutoff'    : {'_p':[1e-11, 2.0, 2e3 ,1.],    'param_names':['Norm','Index','Cutoff', 'b']},
        'Constant'            : {'_p':[1.],                     'param_names':['Scale']},
        'InterpConstants'     : {'_p':[1.]*5,                   'param_names':['Scale_Vector'],'e_breaks':N.log10([100,300,1000,3000,3e5])},
        'FileFunction'        : {'_p':[1.],                     'param_names':['Norm']},
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
        the_model._p = None
        the_model.param_names = ['Undefined Model']
        the_model.background  = False

    @staticmethod
    def finish(the_model):
        """Common values that can be written once the model type has been sussed out."""
        npar = len(the_model._p)
        the_model.cov_matrix = N.zeros([npar,npar]) #default covariance matrix
        the_model.free = N.asarray([True] * npar)
        the_model._p = N.asarray(the_model._p) #redundant now
        

#===============================================================================================#

class Model(object):
    """Spectral model giving dN/dE for a point source.  Default units are ph/cm^2/s/MeV.
        Note that parameters are stored internally in logarithmic format.  This allows
        unconstrained minimization of the naturally-positive parameters."""
    def __init__(self,  **kwargs):
        """
        pars: list of parameter values: Model.XX( a,b) equivalent to Model.XX(p=[a,b])
Optional keyword arguments:

  =========    =======================================================
  Keyword      Description
  =========    =======================================================
  e0             [1000] value in MeV at which normalization is given
  flux_scale  [1e7] multiplier for actual value to make units more convenient
  p              [p1,p2,...] default values of spectral parameters; see docstring individual model classes
  simple_models [] for MixedModel, a list of simple model names composing MixedModel
  free          [True, True,...] a boolean list the same length as p giving the free (True) and fixed (False) parameters
  =========    =======================================================
      """
        iscopy = kwargs.pop('iscopy', False)
        DefaultModelValues.setup(self,**kwargs) # if called from copy method, will set p
        self.__dict__['_p'] = np.asarray(kwargs.pop('p', self._p),float)
        assert len(self._p)==len(self.param_names), 'Model: wrong number of parameters set: %s' % self._p
        self.__dict__.update(**kwargs)
        if not iscopy:
          assert np.all(self._p>0), 'fail parameter positivity constraint' 
          self._p = np.log10(self._p)
        self.free = np.asarray(self.free)

    def len(self):
        return len(self._p)
        
    def __getitem__(self, index):
        return self.getp(index)
        
    def __setitem__(self, index, value):
        self.setp(index,value)
        
    def getp(self, i, internal=False):
        """ get external value for parameter # i
        """
        i=self.mapper(i)
        return self._p[i] if internal else 10**(self._p[i])
    
    def get_all_parameters(self, internal=False):
        """ get a copy of the full set of parameters (external representation)"""
        return self._p.copy() if internal else 10**self._p
    
    def set_all_parameters(self, pars, internal=False):
        """ set all parameters (external representation)"""
        assert len(pars)== len(self._p)
        t = np.asarray(pars, float)
        self._p = t if internal else np.log10(t)

    def setp(self, i, par, internal=False):
        """ set internal value, convert unless inte
        """
        if not internal: 
            assert par>0, 'Model external parameter cannont be negative'
        i=self.mapper(i)
        self._p[i] = par if internal else  np.log10(par)
        
    def get_parameters(self):
        """Return FREE internal parameters ; used for spectral fitting."""
        return self._p[self.free]

    def set_parameters(self,new_vals):
        """Set FREE internal parameters; new_vals should have length equal to number of free parameters."""
        assert(len(new_vals)==(self.free).sum())
        self._p[self.free] = new_vals.astype('f') # downcast to float needed?

    def mapper(self,i):
        """ This object takes in a parameter and maps it to an index.
            Currently, this function is not particularly useful, but it
            could be generalized to allow lazier parameter selection. """
        if isinstance(i,str):
            if i not in self.param_names:
                raise Exception("Unknown parameter name %s" % i)
            return self.param_names.index(i)
        else:
            return i

    def freeze(self,i,freeze=True):
        """Freeze one of the spectral parameters from fitting.
        
            i: a parameter name or index.
            freeze    : if True, freeze parameter; if False, free it
            """
        i=self.mapper(i)
        self.free[i] = not freeze

    def thaw(self,i): self.freeze(i,freeze=False)

    def set_cov_matrix(self,new_cov_matrix):
        self.cov_matrix[N.outer(self.free,self.free)] = N.ravel(new_cov_matrix)

    def get_cov_matrix(self,absolute=True):
        """Return covariance matrix transformed out of log space."""
        p = 10**self._p if absolute else N.ones_like(self._p)
        jac = N.log10(N.exp(1)) #log_10(e)
        pt=p.reshape((p.shape[0],1)) #transpose
        return p*self.cov_matrix*pt/jac**2

    def get_free_errors(self):
        """Return the diagonal elements of the (log space) covariance matrix for free parameters."""
        return N.diag(self.cov_matrix)[self.free]**0.5

    def corr_coef(self):
        """Return the linear correlation coefficients for the estimated covariance matrix."""
        sigmas = N.diag(self.cov_matrix)**0.5
        return self.cov_matrix / N.outer(sigmas,sigmas)

    def statistical(self,absolute=False,two_sided=False):
        """Return the parameter values and fractional statistical errors.
            If no error estimates are present, return 0 for the fractional error.
        two_sided : bool
            if set, return 3 tuples: values, +errors, -errors    
            """
        p = 10**self._p #convert from log format
        z = N.zeros_like(p)
        vars = N.diag(self.cov_matrix)
        # this check for valid covariance matrix
        if N.all(self.cov_matrix==0) :
            return (p,z,z) if two_sided else (p,z) 
        try: #see if error estimates are present
            if not two_sided:
                vars = N.diag(self.get_cov_matrix(absolute=False))
                errs = vars**0.5 
                return p,errs*(p if absolute else N.ones_like(p))
            else:
                errs = vars**0.5
                lo_rat = (p-10**(self._p-errs))/(1. if absolute else p)
                hi_rat = (10**(self._p+errs)-p)/(1. if absolute else p)
                return p,hi_rat,lo_rat
        except:
            return (p,z,z) if two_sided else (p,z) 

    def __str__(self,absolute=False, indent=''):
        """Return a pretty print version of parameter values and errors.
          indent: string to prepend to each line (must be called explicitly)
        """
        #p,avg         = self.statistical(absolute=absolute,two_sided=False)
        p,hi_p,lo_p = self.statistical(absolute=absolute,two_sided=True)
        if hasattr(self,'index_offset'):
            p[1]=p[1]-self.index_offset #Index is parameter 1
        if not self.background:
            if not N.all(self.cov_matrix==0):
                f,fhi,flo    = self.i_flux(e_weight=0,two_sided=True,cgs=True,error=True)
                e,ehi,elo    = self.i_flux(e_weight=1,emax=3e5,two_sided=True,cgs=True,error=True)
                if not absolute:
                    fhi /= f; flo /= f; ehi /= e; elo /= e;
            else:
                f = self.i_flux(e_weight=0,cgs=True,error=False)
                e = self.i_flux(e_weight=1,cgs=True,error=False)
                fhi = flo = ehi = elo = 0
            p              = N.append(p, [f,e])
            hi_p          = N.append(hi_p, N.abs([fhi,ehi]))
            lo_p          = N.append(lo_p, N.abs([flo,elo]))
            pnames        = self.param_names + ['Ph. Flux','En. Flux']
        else: pnames = self.param_names

        l=[]
        if (not self.background and N.any(lo_p[0:-2]!=0)) or \
                 (self.background and N.any(lo_p!=0)): #if statistical errors are present    
            for i in xrange(len(pnames)):
                t_n = '%-10s' % pnames[i]
                if i < len(self._p):
                    # if free is empty (shouldn't happen normally) treat as all False
                    frozen = '' if len(self.free)>i and self.free[i] else '(FROZEN)'
                else:
                    frozen = '(DERIVED)'
                if not absolute:
                    low, high = max(0, lo_p[i]), max(0,hi_p[i]) 
                    if low>1e2 or high >1e2: 
                        low=high=0
                        frozen= '(Failed fit)'
                    q = [high,low,(high*low)**0.5]
                    if N.any(N.isnan(q)): q = 3*[0]
                    l+=[t_n+': (1 + %.3f - %.3f) (avg = %.3f) %-10.3g %s'% (tuple(q) + (p[i],frozen))]
                else:
                    l+=[t_n+': %.3g + %.3g - %.3g (avg = %.3g) %s'%(p[i],hi_p[i],lo_p[i],(hi_p[i]*lo_p[i])**0.5,frozen)]
            return ('\n'+indent).join(l)
        else: #if no errors are present
            for i in xrange(len(pnames)):
                t_n = '%-10s' % pnames[i]
                if i < len(self._p):
                    frozen = '' if self.free[i] else '(FROZEN)'
                else:
                    frozen = '(DERIVED)'
                l+=[t_n+': %.3g %s'%(p[i],frozen)]
            return ('\n'+indent).join(l)
            
    def __repr__(self): return self.__str__()

    def i_flux(self,emin=100,emax=N.inf,e_weight=0,cgs=False,error=False,two_sided=False):
        """Return the integral flux, \int_{emin}^{emax} dE E^{e_weight} dN/dE.
           e_weight = 0 gives the photon flux (ph cm^-2 s^-1)
           e_weight = 1 gives the energy flux (MeV cm^-2 s^-1) (see kwargs)

Optional keyword arguments:

  =========    =======================================================
  Keyword      Description
  =========    =======================================================
  emin         [100] lower bound in MeV
  emax         [N.inf] upper bound in MeV
  e_weight     [0] energy power by which to scale dN/dE
  cgs          [False] if True, energy in ergs, else in MeV
  error        [False] if True, return value is a tuple with flux and estimated error
  two_sided    [False] if True, return value is a triple with flux, estimated high and low errors
  =========    =======================================================
        """

        # check for a divergent flux
        if 100*self(100) <= 1e5*self(1e5): emax = min(5e5,emax)

        try:
            func    = self if e_weight == 0 else lambda e: self(e)*e**e_weight
            units  = 1.60218e-6**(e_weight) if cgs else 1. #extra factor from integral!
            epsabs = func(emin)*1e-4 # needed since epsrel does not seem to work
            flux    =  units*quad(func,emin,emax,epsabs=epsabs)[0]
            if not cgs: flux*=self.flux_scale #remove this?
            if error:
                args = (emin,emax,e_weight,cgs,False)
                d    = self.__flux_derivs__(*args)[self.free]
                dt  = d.reshape( (d.shape[0],1) ) #transpose
                err = (d * self.cov_matrix[self.free].transpose()[self.free] * dt).sum()**0.5
                if not two_sided:
                    return (flux,err)
                else: #use log transform to estimate two-sided errors
                    log_err  = err/flux
                    log_flux = N.log(flux)
                    return (flux,N.exp(log_flux+log_err)-flux,flux-N.exp(log_flux-log_err))

            return flux
        except Exception:
            print 'Encountered a numerical error when attempting to calculate integral flux.'
            return np.nan if not error else ([np.nan]*(3 if two_sided else 2))

    def set_flux(self,flux,*args,**kwargs):
        """ Set the flux of the source. 
                
            This function ensures that after the function, call,
                flux == model.i_flux(**kwargs)
            where args and kwargs is consistently passed into i_flux and set_flux
        """
        self._p[0] += N.log10(flux/self.i_flux(*args,**kwargs))

    def copy(self):
        
        a = eval(self.name+'(iscopy=True, **self.__dict__)') #create instance of same spectral model type
        
        a._p = N.asarray(self._p).copy() #copy in log values
        a.free = N.asarray(self.free).copy()
        try: a.cov_matrix = self.cov_matrix.__copy__()
        except: pass
        return a

    def fast_iflux(self,emin=100,emax=1e6):
        """Return a quick calculation for photon flux for models where it is analytically available."""
        return self.i_flux(emin=emin,emax=emax)


    def expected(self,emin,emax,exposure,skydir,event_class=-1,weighting_function=None):
        """Calculate the expected counts under a particular model.
            Include an optional weighting_function to calculate, e.g., the average PSF
            parameters under this spectral model."""
        
        from pointlike import DoubleVector
        lemin,lemax = N.log10([emin,emax])
        simpsn = max(16,(int(round((lemax-lemin)/0.1)) >> 1) << 1) #10 per decade
        points = N.logspace(lemin,lemax,simpsn+1)
        simpsf = points*N.log(emax/emin)*N.asarray([1.] + ([4.,2.]*(simpsn/2))[:-1] + [1.])/(3.*simpsn)
        
        if event_class < 0:
            exp     = N.asarray(exposure[0].vector_value(skydir,DoubleVector(points))) +\
                        N.asarray(exposure[1].vector_value(skydir,DoubleVector(points)))
        else:
            exp     = N.asarray(exposure[event_class].vector_value(skydir,DoubleVector(points)))
        
        expec = (self(points)*exp*simpsf).sum()
        if weighting_function is not None:
            return (weighting_function(points)*self(points)*exp*simpsf).sum()/expec
        return expec
        
    def __flux_derivs__(self,*args):
        """Use finite differences to estimate the gradient of the integral flux wrt the spectral parameters."""

        # note -- since spectral parameters are log transformed, just add/subtract a small amount in log space
        delta = 1e-5
        errs = N.asarray([delta] * len(self._p) )

        hi,lo = self.copy(),self.copy()
        derivs = []
        for i in xrange(len(self._p)):
            hi._p[i] += errs[i]
            lo._p[i] -= errs[i]
            derivs  += [(hi.i_flux(*args) - lo.i_flux(*args))/(2*errs[i])]
            hi._p[i] -= errs[i]
            lo._p[i] += errs[i]

        return N.asarray(derivs)
        
    def full_name(self):
        return self.pretty_name
    
    def set_e0(self, e0p):
        """ set a new reference energy, adjusting the norm parameter """
        # TODO: protect this
        gamma = 10** self._p[1]
        self._p[0] += gamma * N.log10(self.e0/e0p)
        self.e0 = e0p
        
    def pivot_energy(self):
        """ default to indicate no such """
        return None 
        

#===============================================================================================#

class PowerLaw(Model):
    """Implement a power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0            differential flux at e0 MeV
  gamma        (absolute value of) spectral index
        """
    def __call__(self,e):
        n0,gamma=10**self._p
        return (n0/self.flux_scale)*(self.e0/e)**(gamma-self.index_offset)

    def fast_iflux(self,emin=100,emax=1e6):
        n0,gamma = 10**self._p
        gamma -= self.index_offset
        return n0/(1-gamma)*self.e0**gamma*(emax**(1-gamma)-emin**(1-gamma))

    def gradient(self,e):
        n0,gamma = 10**self._p
        f = (n0/self.flux_scale)*(self.e0/e)**(gamma-self.index_offset)
        return N.asarray([f/n0,f*N.log(self.e0/e)])

    def pivot_energy(self):
        """ assuming a fit was done, estimate the pivot energy 
              (This only implemented for this model)
        """
        A  = 10**self._p[0]
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('PowerLaw fit required before calculating pivot energy')
        return self.e0*N.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ set a new reference energy, adjusting the norm parameter """
        # TODO: move this upstream
        gamma = 10** self._p[1]
        self._p[0] += gamma * N.log10(self.e0/e0p)
        self.e0 = e0p
        
    def full_name(self):
        return '%s, e0=%.0f'% (self.pretty_name,self.e0)

#===============================================================================================#

class PowerLawFlux(Model):
    """Implement a power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  flux         integrated flux from emin to emax MeV
  gamma        (absolute value of) spectral index
        """
    def __call__(self,e):
        flux,gamma=10**self._p
        return ((flux/self.flux_scale)*(1-gamma)/(self.emax**(1-gamma)-self.emin**(1-gamma)))*e**(-gamma)

    def fast_iflux(self,emin=100,emax=N.inf):
        n0,gamma = 10**self._p
        return n0*(emax**(1-gamma) - emin**(1-gamma)) / (self.emax**(1-gamma) - self.emin**(1-gamma))

    def gradient(self,e):
        flux,gamma = 10**self._p
        e1 = self.emax; e0 = self.emin
        d1 = e1**(1-gamma)
        d0 = e0**(1-gamma)
        t  = (M.log(e0)*d0 - M.log(e1)*d1)/(d1 - d0)
        f  = ((flux/self.flux_scale)*(1-gamma)/(d1-d0))*e**(-gamma)
        return N.asarray([f/flux,-f*(N.log(e) + (1./(1-gamma) + t))])

    def full_name(self):
        return '%s, emin=%.0f emax=%.0f'% (self.pretty_name,self.emin,self.emax)

#===============================================================================================#

class BrokenPowerLaw(Model):
    """Implement a broken power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0            differential flux at e0 MeV
  gamma1      (absolute value of) spectral index for e < e_break
  gamma2      (absolute value of) spectral index for e > e_break
  e_break     break energy
        """
    def __call__(self,e):
        n0,gamma1,gamma2,e_break=10**self._p
        return (n0/self.flux_scale)*(e_break/e)**np.where(e<e_break,gamma1,gamma2)

    def gradient(self,e):
        n0,gamma1,gamma2,e_break=10**self._p
        mask = e < e_break
        x = e_break/e
        lx = np.log(x)
        g = np.where(mask,gamma1,gamma2)
        f = (n0/self.flux_scale)*x**g
        return np.asarray([f/n0,f*lx*mask,f*lx*(~mask),f/e_break*g])
#===============================================================================================#

class BrokenPowerLawFlux(Model):
    """ Similar to PowerLawFlux for BrokenPowerLaw spectrum, the integral 
         flux is the free parameter rather than the Prefactor.

Spectral parameters:

     flux      integrated flux from emin to emax MeV
     gamma1  (absolute value of) spectral index for e < e_break
     gamma2  (absolute value of) spectral index for e > e_break
     e_break break energy
        """
    def __call__(self,e):
        flux,gamma1,gamma2,e_break=10**self._p
        if self.emax < e_break:
             norm=(1-gamma1)*e_break**(-gamma1)/(self.emax**(1-gamma1)-self.emin**(1-gamma1))
        elif self.emin > e_break:
             norm=(1-gamma2)*e_break**(-gamma2)/(self.emax**(1-gamma2)-self.emin**(1-gamma2))
        else:
             norm=1/(e_break**(gamma1)*(e_break**(1-gamma1)-self.emin**(1-gamma1))/(1-gamma1) + \
                        e_break**(gamma2)*(self.emax**(1-gamma2)-e_break**(1-gamma2))/(1-gamma1))

        return (flux/self.flux_scale)*norm*N.where( e < e_break, (e_break/e)**gamma1, (e_break/e)**gamma2 )

    def full_name(self):
        return '%s, emin=%.0f emax=%.0f'% (self.pretty_name,self.emin,self.emax)

#===============================================================================================#

class BrokenPowerLawCutoff(Model):
    """Implement a broken power law.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0            differential flux at e0 MeV
  gamma1      (absolute value of) spectral index for e < e_break
  gamma2      (absolute value of) spectral index for e > e_break
  e_break     break energy
        """
    def __call__(self,e):
        n0,gamma1,gamma2,e_break,cutoff=10**self._p
        return (n0/self.flux_scale)*N.where( e < e_break, (e_break/e)**gamma1, (e_break/e)**gamma2 )*N.exp(-e/cutoff)

#===============================================================================================#

class SmoothBrokenPowerLaw(Model):
    """Implement a smoothed broken power law. This is similar to a broken 
       powerlaw but uses the parameter beta to smoothly interpolate between 
       the two powerlaws.

       Everything is defined exactly the same as gtlike except for the 
       usual replacement that gamma1 and gamma2 are defined to be the 
       negative of the gtlike definition. 

       Note that this function has a somewhat confusing feature that when
       beta>0, the function always becomes softer for energies greater then
       the break. For beta<0, the function always becomes harder for energies
       greater than the break. So you need to pick beta in advance depending
       upon whether you want your function to slope up or down.

       For now, I am not allowing beta to be fit. It would be somewhat
       ackward trying to fit it since it may need negative values.

       Also note that e0 is used in the spectral function and can be set (but 
       not fit) by passing in e0 when initialized. """
    def __call__(self,e):
        n0,gamma1,gamma2,e_break=10**self._p
        return (n0/self.flux_scale)*(self.e0/e)**gamma1*\
               (1+(e_break/e)**((gamma1-gamma2)/self.beta))**(-self.beta)

    def gradient(self,e):
        """ lots of derivatives. You can derive them all by hitting the 
            derivative on log(dN/dE). """

        n0,gamma1,gamma2,e_break=10**self._p

        f=self(e)
        df_n0=f/n0

        bottom=1.0+(e_break/e)**(-(gamma1-gamma2)/self.beta)
        df_gamma1=f*(N.log(self.e0/e)-N.log(e_break/e)/bottom)
        df_gamma2=f*N.log(e_break/e)/bottom
        df_break=-f*(gamma1-gamma2)/e_break/bottom

        # for completeness, here is the gradient with respect to beta.
        # For now, no compelling reason to fit it.
        # df_beta=f*(-N.log(1.0+(e_break/e)**((gamma1-gamma2)/self.beta))
        #            +((gamma1-gamma2)/self.beta)*N.log(e_break/e)/bottom)

        return N.asarray([df_n0,df_gamma1,df_gamma2,df_break])

    def full_name(self):
        return '%s, e0=%.0f, beta=%.3g'% (self.pretty_name,self.e0,self.beta)

#===============================================================================================#

class DoublePowerLaw(Model):
    """Spectral model is the sum of two indepedendent power laws.  E.g., the Crab Nebula = IC + synch.

Spectral parameters:

  n0             differential flux at e0 MeV for first power law
  gamma1        (absolute value of) spectral index for first power law
  gamma2        (absolute value of) spectral index for second power law
  ratio         ratio of the differential fluxes of first and second power law at e0
        """
    def __call__(self,e):
        n0,gamma1,gamma2,ratio=10**self._p
        return (n0/self.flux_scale)*((self.e0/e)**gamma1 + ratio*(self.e0/e)**gamma2)

#===============================================================================================#

class DoublePowerLawCutoff(Model):
    """Spectral model is the sum of two indepedendent power laws, one with a cutoff.  E.g., a pulsar + PWN.

Spectral parameters:

  n0             differential flux at e0 MeV for first power law
  gamma1        (absolute value of) spectral index for first power law
  gamma2        (absolute value of) spectral index for second power law
  cutoff        cutoff -- note goes with gamma!
  ratio         ratio of the differential fluxes of first and second power law at e0
        """
    def __call__(self,e):
        n0,gamma1,gamma2,cutoff,ratio=10**self._p
        return (n0/self.flux_scale)*((self.e0/e)**gamma1*N.exp(-e/cutoff) + ratio*(self.e0/e)**gamma2)


#===============================================================================================#


class LogParabola(Model):
    """Implement a log parabola (for blazars.)  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0            differential flux at e0 MeV
  alpha        (absolute value of) constant spectral index
  beta         co-efficient for energy-dependent spectral index
  e_break     break energy
        """
    def __call__(self,e):
        n0,alpha,beta,e_break=10**self._p
        #alpha -= self.index_offset
#        return (n0/self.flux_scale)*(e_break/e)**(alpha - beta*N.log(e_break/e))
        x = N.log(e_break/e)
        y = (alpha - beta*x)*x
        return (n0/self.flux_scale) * N.exp(y) #N.clip(y, -10, 100)) #protect over, underflows

    def gradient(self,e):
        n0,alpha,beta,e_break=10**self._p
        #alpha -= self.index_offset
        #f  = (n0/self.flux_scale)*(e_break/e)**(alpha - beta*N.log(e_break/e))
        #log_term = N.log(e_break/e)
        #return N.asarray([f/n0,f*log_term,-f*log_term**2,f*alpha/e_break])
        x =N.log(e_break/e)
        y = (alpha - beta*x)*x
        f = (n0/self.flux_scale)*N.exp(y) # N.clip(y, -10, 100))
        return N.asarray([f/n0, f*x, -f*x**2, f*alpha/e_break])

    def pivot_energy(self):
        """  estimate the pivot energy, assuming that the beta part is not so important
        """
        #if sum(self.free) != 2 or self._p[2]>-2:
        #    raise Exception( 'Cannot calculate pivot energy for LogParabola unless in power-law mode')
        A  = 10**self._p[0]
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('Models.LogParabola: Fit required before calculating pivot energy')
        ebreak = 10**self._p[3]
        return ebreak*N.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ set a new break energy, adjusting the norm parameter """
        ebreak = 10** self._p[3]
        gamma = 10** self._p[1]
        self._p[0] += gamma * N.log10(ebreak/e0p)
        self._p[3] = N.log10(e0p)
        self.e0 = e0p
 
    def create_powerlaw(self, beta_max=3e-2):
        """ if beta is small and fixed, return an equivalent PowerLaw, otherwise just return self """
        if self[2]>beta_max or self.free[2]: return self
        nm = PowerLaw(p=self[0:2], e0=self[3])
        nm.cov_matrix=self.cov_matrix[:-2,:-2]
        return nm
    
#===============================================================================================#

class ExpCutoff(Model):
    """Implement a power law with exponential cutoff.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0            differential flux at e0 MeV
  gamma        (absolute value of) spectral index
  cutoff      e-folding cutoff energy (MeV)
        """
    def __call__(self,e):
        n0,gamma,cutoff=10**self._p
        return (n0/self.flux_scale) * (self.e0/e)**gamma * N.exp(-e/cutoff)

    def gradient(self,e):
        n0,gamma,cutoff = 10**self._p
        f = (n0/self.flux_scale) * (self.e0/e)**gamma * N.exp(-e/cutoff)
        return N.asarray([f/n0,f*N.log(self.e0/e),f*e/cutoff**2])

    def pivot_energy(self):
        """ assuming a fit was done, estimate the pivot energy 
              
        """
        A  = 10**self._p[0]
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('%s fit required before calculating pivot energy' %self.name)
        return self.e0*N.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ set a new reference energy, adjusting the norm parameter """
        gamma = 10** self._p[1]
        self._p[0] += gamma * N.log10(self.e0/e0p)
        self.e0 = float(e0p) 
        
    def full_name(self):
        return '%s, e0=%.0f'% (self.pretty_name,self.e0)
        
#===============================================================================================#

class ExpCutoffPlusPL(Model):
    """Implement a power law with exponential cutoff + an additional power law.  A la pulsar + PWN.

Spectral parameters:

  n0_1         differential flux at e0 MeV
  gamma_1     (absolute value of) spectral index
  cutoff_1    e-folding cutoff energy (MeV)
  n0_2
  gamma_2
        """
    def __call__(self,e):
        n0_1,gamma_1,cutoff_1,n0_2,gamma_2 = 10**self._p
        return (n0_1/self.flux_scale)*(self.e0/e)**gamma_1*N.exp(-e/cutoff_1) + (n0_2/self.flux_scale)*(self.e0/e)**gamma_2

#===============================================================================================#

class AllCutoff(Model):
    """Implement an exponential cutoff.  This for the case when cutoff too low to constrain index.
        See constructor docstring for further keyword arguments.

Spectral parameters:

  n0            differential flux at e0 MeV
  cutoff      e-folding cutoff energy (MeV)
        """
    def __call__(self,e):
        n0,cutoff=10**self._p
        if cutoff < 0: return 0
        return (n0/self.flux_scale)*N.exp(-e/cutoff)

#===============================================================================================#

class PLSuperExpCutoff(Model):
    """Implement a power law with hyper-exponential cutoff.  See constructor docstring for further keyword arguments.

Spectral parameters:

  n0            differential flux at e0 MeV
  gamma        (absolute value of) spectral index
  cutoff      e-folding cutoff energy (MeV)
  b             additional power in the exponent
        """
    def __call__(self,e):
        n0,gamma,cutoff,b=10**self._p
        return (n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-(e/cutoff)**b)

    def gradient(self,e):
        n0,gamma,cutoff,b = 10**self._p
        f = (n0/self.flux_scale)*(self.e0/e)**gamma*N.exp(-(e/cutoff)**b)
        return N.asarray([f/n0,f*N.log(self.e0/e),
                     f*(b/cutoff)*(e/cutoff)**b,f*(e/cutoff)**b*N.log(cutoff/e)])


#===============================================================================================#

class MixedModel(Model):
    """Implement a composite model.  The value is the sum of the simple models.
        See constructor docstring for further keyword arguments.
        NOTA BENE: specify the simple models via the keyword arguments 'models'
        """
    def __call__(self,e):
        counter = 0
        for i in xrange(len(self.n)):
            self.spec_models[i]._p = self._p[counter:counter+self.n[i]]
            counter += self.n[i]
        return N.array([model(e) for model in self.spec_models]).sum(axis=0)


#===============================================================================================#

class Constant(Model):
    def __call__(self,e):
        return N.ones_like(e)*10**self._p[0]
    
    def fast_iflux(self,emin=100,emax=1e6):
        return (emax-emin)*10**self._p[0]

    def gradient(self,e):
        return N.ones_like(e)

#===============================================================================================#

class InterpConstants(Model):

    def __call__(self,e):
        interp = interp1d(self.e_breaks,10**self._p)
        return interp(N.log10(e))

    def set_flux(self,flux,**kwargs):
        raise NotImplementedError("No way to set flux for InterpConstants spectral model")

#===============================================================================================#

class FileFunction(Model):

    def __make_interp__(self):
        self.interp = interp1d(N.log10(self.energy),N.log10(self.flux),
                bounds_error=False,fill_value=-N.inf)

    def __init__(self,*args,**kwargs):

        super(FileFunction,self).__init__(*args,**kwargs)

        if not hasattr(self,'file'):
            raise Exception("FileFunction must be created with a file.")

        file=N.genfromtxt(self.file,unpack=True)
        self.energy,self.flux=file[0],file[1]

        self.__make_interp__()

    def __call__(self,e):
        return 10**(self._p[0]+self.interp(N.log10(e)))
    
    def gradient(self,e):
        return N.asarray([10**self.interp(N.log10(e))])

    def __getstate__(self):
        """ You cannot pickle an interp1d object. """
        d=copy.copy(self.__dict__)
        del d['interp']
        return d

    def __setstate__(self,state):
        """ recreate interp1d object. """
        self.__dict__ = state
        self.__make_interp__()


def convert_exp_cutoff(model):
    # this function need for XML parsing
    if model.name != 'ExpCutoff':
        raise Exception,'Cannot process %s into PLSuperExpCutoff'%(model.name)
    nm = PLSuperExpCutoff()
    nm._p    = N.append(model._p,0)
    nm.free = N.append(model.free,False)
    nm.cov_matrix[:,:] = 0
    nm.cov_matrix[:-1,:-1] = model.cov_matrix[:,:]
    nm.e0 = model.e0
    return nm
    
    
    
