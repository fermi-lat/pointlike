"""A set of classes to implement spectral models.

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/Models.py,v 1.81 2012/02/06 21:57:08 wallacee Exp $

    author: Matthew Kerr, Joshua Lande
"""
import os
import copy
import collections
import operator
import numpy as N
import numpy as np
from abc import abstractmethod

from scipy.integrate import quad
from scipy.interpolate import interp1d

#===============================================================================================#

class DefaultModelValues(object):
    """ Static methods and class members to assign default values to the spectral models. """

    simple_models = dict(
        PowerLaw            = dict(_p=[1e-11, 2.0],             param_names=['Norm','Index'],index_offset=0, e0=1e3),
        PowerLawFlux        = dict(_p=[1e-7 , 2.0],             param_names=['Int_Flux','Index'],emin=100,emax=1e6),
        BrokenPowerLaw      = dict(_p=[1e-11, 2.0, 2.0 ,1e3],   param_names=['Norm','Index_1','Index_2', 'E_break']),
        BrokenPowerLawFlux  = dict(_p=[1e-7, 2.0, 2.0 ,1e3],    param_names=['Int_Flux','Index_1','Index_2', 'E_break'],emin=100,emax=1e6),
        BrokenPowerLawCutoff= dict(_p=[1e-11,2,2,1e3,3e3],      param_names=['Norm','Index_1','Index_2','E_break','Cutoff']),
        SmoothBrokenPowerLaw= dict(_p=[1e-11,2.0,2.0, 1e3],     param_names=['Norm','Index_1','Index_2','E_break'],beta=0.1, e0=1e3),
        DoublePowerLaw      = dict(_p=[5e-12, 2.0, 2.0, 1],     param_names=['Norm','Index_1','Index_2','Ratio'], e0=1e3),
        DoublePowerLawCutoff= dict(_p=[5e-12,2,2,1e3,1],        param_names=['Norm','Index_1','Index_2','Cutoff','Ratio'], e0=1e3),
        LogParabola         = dict(_p=[1e-11, 2.0, 1e-5,2e3],   param_names=['Norm','Index','beta','E_break']),
        ExpCutoff           = dict(_p=[1e-11, 2.0, 2e3],        param_names=['Norm','Index','Cutoff'], e0=1e3),
        ExpCutoffPlusPL     = dict(_p=[1e-11,2.0,2e3,1e-12,1.5],param_names=['Norm1','Index1','Cutoff1','Norm2','Index2'], e0=1e3),
        AllCutoff           = dict(_p=[1e-11, 1e3],             param_names=['Norm','Cutoff']),
        PLSuperExpCutoff    = dict(_p=[1e-11, 2.0, 2e3 ,1.],    param_names=['Norm','Index','Cutoff', 'b'], e0=1e3),
        Constant            = dict(_p=[1.],                     param_names=['Scale']),
        InterpConstants     = dict(_p=[1.]*5,                   param_names=['Scale_Vector'],e_breaks=np.log10([100,300,1000,3000,3e5])),
        FileFunction        = dict(_p=[],                       param_names=[]),
        DMFitFunction       = dict(_p=[1.0, 100.],              param_names=['sigmav','mass'], norm=1, bratio=1.0, channel0=1, channel1=1,
                                   file='$(INST_DIR)/Likelihood/src/dmfit/gammamc_dif.dat'),
        )

    @staticmethod
    def setup(the_model):
        
        DefaultModelValues.start(the_model)
        classname = the_model.name = the_model.pretty_name = the_model.__class__.__name__
        
        if classname in DefaultModelValues.simple_models:
            for key,val in DefaultModelValues.simple_models[classname].items():
                the_model.__dict__[key]=val

        DefaultModelValues.finish(the_model)

    @staticmethod
    def start(the_model):
        """ Common values independent of the model type."""
        the_model.background  = False

    @staticmethod
    def finish(the_model):
        """Common values that can be written once the model type has been sussed out."""
        npar = len(the_model._p)
        the_model.cov_matrix = np.zeros([npar,npar]) #default covariance matrix
        the_model.free = np.asarray([True] * npar)
        the_model._p = np.asarray(the_model._p) #redundant now
        

#===============================================================================================#

class Model(object):
    """ Spectral model giving dN/dE for a point source.  Default units are ph/cm^2/s/MeV.
        Note that parameters are stored internally in logarithmic format.  This allows
        unconstrained minimization of the naturally-positive parameters."""
    def __init__(self,  **kwargs):
        """ pars: list of parameter values: Model.XX( a,b) equivalent to Model.XX(p=[a,b])

            Optional keyword arguments:

                =========    =======================================================
                Keyword      Description
                =========    =======================================================
                p              [p1,p2,...] default values of spectral parameters; see docstring individual model classes
                free           [True, True,...] a boolean list the same length as p giving the free (True) and fixed (False) parameters
                =========    =======================================================
        """
        DefaultModelValues.setup(self) # if called from copy method, will set p
        self._p = np.asarray(kwargs.pop('p', self._p),float)
        assert len(self._p)==len(self.param_names), 'Model: wrong number of parameters set: %s' % self._p

        # This has to come before the kwargs parsing
        # incase the __getitem__ function is called
        # with an input spectral parameter.
        assert np.all(self._p>0), 'fail parameter positivity constraint' 
        self._p = np.log10(self._p)

        # deal with other kwargs. If they are in the list of spectral
        # parameters, update the corresponding value. Tag them to 
        # the class.
        for k,v in kwargs.items():
            if k in self: 
                self[k]=v
            else:
                self.__dict__[k]=v

        self.free = np.asarray(self.free)

    def len(self):
        return len(self._p)

    def __contains__(self, item):
        try:
            self.__getitem__(item)
            return True
        except Exception, ex:
            return False
        
    def __getitem__(self, index):
        return self.getp(index)
        
    def __setitem__(self, index, value):
        self.setp(index,value)
        
    def getp(self, i, internal=False):
        """ Get external value for parameter # i """
        i=self.mapper(i)
        return self._p[i] if internal else 10**(self._p[i])
    
    def get_all_parameters(self, internal=False):
        """ Get a copy of the full set of parameters (external representation)"""
        return self._p.copy() if internal else 10**self._p
    
    def set_all_parameters(self, pars, internal=False):
        """ set all parameters (external representation)"""
        assert len(pars)== len(self._p)
        t = np.asarray(pars, float)
        self._p = t if internal else np.log10(t)

    def setp(self, i, par, internal=False):
        """ set internal value, convert unless internal. """
        i=self.mapper(i)
        if not internal: 
            assert par>0, 'Model external parameter cannont be negative'
        self._p[i] = par if internal else  np.log10(par)
        
    def get_parameters(self):
        """ Return FREE internal parameters ; used for spectral fitting."""
        return self._p[self.free]

    def set_parameters(self,new_vals):
        """ Set FREE internal parameters; new_vals should have length equal to number of free parameters."""
        assert(len(new_vals)==(self.free).sum())
        self._p[self.free] = new_vals.astype('f') # downcast to float needed?

    def mapper(self,i):
        """ This object takes in a parameter and maps it to an index.
            Currently, this function is not particularly useful, but it
            could be generalized to allow lazier parameter selection. """
        if isinstance(i,str):
            if i.lower() not in np.char.lower(self.param_names):
                raise Exception("Unknown parameter name %s" % i)
            return np.where(np.char.lower(self.param_names)==i.lower())[0][0]
        else:
            return i

    def freeze(self,i,freeze=True):
        """ Freeze one of the spectral parameters from fitting.
        
                i: a parameter name or index.
                freeze    : if True, freeze parameter; if False, free it
            """
        i=self.mapper(i)
        self.free[i] = not freeze

    def thaw(self,i): self.freeze(i,freeze=False)

    def set_cov_matrix(self,new_cov_matrix):
        self.cov_matrix[np.outer(self.free,self.free)] = np.ravel(new_cov_matrix)

    def get_cov_matrix(self,absolute=True):
        """Return covariance matrix transformed out of log space."""
        p = 10**self._p if absolute else np.ones_like(self._p)
        jac = np.log10(np.exp(1)) #log_10(e)
        pt=p.reshape((p.shape[0],1)) #transpose
        return p*self.cov_matrix*pt/jac**2

    def get_free_errors(self):
        """Return the diagonal elements of the (log space) covariance matrix for free parameters."""
        return np.diag(self.cov_matrix)[self.free]**0.5

    def corr_coef(self):
        """Return the linear correlation coefficients for the estimated covariance matrix."""
        sigmas = np.diag(self.cov_matrix)**0.5
        return self.cov_matrix / np.outer(sigmas,sigmas)

    def statistical(self,absolute=False,two_sided=False):
        """ Return the parameter values and fractional statistical errors.
            If no error estimates are present, return 0 for the fractional error.

            two_sided : bool
                if set, return 3 tuples: values, +errors, -errors    
        """
        p = 10**self._p #convert from log format
        z = np.zeros_like(p)
        vars = np.diag(self.cov_matrix)
        vars[vars<0]=0
        # this check for valid covariance matrix
        if np.all(self.cov_matrix==0) :
            return (p,z,z) if two_sided else (p,z) 
        try: #see if error estimates are present
            if not two_sided:
                vars = np.diag(self.get_cov_matrix(absolute=False))
                vars[vars<0]=0
                errs = vars**0.5 
                return p,errs*(p if absolute else np.ones_like(p))
            else:
                errs = vars**0.5
                lo_rat = (p-10**(self._p-errs))/(1. if absolute else p)
                hi_rat = (10**(self._p+errs)-p)/(1. if absolute else p)
                return p,hi_rat,lo_rat
        except:
            return (p,z,z) if two_sided else (p,z) 

    def error(self,i, internal=False):
        """ get error for parameter # i """
        i=self.mapper(i)
        return (np.diag(self.get_cov_matrix(absolute=not internal))**0.5)[i]

    def __str__(self,absolute=False, indent=''):
        """ Return a pretty print version of parameter values and errors.
            indent: string to prepend to each line (must be called explicitly)
        """
        #p,avg         = self.statistical(absolute=absolute,two_sided=False)
        p,hi_p,lo_p = self.statistical(absolute=absolute,two_sided=True)
        if hasattr(self,'index_offset'):
            p[1]=p[1]-self.index_offset #Index is parameter 1
            hi_p[1]-=self.index_offset
            lo_p[1]-=self.index_offset
        if not self.background:
            if not np.all(self.cov_matrix==0):
                f,fhi,flo    = self.i_flux(e_weight=0,emax=3e5,two_sided=True,cgs=True,error=True)
                e,ehi,elo    = self.i_flux(e_weight=1,emax=3e5,two_sided=True,cgs=True,error=True)
                if not absolute:
                    fhi /= f; flo /= f; ehi /= e; elo /= e;
            else:
                f = self.i_flux(e_weight=0,cgs=True,error=False)
                e = self.i_flux(e_weight=1,cgs=True,error=False)
                fhi = flo = ehi = elo = 0
            p              = np.append(p, [f,e])
            hi_p          = np.append(hi_p, np.abs([fhi,ehi]))
            lo_p          = np.append(lo_p, np.abs([flo,elo]))
            pnames        = self.param_names + ['Ph. Flux','En. Flux']
        else: pnames = self.param_names

        l=[]
        if (not self.background and np.any(lo_p[0:-2]!=0)) or \
                 (self.background and np.any(lo_p!=0)): #if statistical errors are present    
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
                    if np.any(np.isnan(q)): q = 3*[0]
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

    def i_flux(self,emin=100,emax=np.inf,e_weight=0,cgs=False,error=False,two_sided=False, quiet=False):
        """ Return the integral flux, \int_{emin}^{emax} dE E^{e_weight} dN/dE.
            e_weight = 0 gives the photon flux (ph cm^-2 s^-1)
            e_weight = 1 gives the energy flux (MeV cm^-2 s^-1) (see kwargs)

            Optional keyword arguments:

                =========    =======================================================
                Keyword      Description
                =========    =======================================================
                emin         [100] lower bound in MeV
                emax         [np.inf] upper bound in MeV
                e_weight     [0] energy power by which to scale dN/dE
                cgs          [False] if True, energy in ergs, else in MeV
                error        [False] if True, return value is a tuple with flux and estimated error
                two_sided    [False] if True, return value is a triple with flux, estimated high and low errors
                quiet        [False] set to suppress error message
                =========    =======================================================
        """

        # check for a divergent flux
        if 100*self(100) <= 1e5*self(1e5): emax = min(5e5,emax)

        try:
            func    = self if e_weight == 0 else lambda e: self(e)*e**e_weight
            units  = 1.60218e-6**(e_weight) if cgs else 1. #extra factor from integral!
            epsabs = func(emin)*1e-4 # needed since epsrel does not seem to work
            flux    =  units*quad(func,emin,emax,epsabs=epsabs,full_output=True)[0]
            if error:
                # will silently ignore 'free' parameters without errors
                mask = (self.free) * (self.cov_matrix.diagonal()>0)
                args = (emin,emax,e_weight,cgs,False)
                d    = self.__flux_derivs__(*args)[mask]
                dt  = d.reshape( (d.shape[0],1) ) #transpose
                err = (d * self.cov_matrix[mask].transpose()[mask] * dt).sum()**0.5
                if not two_sided:
                    return (flux,err)
                else: #use log transform to estimate two-sided errors
                    log_err  = err/flux
                    log_flux = np.log(flux)
                    return (flux,np.exp(log_flux+log_err)-flux,flux-np.exp(log_flux-log_err))

            return flux
        except Exception, msg:
            if not quiet:
                print 'Encountered a numerical error, "%s", when attempting to calculate integral flux.'%msg
            return np.nan if not error else ([flux, np.nan,np.nan] if two_sided else [flux, np.nan])

    def set_flux(self,flux,*args,**kwargs):
        """ Set the flux of the source. 
                
            This function ensures that after the function call,
                flux == model.i_flux(**kwargs)
            where args and kwargs is consistently passed into i_flux and set_flux

            For example:

                >>> model = PowerLaw(index=2)
                >>> model.set_flux(1e-7, emin=1e3, emax=1e5)
                >>> print '%g' % model.i_flux(emin=1e3, emax=1e5)
                1e-07

            Note that this implementation is robust even when the source
            has an initial flux of 0:

                >>> model.setp(0, -np.inf, internal=True)
                >>> model.set_flux(1e-7)
                >>> print '%g' % model.i_flux()
                1e-07

            Just a simple example of a unitless powerlaw:

                >>> model = PowerLaw(norm=1, index=2, e0=1)
                >>> energies = np.logspace(0,6,100)
                >>> np.allclose(model(energies), energies**-2)
                True
        """
        self.setp(0, 0, internal=True)
        new_prefactor = np.log10(flux/self.i_flux(*args,**kwargs))
        self.setp(0,new_prefactor,internal=True)

    def copy(self): return copy.deepcopy(self)

    def fast_iflux(self,emin=100,emax=1e6):
        """Return a quick calculation for photon flux for models where it is analytically available."""
        return self.i_flux(emin=emin,emax=emax)


    def expected(self,emin,emax,exposure,skydir,event_class=-1,weighting_function=None):
        """ Calculate the expected counts under a particular model.
            Include an optional weighting_function to calculate, e.g., the average PSF
            parameters under this spectral model."""
        
        from pointlike import DoubleVector
        lemin,lemax = np.log10([emin,emax])
        simpsn = max(16,(int(round((lemax-lemin)/0.1)) >> 1) << 1) #10 per decade
        points = np.logspace(lemin,lemax,simpsn+1)
        simpsf = points*np.log(emax/emin)*np.asarray([1.] + ([4.,2.]*(simpsn/2))[:-1] + [1.])/(3.*simpsn)
        
        if event_class < 0:
            exp     = np.asarray(exposure[0].vector_value(skydir,DoubleVector(points))) +\
                        np.asarray(exposure[1].vector_value(skydir,DoubleVector(points)))
        else:
            exp     = np.asarray(exposure[event_class].vector_value(skydir,DoubleVector(points)))
        
        expec = (self(points)*exp*simpsf).sum()
        if weighting_function is not None:
            return (weighting_function(points)*self(points)*exp*simpsf).sum()/expec
        return expec
        
    def __flux_derivs__(self,*args):
        """ Use finite differences to estimate the gradient of the integral flux wrt the spectral parameters."""

        # note -- since spectral parameters are log transformed, just add/subtract a small amount in log space
        delta = 1e-5
        errs = np.asarray([delta] * len(self._p) )
        hi,lo = self.copy(),self.copy()
        derivs = []
        for i in xrange(len(self._p)):
            hi.setp(i,hi._p[i] + errs[i],internal=True)
            lo.setp(i,lo._p[i] - errs[i],internal=True)
            derivs  += [(hi.i_flux(*args) - lo.i_flux(*args))/(2*errs[i])]

        return np.asarray(derivs)
        
    def full_name(self):
        return self.pretty_name
        
    def pivot_energy(self):
        """ default to indicate no such """
        return None 

    def zero(self):
        if self.getp(0,internal=True)==-np.inf: return #already zeroed
        self.old_flux = self.getp(0, internal=True) # get the internal value
        self.setp(0,  -np.inf, internal=True)
        self.old_free = self.free.copy()
        self.free[:] = False

    def unzero(self):
        assert self.old_flux!=-np.inf, 'attempt to unzero non-zeroed source %d ' % which
        self.setp(0, self.old_flux, internal=True)
        self.free = self.old_free.copy()

#===============================================================================================#

class PowerLaw(Model):
    """ Implement a power law.  See constructor docstring for further keyword arguments.

        Spectral parameters:

            n0            differential flux at e0 MeV
            gamma        (absolute value of) spectral index

        It is easy to create a PowerLaw and access and set it's values:

            >>> model=PowerLaw(p=[1e-7,2])
            >>> model['Norm'] == 1e-7 and model['norm'] == 1e-7
            True

        It is easy to set values internally
            
            >>> model['Norm'] = 1e-6
            >>> print model['Norm']
            1e-06

        When you create the object, you can specify the parameters by kwargs:

           >>> model=PowerLaw(norm=1e-11, index=1)
           >>> model['index'] == 1 and model['norm'] == 1e-11
           True

        You can also make deep copies of the spectral models:

            >>> copy_model=model.copy()
            >>> print copy_model['norm'] 
            1e-11
            >>> model['norm'] = 1e-10
            >>> print copy_model['norm']
            1e-11
            >>> print model['norm']
            1e-10
    """
    def __call__(self,e):
        n0,gamma=10**self._p
        return n0*(self.e0/e)**(gamma-self.index_offset)

    def fast_iflux(self,emin=100,emax=1e6):
        n0,gamma = 10**self._p
        gamma -= self.index_offset
        return n0/(1-gamma)*self.e0**gamma*(emax**(1-gamma)-emin**(1-gamma))

    def gradient(self,e):
        n0,gamma = 10**self._p
        f = n0*(self.e0/e)**(gamma-self.index_offset)
        return np.asarray([f/n0,f*np.log(self.e0/e)])

    def pivot_energy(self):
        """ assuming a fit was done, estimate the pivot energy 
              (This only implemented for this model)
        """
        A  = 10**self._p[0]
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('PowerLaw fit required before calculating pivot energy')
        return self.e0*np.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ set a new reference energy, adjusting the norm parameter """
        # TODO: move this upstream
        gamma = 10** self._p[1]
        self._p[0] += gamma * np.log10(self.e0/e0p)
        self.e0 = e0p
        
    def full_name(self):
        return '%s, e0=%.0f'% (self.pretty_name,self.e0)
    @property
    def eflux(self):
        return self.e0**2*10**self._p[0]
#===============================================================================================#

class PowerLawFlux(Model):
    """ Implement a power law.  See constructor docstring for further keyword arguments.

        Spectral parameters:

            flux         integrated flux from emin to emax MeV
            gamma        (absolute value of) spectral index

        PowerLawFlux allows the flux to be set:

            >>> model = PowerLawFlux(Int_Flux=15, emin=50, emax=1000)
            >>> np.allclose([model.i_flux(50, 1000),model.fast_iflux(50, 1000)], [15,15])
            True
    """
    def __call__(self,e):
        flux,gamma=10**self._p
        return (flux*(1-gamma)/(self.emax**(1-gamma)-self.emin**(1-gamma)))*e**(-gamma)

    def fast_iflux(self,emin=100,emax=np.inf):
        n0,gamma = 10**self._p
        return n0*(emax**(1-gamma) - emin**(1-gamma)) / (self.emax**(1-gamma) - self.emin**(1-gamma))

    def gradient(self,e):
        flux,gamma = 10**self._p
        e1 = self.emax; e0 = self.emin
        d1 = e1**(1-gamma)
        d0 = e0**(1-gamma)
        t  = (np.log(e0)*d0 - np.log(e1)*d1)/(d1 - d0)
        f  = (flux*(1-gamma)/(d1-d0))*e**(-gamma)
        return np.asarray([f/flux,-f*(np.log(e) + (1./(1-gamma) + t))])

    def full_name(self):
        return '%s, emin=%.0f emax=%.0f'% (self.pretty_name,self.emin,self.emax)

#===============================================================================================#

class BrokenPowerLaw(Model):
    """ Implement a broken power law.  See constructor docstring for further keyword arguments.

        Spectral parameters:

            n0            differential flux at e0 MeV
            gamma1      (absolute value of) spectral index for e < e_break
            gamma2      (absolute value of) spectral index for e > e_break
            e_break     break energy
    """
    def __call__(self,e):
        n0,gamma1,gamma2,e_break=10**self._p
        return n0*(e_break/e)**np.where(e<e_break,gamma1,gamma2)

    def gradient(self,e):
        n0,gamma1,gamma2,e_break=10**self._p
        mask = e < e_break
        x = e_break/e
        lx = np.log(x)
        g = np.where(mask,gamma1,gamma2)
        f = n0*x**g
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

        return flux*norm*np.where( e < e_break, (e_break/e)**gamma1, (e_break/e)**gamma2 )

    def full_name(self):
        return '%s, emin=%.0f emax=%.0f'% (self.pretty_name,self.emin,self.emax)

#===============================================================================================#

class BrokenPowerLawCutoff(Model):
    """ Implement a broken power law.  See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            gamma1      (absolute value of) spectral index for e < e_break
            gamma2      (absolute value of) spectral index for e > e_break
            e_break     break energy
    """
    def __call__(self,e):
        n0,gamma1,gamma2,e_break,cutoff=10**self._p
        return n0*np.where( e < e_break, (e_break/e)**gamma1, (e_break/e)**gamma2 )*np.exp(-e/cutoff)

#===============================================================================================#

class SmoothBrokenPowerLaw(Model):
    """ Implement a smoothed broken power law. This is similar to a broken 
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
        return n0*(self.e0/e)**gamma1*\
               (1+(e_break/e)**((gamma1-gamma2)/self.beta))**(-self.beta)

    def gradient(self,e):
        """ lots of derivatives. You can derive them all by hitting the 
            derivative on log(dN/dE). """

        n0,gamma1,gamma2,e_break=10**self._p

        f=self(e)
        df_n0=f/n0

        bottom=1.0+(e_break/e)**(-(gamma1-gamma2)/self.beta)
        df_gamma1=f*(np.log(self.e0/e)-np.log(e_break/e)/bottom)
        df_gamma2=f*np.log(e_break/e)/bottom
        df_break=-f*(gamma1-gamma2)/e_break/bottom

        # for completeness, here is the gradient with respect to beta.
        # For now, no compelling reason to fit it.
        # df_beta=f*(-np.log(1.0+(e_break/e)**((gamma1-gamma2)/self.beta))
        #            +((gamma1-gamma2)/self.beta)*np.log(e_break/e)/bottom)

        return np.asarray([df_n0,df_gamma1,df_gamma2,df_break])

    def full_name(self):
        return '%s, e0=%.0f, beta=%.3g'% (self.pretty_name,self.e0,self.beta)

#===============================================================================================#

class DoublePowerLaw(Model):
    """ Spectral model is the sum of two indepedendent power laws.  E.g., the Crab Nebula = IC + synch.

        Spectral parameters:
        
            n0             differential flux at e0 MeV for first power law
            gamma1        (absolute value of) spectral index for first power law
            gamma2        (absolute value of) spectral index for second power law
            ratio         ratio of the differential fluxes of first and second power law at e0
    """
    def __call__(self,e):
        n0,gamma1,gamma2,ratio=10**self._p
        return n0*((self.e0/e)**gamma1 + ratio*(self.e0/e)**gamma2)

#===============================================================================================#

class DoublePowerLawCutoff(Model):
    """ Spectral model is the sum of two indepedendent power laws, one with a cutoff.  E.g., a pulsar + PWnp.

        Spectral parameters:
        
            n0             differential flux at e0 MeV for first power law
            gamma1        (absolute value of) spectral index for first power law
            gamma2        (absolute value of) spectral index for second power law
            cutoff        cutoff -- note goes with gamma!
            ratio         ratio of the differential fluxes of first and second power law at e0
    """
    def __call__(self,e):
        n0,gamma1,gamma2,cutoff,ratio=10**self._p
        return n0*((self.e0/e)**gamma1*np.exp(-e/cutoff) + ratio*(self.e0/e)**gamma2)


#===============================================================================================#


class LogParabola(Model):
    """ Implement a log parabola (for blazars.)  See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            alpha        (absolute value of) constant spectral index
            beta         co-efficient for energy-dependent spectral index
            e_break     break energy
    """
    def __call__(self,e):
        n0,alpha,beta,e_break=10**self._p
        #alpha -= self.index_offset
#        return n0*(e_break/e)**(alpha - beta*np.log(e_break/e))
        x = np.log(e_break/e)
        y = (alpha - beta*x)*x
        return n0*np.exp(y) #np.clip(y, -10, 100)) #protect over, underflows

    def gradient(self,e):
        n0,alpha,beta,e_break=10**self._p
        #alpha -= self.index_offset
        #f  = n0*(e_break/e)**(alpha - beta*np.log(e_break/e))
        #log_term = np.log(e_break/e)
        #return np.asarray([f/n0,f*log_term,-f*log_term**2,f*alpha/e_break])
        x =np.log(e_break/e)
        y = (alpha - beta*x)*x
        f = n0*np.exp(y) # np.clip(y, -10, 100))
        return np.asarray([f/n0, f*x, -f*x**2, f*alpha/e_break])

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
        return ebreak*np.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ set a new break energy, adjusting the norm parameter """
        ebreak = 10** self._p[3]
        gamma = 10** self._p[1]
        self._p[0] += gamma * np.log10(ebreak/e0p)
        self._p[3] = np.log10(e0p)
 
    def create_powerlaw(self, beta_max=3e-2):
        """ if beta is small and fixed, return an equivalent PowerLaw, otherwise just return self """
        if self[2]>beta_max or self.free[2]: return self
        nm = PowerLaw(p=self[0:2], e0=self[3])
        nm.cov_matrix=self.cov_matrix[:-2,:-2]
        return nm
    
    @property
    def e0(self):
        return 10**self._p[3]
        
    @property
    def eflux(self):
        n0, alpha,beta, ebreak = 10**self._p
        return n0 * ebreak**2
#===============================================================================================#

class ExpCutoff(Model):
    """ Implement a power law with exponential cutoff.  See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            gamma        (absolute value of) spectral index
            cutoff      e-folding cutoff energy (MeV)
    """
    def __call__(self,e):
        n0,gamma,cutoff=10**self._p
        return n0* (self.e0/e)**gamma * np.exp(-e/cutoff)

    def gradient(self,e):
        n0,gamma,cutoff = 10**self._p
        f = n0* (self.e0/e)**gamma * np.exp(-e/cutoff)
        return np.asarray([f/n0,f*np.log(self.e0/e),f*e/cutoff**2])

    def pivot_energy(self):
        """ assuming a fit was done, estimate the pivot energy 
              
        """
        A  = 10**self._p[0]
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('%s fit required before calculating pivot energy' %self.name)
        return self.e0*np.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ set a new reference energy, adjusting the norm parameter """
        gamma = 10** self._p[1]
        self._p[0] += gamma * np.log10(self.e0/e0p)
        self.e0 = float(e0p) 
        
    def full_name(self):
        return '%s, e0=%.0f'% (self.pretty_name,self.e0)
 
    @property
    def eflux(self):
        n0 = 10**self._p[0]
        return n0 * self.e0**2
    
#===============================================================================================#

class ExpCutoffPlusPL(Model):
    """ Implement a power law with exponential cutoff + an additional power law.  A la pulsar + PWnp.

        Spectral parameters:
        
            n0_1         differential flux at e0 MeV
            gamma_1     (absolute value of) spectral index
            cutoff_1    e-folding cutoff energy (MeV)
            n0_2
            gamma_2
    """
    def __call__(self,e):
        n0_1,gamma_1,cutoff_1,n0_2,gamma_2 = 10**self._p
        return n0_1*(self.e0/e)**gamma_1*np.exp(-e/cutoff_1) + n0_2*(self.e0/e)**gamma_2

#===============================================================================================#

class AllCutoff(Model):
    """ Implement an exponential cutoff.  This for the case when cutoff too low to constrain index.
        See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            cutoff      e-folding cutoff energy (MeV)
    """
    def __call__(self,e):
        n0,cutoff=10**self._p
        if cutoff < 0: return 0
        return n0*np.exp(-e/cutoff)

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
        return n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)

    def gradient(self,e):
        n0,gamma,cutoff,b = 10**self._p
        f = n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)
        return np.asarray([f/n0,f*np.log(self.e0/e),
                     f*(b/cutoff)*(e/cutoff)**b,f*(e/cutoff)**b*np.log(cutoff/e)])

    def pivot_energy(self):
        """ assuming a fit was done, estimate the pivot energy 
              
        """
        if self._p[3]!=0. : print "WARNING: b is not 1, the pivot energy computation might be inaccurate"
        A  = 10**self._p[0]
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('%s fit required before calculating pivot energy' %self.name)
        return self.e0*np.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ set a new reference energy, adjusting the norm parameter """
        gamma = 10** self._p[1]
        self._p[0] += gamma * np.log10(self.e0/e0p)
        self.e0 = float(e0p) 
#===============================================================================================#

class CompositeModel(Model):
    """ A model which joins other models. Subclasses must
        implement the __call__(e) function which says
        how to join the models. """

    def __init__(self,*models, **kwargs):
        iscopy = kwargs.pop('iscopy', False)
        if len(models) < 1:
            raise Exception("CompositeModel must be created with more than one spectral model")
        for m in models:
            if not isinstance(m,Model):
                raise Exception("CompositeModel must be created with a list of models.")

        self.flux_scale = 1.
        self.models = models
        self.cov_matrix = np.zeros([self.npar,self.npar]) #default covariance matrix
        self.e0 = 1000. # not sure why, but sed_plotter needs this

    def __call__(self,e): 
        """ Must be implemented by a subclass. """
        pass

    @property
    def param_names(self):
        return reduce(operator.add,[i.param_names for i in self.models])

    @property
    def pretty_name(self):
        return self.operator.join([i.pretty_name for i in self.models])

    @property
    def background(self):
        """ Seems like a reasonable test. """
        return np.any([i.background for i in self.models])

    @background.setter
    def background(self,b):
        for m in self.models: m.background = b

    @property
    def n(self): 
        return np.asarray([len(i._p) for i in self.models])

    @property
    def npar(self): 
        return sum([len(i._p) for i in self.models])

    @property
    def _p(self): 
        return np.append(*[i._p for i in self.models])

    @_p.setter
    def _p(self, value):
        assert(len(self._p) == len(value))
        counter=0
        for i in xrange(len(self.n)):
            self.models[i]._p = value[counter:counter+self.n[i]]
            counter += self.n[i]

    @property
    def free(self): 
        return np.append(*[i.free for i in self.models])

    @free.setter
    def free(self, value):
        assert(len(self.free) == len(value))
        counter=0
        for i in xrange(len(self.n)):
            self.models[i].free = value[counter:counter+self.n[i]]
            counter += self.n[i]
      
    def setp(self, i, par, internal=False):
        """ set internal value, convert unless internal
        """
        i=self.mapper(i) #?
        if not internal: 
            assert par>0, 'Model external parameter cannont be negative'
        counter=0
        for k in xrange(len(self.n)):
            if counter<i:
                counter += self.n[k]
                continue
            self.models[k]._p[i-counter] = par if internal else  np.log10(par)
            return       

    def set_parameters(self,new_vals):
        """Set FREE internal parameters; new_vals should have length equal to number of free parameters."""
        assert len(new_vals)==(self.free).sum(), 'attempt to set wrong number of free parameters, model %s' %self
        pars = self._p
        pars[self.free] = new_vals.astype('f')
        self._p = pars
 

class FrontBackConstant(CompositeModel):
    """ Composite model that is either/or, for front or back
        select which constant based on value (0 or 1) of ct
    """
    name = 'FrontBackConstant'
    operator='+'
    def __init__(self, f=1, b=1, **kwargs):
        super(FrontBackConstant, self).__init__(Constant(),Constant(), **kwargs)
        self.models[0].param_names=['Scale_front']
        self.models[1].param_names=['Scale_back']
        self.models[0][0]=f
        self.models[1][0]=b
        self.ct = 0
        
    def __call__(self, e):
        return self.models[self.ct](e)
        
    def gradient(self, e):
        return np.hstack([(1-self.ct)*self.models[0].gradient(e).T, 
                          self.ct*self.models[1].gradient(e).T])

#===============================================================================================#

class SumModel(CompositeModel):
    """ Model is the sum of other models. 
    
        Easy to create: 
            >>> m1=PowerLaw(index=2);m1.set_flux(1)
            >>> m2=LogParabola(beta=2); m1.set_flux(1)

        The model is just the sum of the simpler modesl
            >>> sum_model=SumModel(m1,m2)
            >>> energy = np.asarray([1e2,1e3,1e4])
            >>> np.allclose(sum_model(energy), m1(energy) + m2(energy))
            True
    """
    operator = '+'
    name = 'SumModel'

    def __call__(self,e):
        return np.array([model(e) for model in self.models]).sum(axis=0)

    def gradient(self,e):
        """ Assume all models have a gradient! """
        return np.append([i.gradient(e) for i in self.models])

    def set_flux(self,flux,*args,**kwargs):
        """ Set the flux of the source. 

            Makes sure relative normalization of each component is unchanged
        """
        change = np.log10(flux/self.i_flux(*args,**kwargs))
        for m in self.models: m._p[0] += change

#===============================================================================================#

class ProductModel(CompositeModel):
    """ Model is the product of other Models.

        Easy to create: 

            >>> m1=PowerLaw(index=2);m1.set_flux(1)
            >>> m2=LogParabola(beta=2); m1.set_flux(1)
            >>> prod_model=ProductModel(m1,m2)

        The model is just the product of the simpler modesl

            >>> energy = np.asarray([1e2,1e3,1e4])
            >>> np.allclose(prod_model(energy), m1(energy) * m2(energy))
            True

        Note: set_flux function just adjusts normalziation of 
              first parameter, which should be correct! """
    operator = '*'
    name = 'ProductModel'

    def __call__(self,e):
        return np.array([model(e) for model in self.models]).prod(axis=0)

    def gradient(self,e):
        """ Assume all models have a gradient! """
        return np.append([i.gradient(e)/i.__call__(e) for i in self.models])*self.__call__(e)

#===============================================================================================#

class Constant(Model):
    def __call__(self,e):
        return np.ones_like(e)*10**self._p[0]
    
    def fast_iflux(self,emin=100,emax=1e6):
        return (emax-emin)*10**self._p[0]

    def gradient(self,e):
        return  np.array([np.ones_like(e)])

#===============================================================================================#

class InterpConstants(Model):

    def __call__(self,e):
        interp = interp1d(self.e_breaks,10**self._p)
        return interp(np.log10(e))

    def set_flux(self,flux,**kwargs):
        raise NotImplementedError("No way to set flux for InterpConstants spectral model")

#===============================================================================================#

class FileFunction(Model):
    r""" Defines a spectral model from an ascii file with the same
        convention as gtlike. Each point is connected with a power-law
        interpolation.
        
        For example, we could create a file which defines a PowerLaw,
        but very coarsely:

            >>> model = PowerLaw()

        Save out the model to a file
            
            >>> from tempfile import NamedTemporaryFile
            >>> temp = NamedTemporaryFile()
            >>> filename = temp.name
            >>> e = np.logspace(0, 5, 3)
            >>> pdf = model(e)
            >>> temp.write('\n'.join('%s\t%s' % (i,j) for i,j in zip(e,pdf)))
            >>> temp.seek(0)

        Load the file into the FileFunction object

            >>> file_function = FileFunction(file=filename)

        And the power-law is restored, even between the saved points.

            >>> energies = np.logspace(0.122, 4.83, 42)
            >>> np.allclose(model(e), file_function(e))
            True
    """
    def __make_interp__(self):
        self.interp = interp1d(np.log10(self.energy),np.log10(self.flux),
                bounds_error=False,fill_value=-np.inf)

    def __init__(self,*args,**kwargs):

        super(FileFunction,self).__init__(*args,**kwargs)

        if not hasattr(self,'file'):
            raise Exception("FileFunction must be created with a file.")

        file=np.genfromtxt(self.file,unpack=True)
        self.energy,self.flux=file[0],file[1]

        self.__make_interp__()

    def __call__(self,e):
        return 10**self.interp(np.log10(e))
    
    def gradient(self,e):
        return np.asarray([])

    def __getstate__(self):
        """ You cannot pickle an interp1d object. """
        d=copy.copy(self.__dict__)
        del d['interp']
        return d

    def __setstate__(self,state):
        """ recreate interp1d object. """
        self.__dict__ = state
        self.__make_interp__()

class DMFitFunction(Model):
    """ Wrap gtlike's DMFitFunction interface. 
    
        N.B. The bug Sheridan reported that the set_flux function 
        was not working should now be fixed:
        
            >>> model = DMFitFunction()
            >>> model.set_flux(1e-7, emin=1e3, emax=1e5)
            >>> print '%g' % model.i_flux(emin=1e3, emax=1e5)
            1e-07

        Test the getters and setters

            >>> model['sigmav']=3.14
            >>> print '%g' % model['sigmav']
            3.14

        There was previously a bug in set_parameters, 
        lets see if its fixed:

            >>> model.set_parameters(np.log10([5,500]))
            >>> print '%g' % model['sigmav']
            5
            >>> print '%g' % model['mass']
            500

        Note, the parameters which are not directly fit (like bratio) get set correctly:

            >>> model = DMFitFunction(bratio=2)
            >>> print model.dmf.getParam('bratio').getTrueValue()
            2.0
            >>> model = DMFitFunction(bratio=3)
            >>> print model.dmf.getParam('bratio').getTrueValue()
            3.0

        Test a few hard coded values, to make sure the function values are correct:

            >>> model = DMFitFunction(sigmav=1e-26, mass=100,
            ... channel0=4, channel1=1, bratio=1, norm=2.5e17)

            >>> model =DMFitFunction(norm=2.5e17, sigmav=1e-26, channel0=4,channel1=1,mass=100,bratio=1.0)

        These points agree with the fortran code.

            >>> e = [1, 10, 100, 1000, 10000, 100000 , 1000000]
            >>> dnde = [ 9.55801576e-18, 2.04105211e-16,  4.43719263e-16, 1.00123992e-16, 1.44911940e-18, 0.0, 0.0 ]
            >>> print np.allclose(model(e), dnde)
            True
    """
    def full_name(self):
        return '%s, norm=%.1f, bratio=%.1f channel0=%d, channel1=%d' % (self.pretty_name,
                                                                        self.norm, self.bratio, 
                                                                        self.channel0, self.channel1)

    def __getstate__(self):
        d=copy.copy(self.__dict__)
        del d['dmf']
        return d

    def __setstate__(self,state):
        self.__dict__ = state
        self._update()

    def _update(self):
        """ Update the DMFitFunction internally.
            This function should be called
            automatically when necessary.
        """
        if not hasattr(self,'dmf'):
            import pyLikelihood
            self.dmf=pyLikelihood.DMFitFunction()

        for i,param_name in enumerate(self.param_names):
            self.dmf.setParam(param_name,self[param_name])

        # Set the parameters which are not fixed explicitly
        self.dmf.setParam('norm',self.norm)
        self.dmf.setParam('bratio',self.bratio)
        self.dmf.setParam('channel0',self.channel0)
        self.dmf.setParam('channel1', self.channel1)

    def __init__(self,  *args, **kwargs):
        import pyLikelihood

        # the DMFitFunction must exist before __init__ is called because
        # the __init__ will call setp().
        self.dmf=pyLikelihood.DMFitFunction()
        super(DMFitFunction,self).__init__(*args,**kwargs)

        # unbound all parameters in gtlike
        for n in np.append(self.param_names,['norm','bratio','channel0','channel1']):
            self.dmf.getParam(n).setBounds(-float('inf'),float('inf'))

        from SpatialModels import SpatialMap
        self.dmf.readFunction(SpatialMap.expand(self.file))
        self._update() # update all parameters in DMFitFunction

    def setp(self, *args, **kwargs):
        super(DMFitFunction,self).setp(*args, **kwargs)
        self._update()

    def set_parameters(self, *args, **kwargs):
        super(DMFitFunction,self).set_parameters(*args, **kwargs)
        self._update()

    def set_all_parameters(self, *args, **kwargs):
        super(DMFitFunction,self).set_all_parameters(*args, **kwargs)
        self._update()

    def __call__(self,e):
        """ Return energy in MeV. This could be vectorized. """
        from pyLikelihood import dArg
        if isinstance(e,collections.Iterable):
            return np.asarray([self.dmf(dArg(i)) for i in e])
        else:
            return self.dmf(dArg(e))


def convert_exp_cutoff(model):
    """ this function is need for XML parsing. """
    if model.name != 'ExpCutoff':
        raise Exception,'Cannot process %s into PLSuperExpCutoff'%(model.name)
    nm = PLSuperExpCutoff()
    nm._p    = np.append(model._p,0)
    nm.free = np.append(model.free,False)
    nm.cov_matrix[:,:] = 0
    nm.cov_matrix[:-1,:-1] = model.cov_matrix[:,:]
    nm.e0 = model.e0
    return nm
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
