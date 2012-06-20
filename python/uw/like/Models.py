"""A set of classes to implement spectral models.

    $Header: /nfs/slac/g/glast/ground/cvs/ScienceTools-scons/pointlike/python/uw/like/Models.py,v 1.103 2012/06/15 01:29:38 lande Exp $

    author: Matthew Kerr, Joshua Lande
"""
import os
import copy
import operator
import numpy as N
import numpy as np
from abc import abstractmethod

from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy import roots, optimize

from uw.utilities.parmap import LinearMapper, LogMapper


class Model(object):
    """ Spectral model giving dN/dE for a point source.  

        Default units are ph/cm^2/s/MeV.
    """
    def __init__(self, p=None, free=None, **kwargs):
        """ 

            Optional keyword arguments:

                =========    =======================================================
                Keyword      Description
                =========    =======================================================
                p              [p1,p2,...] default values of spectral parameters; see docstring individual model classes
                free           [True, True,...] a boolean list the same length as p giving the free (True) and fixed (False) parameters
                =========    =======================================================

            The input to this object should be fairly flexible:
                
                >>> model = PowerLaw(e0=100)
                >>> print model.e0
                100

            Note that pointlike internally stores parameters using a mapping:

                >>> model = Constant(scale=1, mappers=[LogMapper])
                >>> model['scale']
                1.0
                >>> model.getp('scale',internal=True)
                0.0
                >>> model.setp('scale',1,internal=True)
                >>> model.getp('scale',internal=True)
                1.0
                >>> model.getp('scale')
                10.0
                >>> model = Constant(scale=1, mappers=[LinearMapper])
                >>> model['scale']
                1.0
                >>> model.getp('scale',internal=True)
                1.0
        """
        self.background = False
        self.name = self.pretty_name = self.__class__.__name__

        self.npar = len(self.param_names)
        self.internal_cov_matrix = np.zeros([self.npar,self.npar]) #default covariance matrix

        if free is None:
            free = np.ones(self.npar)
        self.free = np.asarray(free, dtype=bool)

        if kwargs.has_key('mappers'):
            self.mappers = kwargs.pop('mappers')
        else:
            self.mappers = self.default_mappers

        for k,v in self.default_extra_params.items(): 
            setattr(self,k,v)

        self._p = np.empty(self.npar, dtype=float)
        if p is None:
            p = self.default_p
        self.set_all_parameters(p)

        # deal with other kwargs. If they are in the list of spectral
        # parameters, update the corresponding value. Tag them to 
        # the class.
        for k,v in kwargs.items():
            if k in self:
                # if k is a param name
                self[k]=v
            else:
                raise Exception("Unable to set parameter unknown parameter %s"  % k)
    
    def len(self): 
        """ Exmaple:
                >>> model = PowerLaw()
                >>> print model.len()
                2
        """
        return self.npar

    def __contains__(self, item):
        """ Example:
                >>> model = PowerLaw()
                >>> print 'Norm' in model
                True
                >>> print 'Normalization' in model
                False
        """
        try:
            self.__getitem__(item)
            return True
        except Exception, ex:
            return False
        
    def __getitem__(self, index):
        """ Wrapper around getp:
                
                >>> model=PowerLaw()
                >>> print model['index']
                2.0
        """
        return self.getp(index)
        
    def __setitem__(self, index, value):
        self.setp(index,value)
        
    def getp(self, i, internal=False):
        """ Get external value for parameter # i 

                >>> model=PowerLaw()
                >>> print model['index']
                2.0
        
            You can also get the non-fit extra_params using getp:
                
                >>> print model['e0']
                1000.0

            """
        if i in self.default_extra_params.keys():
            return getattr(self,i)
        i=self.name_mapper(i)
        return self._p[i] if internal else self.mappers[i].toexternal(self._p[i])
    
    def get_all_parameters(self, internal=False):
        """ Get a copy of the full set of parameters (external representation)"""
        if internal:
            return self._p.copy() 
        else:
            return np.asarray(map(self.getp,np.arange(self.npar)))
    
    def set_all_parameters(self, pars, internal=False):
        """ set all parameters (external representation)"""
        assert len(pars)==self.npar
        pars  = np.asarray(pars, dtype=float)
        if internal:
            self._p = pars
        else:
            for i,p in enumerate(pars):
                self.setp(i,p)

    def setp(self, i, par, internal=False):
        """ set internal value, convert unless internal. 
            
                >>> model=PowerLaw()
                >>> model['index'] = 1
                >>> print model['index'] 
                1.0

            You can also set the non-fit extra_params using setp:

                >>> model['e0'] = 100
                >>> print model.e0
                100
        
        """
        if i in self.default_extra_params.keys():
            return setattr(self,i,par)
        i=self.name_mapper(i)
        self._p[i] = par if internal else self.mappers[i].tointernal(par)
        
    def get_parameters(self):
        """ Return FREE internal parameters ; used for spectral fitting."""
        return self._p[self.free]

    def set_parameters(self,new_vals):
        """ Set FREE internal parameters; new_vals should have length equal to number of free parameters."""
        assert(len(new_vals)==(self.free).sum())
        self._p[self.free] = new_vals.astype(float) # downcast to float needed?

    def name_mapper(self,i):
        """ This object takes in a parameter and maps it to an index.
            Currently, this function is not particularly useful, but it
            could be generalized to allow lazier parameter selection. 

            For example:

                >>> model = PowerLaw()
                >>> model.name_mapper('norm')
                0
                >>> model.name_mapper('index')
                1
            """
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

            For example:

                >>> model = PowerLaw()
                >>> print model.free
                [ True  True]
                >>> model.freeze('index')
                >>> print model.free
                [ True False]
            """
        i=self.name_mapper(i)
        self.free[i] = not freeze

    def thaw(self,i): self.freeze(i,freeze=False)

    def set_cov_matrix(self,new_cov_matrix):
        self.internal_cov_matrix[np.outer(self.free,self.free)] = np.ravel(new_cov_matrix)

    def dexternaldinternal(self):
        """ Return the derivative of the external parameters with repsect to the internal parametesr. """
        return np.asarray([m.dexternaldinternal(p) for m,p in zip(self.mappers,self.get_all_parameters())])

    def get_cov_matrix(self):
        """Return covariance matrix transformed out of log space."""
        p = self.dexternaldinternal()
        pt=p.reshape((p.shape[0],1))
        return p*self.internal_cov_matrix*pt

    def get_free_errors(self):
        """Return the diagonal elements of the (log space) covariance matrix for free parameters."""
        return np.diag(self.internal_cov_matrix)[self.free]**0.5

    def corr_coef(self):
        """Return the linear correlation coefficients for the estimated covariance matrix."""
        sigmas = np.diag(self.internal_cov_matrix)**0.5
        return self.internal_cov_matrix / np.outer(sigmas,sigmas)

    def statistical(self,absolute=False,two_sided=False):
        """ Return the parameter values and fractional statistical errors.
            If no error estimates are present, return 0 for the fractional error.

            two_sided : bool
                if set, return 3 tuples: values, +errors, -errors    
        """
        p = self.get_all_parameters()
        if np.all(self.internal_cov_matrix==0) :
            # covariance matrix is invalid
            z = np.zeros_like(p)
            return (p,z,z) if two_sided else (p,z) 
        elif not two_sided:
            vars = np.diag(self.get_cov_matrix())
            vars[vars<0]=0
            errs = vars**0.5 
            return p,errs/(1 if absolute else p)
        else:
            vars = np.diag(self.internal_cov_matrix)
            vars[vars<0]=0
            errs = vars**0.5
            lo_rat = p - np.asarray([self.mappers[i].toexternal(self.getp(i,internal=True) - errs[i]) 
                                     for i in range(self.npar)])
            hi_rat = np.asarray([self.mappers[i].toexternal(self.getp(i,internal=True) + errs[i]) 
                                 for i in range(self.npar)]) - p
            return p,hi_rat/(1 if absolute else p),lo_rat/(1 if absolute else p)

    @abstractmethod
    def external_gradient(self,e):
        """ Return gradient with respect to external parameters. """
        pass

    def gradient(self,e):
        """ Return gradient with respect to internal parameters. """
        dextdint=self.dexternaldinternal()
        dextdint=dextdint.reshape((len(dextdint),1))
        return self.external_gradient(e)*dextdint

    def error(self,i):
        """ Get the EXTERNAL error for parameter i """
        i=self.name_mapper(i)
        return (np.diag(self.get_cov_matrix())**0.5)[i]

    def set_error(self,i,error):
        """ Set the INTERNAL error for parameter i 

                >>> model = PowerLaw(mappers=[LogMapper,LogMapper])
                >>> print model.error('norm')
                0.0
                >>> model.set_error('norm',1)
                >>> print model.error('norm')
                1.0

                >>> model = PowerLaw(mappers=[LinearMapper,LinearMapper])
                >>> model.set_error('norm',1)
                >>> print model.error('norm')
                1.0
        """
        i=self.name_mapper(i)

        # delete any covariance
        self.internal_cov_matrix[i,:]=self.internal_cov_matrix[:,i]=0
        dextdint=self.mappers[i].dexternaldinternal(self.getp(i))
        dintdext=1/dextdint
        self.internal_cov_matrix[i,i]=(error*dintdext)**2

    def __str__(self,absolute=False, indent=''):
        """ Return a pretty print version of parameter values and errors.
            indent: string to prepend to each line (must be called explicitly)

                >>> model = PowerLaw()
                >>> print model
                Norm      : 1e-11 
                Index     : 2 
                Ph. Flux  : 9.99e-08 (DERIVED)
                En. Flux  : 1.11e-10 (DERIVED)
        """
        p,hi_p,lo_p = self.statistical(absolute=absolute,two_sided=True)
        if not self.background:
            if not np.all(self.internal_cov_matrix==0):
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
                if i < self.npar:
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
                if i < self.npar:
                    frozen = '' if self.free[i] else '(FROZEN)'
                else:
                    frozen = '(DERIVED)'
                l+=[t_n+': %.3g %s'%(p[i],frozen)]
            return ('\n'+indent).join(l)
            
    def __repr__(self): return self.__str__()

    def i_flux(self,emin=100,emax=1e5,e_weight=0,cgs=False,error=False,two_sided=False, quiet=False):
        """ Return the integral flux, \int_{emin}^{emax} dE E^{e_weight} dN/dE.
            e_weight = 0 gives the photon flux (ph cm^-2 s^-1)
            e_weight = 1 gives the energy flux (MeV cm^-2 s^-1) (see kwargs)

            Optional keyword arguments:

                =========    =======================================================
                Keyword      Description
                =========    =======================================================
                emin         [100] lower bound in MeV
                emax         [1e5] upper bound in MeV
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
            epsabs = min(func(emin),func(emax))*1e-10 # needed since epsrel does not seem to work
            flux    =  units*quad(func,emin,emax,epsabs=epsabs,full_output=True)[0]
            if error:
                # will silently ignore 'free' parameters without errors
                mask = (self.free) * (self.internal_cov_matrix.diagonal()>0)
                args = (emin,emax,e_weight,cgs,False)
                d    = self.__flux_derivs__(*args)[mask]
                dt  = d.reshape( (d.shape[0],1) ) #transpose
                err = (d * self.internal_cov_matrix[mask].transpose()[mask] * dt).sum()**0.5
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

    def save_profile(self, filename, emin, emax, numpoints=200, clip_ends=True, clip_fraction=1e-20):
        """ Create a textfile of the spectral model. Useful for running
            gtobssim or gtlike with weird spectra. Also useful
            for creating a numeric representation of a spectral
            model using the FileFunction object. 

            if clip_ends is True, the file will clip spectral points a fraction clip_fraction less than
            the maximum at either end of the spectrum. This is to avoid having numerically small (or 0)
            values in the spectral file which can cause trouble in gtobssim, which does a powerlaw
            interpolation of the points.
        """
        energies=np.logspace(np.log10(emin),np.log10(emax),numpoints)
        fluxes=self(energies)

        if clip_ends:
            max_flux = max(fluxes)
            clip_flux = clip_fraction*max_flux
            # Clip out 0s on either end of spectral file. Otherwise gtobbsim gets angry.
            while fluxes[0] < clip_flux: energies,fluxes=energies[1:],fluxes[1:]
            while fluxes[-1] < clip_flux: energies,fluxes=energies[:-1],fluxes[:-1]

        # Paranoid check:
        if np.any(fluxes==0): raise Exception("Error: 0s found in differential flux")

        open(filename,'w').write('\n'.join(['%g\t%g' % (i,j) for i,j in zip(energies,fluxes)]))

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
        self.setp(0, 1) # First, set prefactor to 1
        new_prefactor = flux/self.i_flux(*args,**kwargs)
        self.setp(0,new_prefactor)

    def copy(self): return copy.deepcopy(self)

    def fast_iflux(self,emin=100,emax=1e5):
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
            # reset parameters
            hi.setp(i,hi._p[i] - errs[i],internal=True)
            lo.setp(i,lo._p[i] + errs[i],internal=True)

        return np.asarray(derivs)
        
    def full_name(self):
        return self.pretty_name
        
    @abstractmethod
    def pivot_energy(self):
        pass

    def zero(self):
        """ Set source's flux to 0:

                >>> model = PowerLaw()
                >>> print model.iszero()
                False
                >>> print model['norm']
                1e-11
                >>> model.zero()
                >>> print model.iszero()
                True
                >>> print model['norm']
                0.0
                >>> model.unzero()
                >>> print model['norm']
                1e-11
                >>> print model.iszero()
                False
        """
        if self.iszero(): return #already zeroed
        self.old_flux = self.getp(0) # get the internal value
        self.setp(0, 0)
        self.old_free = self.free.copy()
        self.free[:] = False
    
    def iszero(self):
        return self.getp(0)==0

    def unzero(self):
        self.setp(0, self.old_flux)
        self.free = self.old_free.copy()

    def set_prefactor(self,  prefactor, energy):
        """ Set the prefactor of the source (at a given energy).
            Useful if you want to specify dN/dE for a source at a 
            particular energy. 

                >>> model = PowerLaw(e0=1000)
                >>> model.set_prefactor(1e-10, 100)
                >>> print model(100)
                1e-10
        """
        self.setp(0,1) #set prefactor to 1
        self.setp(0, prefactor/self(energy))

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
    default_p=[1e-11, 2.0]
    default_extra_params=dict(e0=1e3)
    param_names=['Norm','Index']
    default_mappers=[LogMapper,LinearMapper]

    def __call__(self,e):
        n0,gamma=self['Norm'],self['Index']
        return n0*(self.e0/e)**gamma

    def fast_iflux(self,emin=100,emax=1e6):
        n0,gamma=self['Norm'],self['Index']
        return n0/(1-gamma)*self.e0**gamma*(emax**(1-gamma)-emin**(1-gamma))

    def external_gradient(self,e):
        n0,gamma=self['Norm'],self['Index']
        f = n0*(self.e0/e)**gamma
        return np.asarray([f/n0,f*np.log(self.e0/e)])

    def pivot_energy(self):
        """ assuming a fit was done, estimate the pivot energy 
              (This only implemented for this model)
        """
        A  = self['Norm']
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('PowerLaw fit required before calculating pivot energy')
        return self.e0*np.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ Set a new reference energy, adjusting the norm parameter 
        
            >>> model1=PowerLaw(e0=1e2)
            >>> model2=model1.copy()
            >>> model2.set_e0(1e3)
            >>> np.allclose(model1(np.logspace(1,7,8)),model2(np.logspace(1,7,8)))
            True
        """
        gamma = self['Index']
        self['norm'] *= (e0p/self.e0)**-gamma
        self.e0 = e0p
        
    def full_name(self):
        return '%s, e0=%.0f'% (self.pretty_name,self.e0)
    @property
    def eflux(self):
        return self.e0**2*self['Norm']

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
    default_p=[1e-7 , 2.0]
    default_extra_params=dict(emin=100, emax=1e6)
    param_names=['Int_Flux','Index']
    default_mappers=[LogMapper,LinearMapper]

    def __call__(self,e):
        flux,gamma=self.get_all_parameters()
        return (flux*(1-gamma)/(self.emax**(1-gamma)-self.emin**(1-gamma)))*e**(-gamma)

    def fast_iflux(self,emin=100,emax=np.inf):
        flux,gamma=self.get_all_parameters()
        return flux*(emax**(1-gamma) - emin**(1-gamma)) / (self.emax**(1-gamma) - self.emin**(1-gamma))

    def external_gradient(self,e):
        flux,gamma=self.get_all_parameters()
        e1 = self.emax; e0 = self.emin
        d1 = e1**(1-gamma)
        d0 = e0**(1-gamma)
        t  = (np.log(e0)*d0 - np.log(e1)*d1)/(d1 - d0)
        f  = (flux*(1-gamma)/(d1-d0))*e**(-gamma)
        return np.asarray([f/flux,-f*(np.log(e) + (1./(1-gamma) + t))])

    def full_name(self):
        return '%s, emin=%.0f emax=%.0f'% (self.pretty_name,self.emin,self.emax)

class BrokenPowerLaw(Model):
    """ Implement a broken power law.  See constructor docstring for further keyword arguments.

        Spectral parameters:

            n0            differential flux at e0 MeV
            gamma1      (absolute value of) spectral index for e < e_break
            gamma2      (absolute value of) spectral index for e > e_break
            e_break     break energy
    """
    default_p=[1e-11, 2.0, 2.0 ,1e3]
    default_extra_params=dict()
    param_names=['Norm','Index_1','Index_2', 'E_break']
    default_mappers=[LogMapper,LinearMapper,LinearMapper,LogMapper]

    def __call__(self,e):
        n0,gamma1,gamma2,e_break=self.get_all_parameters()
        return n0*(e_break/e)**np.where(e<e_break,gamma1,gamma2)

    def external_gradient(self,e):
        n0,gamma1,gamma2,e_break=self.get_all_parameters()
        mask = e < e_break
        x = e_break/e
        lx = np.log(x)
        g = np.where(mask,gamma1,gamma2)
        f = n0*x**g
        return np.asarray([f/n0,f*lx*mask,f*lx*(~mask),f/e_break*g])

class BrokenPowerLawFlux(Model):
    """ Similar to PowerLawFlux for BrokenPowerLaw spectrum, the integral 
        flux is the free parameter rather than the Prefactor.

        Spectral parameters:

           flux      integrated flux from emin to emax MeV
           gamma1  (absolute value of) spectral index for e < e_break
           gamma2  (absolute value of) spectral index for e > e_break
           e_break break energy
    """
    default_p=[1e-7, 2.0, 2.0 , 1e3]
    default_extra_params=dict(emin=100, emax=1e6)
    param_names=['Int_Flux','Index_1','Index_2', 'E_break']
    default_mappers=[LogMapper,LinearMapper,LinearMapper,LogMapper]

    def __call__(self,e):
        flux,gamma1,gamma2,e_break=self.get_all_parameters()
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

class BrokenPowerLawCutoff(Model):
    """ Implement a broken power law.  See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            gamma1      (absolute value of) spectral index for e < e_break
            gamma2      (absolute value of) spectral index for e > e_break
            e_break     break energy
    """
    default_p=[1e-11,2,2,1e3,3e3]
    default_extra_params=dict()
    param_names=['Norm','Index_1','Index_2','E_break','Cutoff']
    default_mappers=[LogMapper,LinearMapper,LinearMapper,LogMapper]

    def __call__(self,e):
        n0,gamma1,gamma2,e_break,cutoff=self.get_all_parameters()
        return n0*np.where( e < e_break, (e_break/e)**gamma1, (e_break/e)**gamma2 )*np.exp(-e/cutoff)

class SmoothBrokenPowerLaw(Model):
    """ Implement a smoothed broken power law. This is similar to a broken 
        powerlaw but uses the parameter beta to smoothly interpolate between 
        the two powerlaws.

        The implemenation is exactly the same as the gtlike implementation:
            http://glast.stanford.edu/cgi-bin/viewcvs/Likelihood/src/SmoothBrokenPowerLaw.cxx

        but with the indices the negative of the gtlike indices for internal consistency.
            

        Everything is defined exactly the same as gtlike except for the 
        usual replacement that gamma1 and gamma2 are defined to be the 
        negative of the gtlike definition. 

        Note that this function has a somewhat confusing feature that when
        beta>0, the function always becomes softer for energies greater then
        the break. For beta<0, the function always becomes harder for energies
        greater than the break. So you need to pick beta in advance depending
        upon whether you want your function to slope up or down.

        For now, I am not allowing beta to be fit. I am not sure anybody ever fits its.

        Also note that e0 is used in the spectral function and can be set (but 
        not fit) by passing in e0 when initialized. """

    default_p=[1e-11,2.0,2.0,1e3]
    default_extra_params=dict(beta=0.1, e0=1e3)
    param_names=['Norm','Index_1','Index_2','E_break']
    default_mappers=[LogMapper,LinearMapper,LinearMapper,LogMapper]

    def __call__(self,e):
        n0,gamma1,gamma2,e_break=self.get_all_parameters()

        # switch from pointlike to gltike convention
        gamma1=-gamma1
        gamma2=-gamma2

        # This is exactly the gtlike formula
        return n0*(e/self.e0)**gamma1*\
               (1+(e/e_break)**((gamma1-gamma2)/self.beta))**(-self.beta)

    def external_gradient(self,e):
        """ lots of derivatives. You can derive them all by hitting the 
            derivative on log(dN/dE). """

        n0,gamma1,gamma2,e_break=self.get_all_parameters()

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


class LogParabola(Model):
    """ Implement a log parabola (for blazars.)  See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            alpha        (absolute value of) constant spectral index
            beta         co-efficient for energy-dependent spectral index
            e_break     break energy
    """
    default_p=[1e-11, 2.0, 1e-5,2e3]
    default_extra_params=dict()
    param_names=['Norm','Index','beta','E_break']
    default_mappers=[LogMapper,LinearMapper,LinearMapper,LogMapper]

    def __call__(self,e):
        n0,alpha,beta,e_break=self.get_all_parameters()
        x = np.log(e_break/e)
        y = (alpha - beta*x)*x
        return n0*np.exp(y)

    def external_gradient(self,e):
        n0,alpha,beta,e_break=self.get_all_parameters()
        x =np.log(e_break/e)
        y = (alpha - beta*x)*x
        f = n0*np.exp(y) # np.clip(y, -10, 100))
        return np.asarray([f/n0, f*x, -f*x**2, f*alpha/e_break])

    def pivot_energy(self):
        """  
        Estimate the pivot energy
        Warning, there are potentially 3 real solutions to the pivot energy equation, and
        the code does not try to choose one : it returns all the real roots.
        """
        A  = 10**self._p[0]
        C = self.get_cov_matrix()
        if not self.free[2]:
            # the case beta fixed is equivalent to a PWL for the determination of the pivot energy.
            if C[1,1]==0:
                raise Exception('Models.LogParabola: Fit required before calculating pivot energy')
            ebreak = 10**self._p[3]
            return ebreak*np.exp( C[0,1]/(A*C[1,1]) )
        #the solution in the general case comes from the resolution of a cubic equation.
        #We assume here that Norm, alpha, and beta are the free parameter of the fit.
        #Indeed E_break should never be set to free anyway
        a= 2.*C[2,2]
        b= 3.*C[1,2]
        c= C[1,1]-2.*C[0,2]/A
        d= -C[0,1]/A
        results = roots([a,b,c,d])
        results = 10**self._p[3]*np.exp(results)
        print "Pivot energy solutions for the LogParabola : ",results
        return np.real(results[np.isreal(results)])

    def set_e0(self, e0p):
        """ set a new break energy, adjusting the Norm and Index parameter,
            so that the differential flux remains the same. Beta remains unchanged in this
            transformation.

                >>> raise Exception("...")
        """
        ebreak = 10** self._p[3]
        gamma = 10** self._p[1]
        beta = 10** self._p[2] 
        # there is a non-negative condition here to check on Index
        if e0p<ebreak*10**(-gamma/2./beta):
            e0p=ebreak*10**(-gamma/2./beta)
            print "Index>0 is not preserved by the transformation with input e0. Setting e0 to %g"%e0p
        x=np.log10(ebreak/e0p)
        self._p[0] += gamma * x + beta*x*x
        gamma -= 2*beta*x
        self._p[1] = np.log10(gamma)
        # beta is unchanged in this e0p->ebreak transformation
        self._p[3] = np.log10(e0p)
 
    def create_powerlaw(self, beta_max=3e-2):
        """ if beta is small and fixed, return an equivalent PowerLaw, otherwise just return self 

                >>> raise Exception("...")
        
        """
        if self[2]>beta_max or self.free[2]: return self
        nm = PowerLaw(p=self[0:2], e0=self[3], mappers=self.mappers[:-2])
        nm.internal_cov_matrix=self.internal_cov_matrix[:-2,:-2]
        return nm
    
    @property
    def e0(self): return self['e_break']
        
    @property
    def eflux(self):
        n0, alpha,beta, ebreak = self.get_all_parameters()
        return n0 * ebreak**2

class ExpCutoff(Model):
    """ Implement a power law with exponential cutoff.  See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            Norm       differential flux at e0 MeV
            Index      (absolute value of) spectral index
            Cutoff      e-folding cutoff energy (MeV)
    """
    default_p=[1e-11, 2.0, 2e3]
    default_extra_params=dict(e0=1e3)
    param_names=['Norm','Index','Cutoff']
    default_mappers=[LogMapper,LinearMapper,LogMapper]

    def __call__(self,e):
        n0,gamma,cutoff=self.get_all_parameters()
        return n0* (self.e0/e)**gamma * np.exp(-e/cutoff)

    def external_gradient(self,e):
        n0,gamma,cutoff=self.get_all_parameters()
        f = n0* (self.e0/e)**gamma * np.exp(-e/cutoff)
        return np.asarray([f/n0,f*np.log(self.e0/e),f*e/cutoff**2])

    def pivot_energy(self):
        """ assuming a fit was done, estimate the pivot energy 
              
        """
        A  = self['Norm']
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('%s fit required before calculating pivot energy' %self.name)
        return self.e0*np.exp( C[0,1]/(A*C[1,1]) )
        
    def set_e0(self, e0p):
        """ Set a new reference energy, adjusting the norm parameter 

                >>> model1 = ExpCutoff(e0=1e2)
                >>> model2 = model1.copy()
                >>> model2.set_e0(1e3)
                >>> model2.e0
                1000.0
                >>> np.allclose(model1(np.logspace(1,7,8)),model2(np.logspace(1,7,8)))
                True
        """
        gamma = self['Index']
        self['norm'] *= (e0p/self.e0)**-gamma
        self.e0 = e0p

        
    def full_name(self):
        return '%s, e0=%.0f'% (self.pretty_name,self.e0)
 
    @property
    def eflux(self):
        n0 = self['n0']
        return n0 * self.e0**2
    

class AllCutoff(Model):
    """ Implement an exponential cutoff.  This for the case when cutoff too low to constrain index.
        See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            cutoff      e-folding cutoff energy (MeV)
    """
    default_p=[1e-11, 1e3]
    default_extra_params=dict()
    param_names=['Norm','Cutoff']
    default_mappers=[LogMapper,LogMapper]

    def __call__(self,e):
        n0,cutoff=self.get_all_parameters()
        if cutoff < 0: return 0
        return n0*np.exp(-e/cutoff)

class PLSuperExpCutoff(Model):
    """Implement a power law with hyper-exponential cutoff.  See constructor docstring for further keyword arguments.

        Spectral parameters:
        
            n0            differential flux at e0 MeV
            gamma        (absolute value of) spectral index
            cutoff      e-folding cutoff energy (MeV)
            b             additional power in the exponent
        """
    default_p=[1e-11, 2.0, 2e3 ,1.]
    default_extra_params=dict(e0=1e3)
    param_names=['Norm','Index','Cutoff', 'b']
    default_mappers=[LogMapper,LinearMapper,LogMapper,LinearMapper]

    def __call__(self,e):
        n0,gamma,cutoff,b=self.get_all_parameters()
        return n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)

    def external_gradient(self,e):
        n0,gamma,cutoff,b=self.get_all_parameters()
        f = n0*(self.e0/e)**gamma*np.exp(-(e/cutoff)**b)
        return np.asarray([f/n0,f*np.log(self.e0/e),
                     f*(b/cutoff)*(e/cutoff)**b,f*(e/cutoff)**b*np.log(cutoff/e)])

    def pivot_energy(self):
        """  
        Assuming a fit was done, estimate the pivot energy. The parameter b is assumed to be fixed at 1.
        """
        if self._p[3]!=0. : print "WARNING: b is not 1, the pivot energy computation might be inaccurate"
        N0, gamma, Ecut, beta = self.get_all_parameters()
        E0 = self.e0
        C = self.get_cov_matrix()
        if C[1,1]==0:
            raise Exception('%s fit required before calculating pivot energy' %self.name)
        if not self.free[2]:
            # The case b and Ecut fixed is equivalent to a PWL for the determination of the pivot energy.
            return E0*np.exp( C[0,1]/(N0*C[1,1]) )
        # We assume here that Norm, gamma, and Ecut are the free parameter of the fit.
        # To estimate Epivot we need to minimize the relative uncertainty on the flux, i.e. we want that the function "func"=0. 
        # x=log(Epivot/E0)
        a = C[0,1]/N0
        b = C[1,1]
        c = E0**2./Ecut**(4.)*C[2,2]
        d = E0/N0/Ecut**2.*C[0,2]
        f = E0/Ecut**2.*C[1,2]
        func = lambda x : a + b*x + c*np.exp(2.*x) + np.exp(x)*(d-f*(1+x))
        result = optimize.newton(func, 1)
        result = E0*np.exp(result)
        print "Pivot energy solution (in MeV) for the PLSuperExpCutoff : ",result
        return result
        
    def set_e0(self, e0p):
        """ set a new reference energy, adjusting the norm parameter:

                >>> model1=PLSuperExpCutoff(e0=1e2)
                >>> model2=model1.copy()
                >>> model2.set_e0(1e3)
                >>> model1['Index'] = model2['Index']
                >>> np.allclose(model1(np.logspace(1,7,8)),model2(np.logspace(1,7,8)))
                True
        """
        gamma = self['Index']
        self['norm'] *= (e0p/self.e0)**-gamma
        self.e0 = e0p


class CompositeModel(Model):
    """ A model which joins other models. Subclasses must
        implement the __call__(e) function which says
        how to join the models. 
        
        CompositeModels should act just like regular models, implementing all
        of the parameters you expect for Models. 
        
        Here, we show how it works for a SumModel, the simplest CompositeModel.
        
        First, we create a SumModel

            >>> m1=PowerLaw(index=2);m1.set_flux(1)
            >>> m2=LogParabola(beta=2); m1.set_flux(1)
            >>> sum_model=SumModel(m1,m2)

        The sum_model has a param_names, mappers

            >>> sum_model.param_names == m1.param_names + m2.param_names
            True
            >>> sum_model.mappers == m1.mappers + m2.mappers
            True

        The CompositeModel.npar is the sum of parameters
            >>> sum_model.npar == m1.npar + m2.npar
            True

        And CompositeModel.npars is an array of the length of each parameter
            >>> np.all(sum_model.npars == [m1.npar,m2.npar])
            True

        And the parameters for the model is just the sum of the parameters for each model individually

            >>> np.all(sum_model.get_all_parameters() == np.append(m1.get_all_parameters(),m2.get_all_parameters()))
            True
            >>> np.all(sum_model._p == np.append(m1._p,m2._p))
            True

            >>> np.all(sum_model.free==np.append(m1.free,m2.free))
            True

        Also, getting and setting parameters individually acts just as one would expect

            >>> sum_model.mappers[1]
            <class 'uw.utilities.parmap.LinearMapper'>
            >>> sum_model[1] 
            2.0
            >>> sum_model.getp(1)
            2.0
            >>> print sum_model._p[1]
            2.0
            >>> sum_model[1] = -1
            >>> sum_model[1] 
            -1.0
            >>> sum_model.getp(1)
            -1.0
            >>> print sum_model._p[1]
            -1.0

        Also, the free parameters should be totally transparent:

            >>> sum_model.free[1]
            True
            >>> sum_model.freeze(1)
            >>> sum_model.free[1]
            False

        And after the individual parameter modifications, the CompositeModel._p
        and CompositeModel.free arrays are updated consistently in the individual
        Model objects and in the Composite object:

            >>> np.all(sum_model._p == np.append(m1._p,m2._p))
            True

            >>> np.all(sum_model.free==np.append(m1.free,m2.free))
            True

        Previously there was a bug with having 3 spatial parameters: 
            >>> c = SumModel(PowerLaw(),PowerLaw(),PowerLaw())
            >>> print c._p
            [-11.   2. -11.   2. -11.   2.]
            >>> print c.free
            [ True  True  True  True  True  True]
    """
    default_extra_params=dict() # Don't allow any of these

    def __init__(self, *models):
        if len(models) < 1:
            raise Exception("CompositeModel must be created with more than one spectral model")
        for m in models:
            if not isinstance(m,Model):
                raise Exception("CompositeModel must be created with a list of models.")

        self.name = self.__class__.__name__
        self.models = models
        self.internal_cov_matrix = np.zeros([self.npar,self.npar]) #default covariance matrix

    @abstractmethod
    def __call__(self,e): 
        """ Must be implemented by a subclass. """
        pass

    @property
    def param_names(self):
        return reduce(operator.add,[i.param_names for i in self.models])

    @property
    def mappers(self):
        return reduce(operator.add,[i.mappers for i in self.models])


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
    def npars(self): 
        return np.asarray([i.npar for i in self.models])

    @property
    def npar(self): 
        return sum([i.npar for i in self.models])

    @property
    def _p(self): 
        return np.concatenate([i._p for i in self.models])

    @_p.setter
    def _p(self, value):
        assert(len(self._p) == len(value))
        counter=0
        for i,n in enumerate(len(self.npars)):
            self.models[i]._p = value[counter:counter+n]
            counter += n

    @property
    def free(self): 
        return np.concatenate([i.free for i in self.models])

    @free.setter
    def free(self, value):
        assert(len(self.free) == len(value))
        counter=0
        for i,n in enumerate(len(self.npars)):
            self.models[i].free = value[counter:counter+n]
            counter += n

    def get_model_and_parameter_index(self,i):
        """ For a given parameter 'i', get the 
            model index (in the list of models)
            and parameter index (in the list of
            parameters for the selected model. """
        i=self.name_mapper(i)
        counter=0
        for k,n in enumerate(self.npars):
            if i < counter+n:
                break
            counter+=n
        model_index=k
        parameter_index=i-counter
        return model_index,parameter_index

    def set_parameters(self,new_vals):
        """ 
            >>> m1=PowerLaw()
            >>> m2=LogParabola()
            >>> sum_model=SumModel(m1,m2)
            >>> sum_model.set_parameters(np.asarray([1,2,3,4,5,6]))
            >>> print sum_model.get_parameters()
            [ 1.  2.  3.  4.  5.  6.]
            >>> m1.get_parameters()
            array([ 1.,  2.])
            >>> m2.get_parameters()
            array([ 3.,  4.,  5.,  6.])

        """
        counter=0
        for n,model in zip(self.npars,self.models):
            model.set_parameters(new_vals[counter:counter+n])
            counter+=n

    def setp(self, i, *args, **kwargs):
        model_index,parameter_index=self.get_model_and_parameter_index(i)
        self.models[model_index].setp(parameter_index,*args, **kwargs)

    def freeze(self,i,*args,**kwargs):
        model_index,parameter_index=self.get_model_and_parameter_index(i)
        self.models[model_index].freeze(parameter_index,*args, **kwargs)


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
        
    def external_gradient(self, e):
        return np.hstack([(1-self.ct)*self.models[0].external_gradient(e).T, 
                          self.ct*self.models[1].external_gradient(e).T])

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

    def external_gradient(self,e):
        """ Assume all models have a gradient! """
        return np.append([i.external_gradient(e) for i in self.models])

    def set_flux(self,flux,*args,**kwargs):
        """ Set the flux of the source. 

            Makes sure relative normalization of each component is unchanged
                
                >>> model = SumModel(PowerLaw(),PowerLaw())
                >>> np.allclose(model.i_flux(),PowerLaw().i_flux()*2)
                True
                >>> model.set_flux(1)
                >>> print model.i_flux()
                1.0
        """
        change = flux/self.i_flux(*args,**kwargs)
        for m in self.models: 
            m.setp(0, m.getp(0)*change)

class ProductModel(CompositeModel):
    """ Model is the product of other Models.

        Note: set_flux function just adjusts normalziation of 
              first parameter, which should be correct! 
              
        These things are Easy to create: 

            >>> m1=PowerLaw(index=2);m1.set_flux(1)
            >>> m2=LogParabola(beta=2); m1.set_flux(1)
            >>> prod_model=ProductModel(m1,m2)

        The model is just the product of the simpler modesl

            >>> energy = np.asarray([1e2,1e3,1e4])
            >>> np.allclose(prod_model(energy), m1(energy) * m2(energy))
            True
    """
    operator = '*'
    name = 'ProductModel'

    def __call__(self,e):
        return np.array([model(e) for model in self.models]).prod(axis=0)

    def external_gradient(self,e):
        """ Assume all models have a gradient! """
        return np.append([i.external_gradient(e)/i.__call__(e) for i in self.models])*self.__call__(e)

class Constant(Model):
    default_p=[1.]
    default_extra_params=dict()
    param_names=['Scale']
    default_mappers=[LogMapper]
    def __call__(self,e):
        return np.ones_like(e)*self['scale']
    
    def fast_iflux(self,emin=100,emax=1e6):
        return (emax-emin)*self['scale']

    def external_gradient(self,e):
        return  np.array([np.ones_like(e)])

class InterpConstants(Model):
    default_p=[1.]*5
    default_extra_params=dict(e_breaks=np.log10([100,300,1000,3000,3e5]))
    param_names=['Scale_Vector']
    default_mappers=[LogMapper,LogMapper,LogMapper,LogMapper,LogMapper]

    def __call__(self,e):
        interp = interp1d(self.e_breaks,self.get_all_parameters())
        return interp(np.log10(e))

    def set_flux(self,flux,**kwargs):
        raise NotImplementedError("No way to set flux for InterpConstants spectral model")

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
            >>> model.save_profile(filename, emin=1, emax=1e5)
            >>> temp.seek(0)

        Load the file into the FileFunction object

            >>> file_function = FileFunction(normalization=1, file=filename)

        And the power-law is restored, even between the saved points.

            >>> energies = np.logspace(0.122, 4.83, 42)
            >>> np.allclose(model(energies), file_function(energies))
            True

        Note, you cannot set the flux of a file function:

            >>> file_function.set_flux(1)
            >>> print file_function.i_flux()
            1.0

        Previously, pickling file_function objects failed:
            
            >>> import pickle
            >>> import os
            >>> pickle.dump(file_function, open(os.devnull,'w'))
    """
    default_p=[1]
    default_extra_params=dict(file=None)
    param_names=['Normalization']
    default_mappers=[LogMapper]

    def __init__(self,**kwargs):

        super(FileFunction,self).__init__(**kwargs)

        if self.file is None:
            raise Exception("FileFunction must be created with a file.")

        file=np.genfromtxt(os.path.expandvars(self.file),unpack=True)
        self.energy,self.flux=file[0],file[1]

        self.__make_interp__()

    def __make_interp__(self):
        self.interp = interp1d(np.log10(self.energy),np.log10(self.flux),
                bounds_error=False,fill_value=-np.inf)

    def __call__(self,e):
        return self['Normalization']*10**self.interp(np.log10(e))

    def external_gradient(self,e):
        return  np.array([np.ones_like(e)])

    def __getstate__(self):
        """ You cannot pickle an interp1d object. """
        d=copy.copy(self.__dict__)
        del d['interp']
        return d

    def __setstate__(self,state):
        """ recreate interp1d object. """
        self.__dict__ = state
        self.__make_interp__()

class SmoothDoubleBrokenPowerLaw(Model):
    """ Spectral Model Taken to be the same as:

            http://glast.stanford.edu/cgi-bin/viewcvs/Likelihood/src/SmoothDoubleBrokenPowerLaw.cxx

        but with the indices the negative of the gtlike indices for internal consistency.
            
            >>> pointlike_model = SmoothDoubleBrokenPowerLaw(Beta12=0.5, Beta23=0.25)
            >>> import pyLikelihood
            >>> _funcFactory = pyLikelihood.SourceFactory_funcFactory()
            >>> gtlike_model = _funcFactory.create('SmoothDoubleBrokenPowerLaw')
            >>> for n in pointlike_model.param_names:
            ...     gtlike_model.setParam(n,(-1 if 'index' in n.lower() else 1)*pointlike_model[n])
            >>> for n in ['Beta23', 'Beta12', 'Scale']:
            ...     gtlike_model.setParam(n,getattr(pointlike_model,n))
            >>> energies = np.logspace(1, 6, 10000)
            >>> from uw.darkmatter.spectral import DMFitFunction
            >>> np.allclose(DMFitFunction.call_pylike_spectrum(gtlike_model, energies), pointlike_model(energies), rtol=1e-20, atol=1e-20)
            True
    """
    default_p=[ 1e-11, 1.0, 1e3, 2.0, 1e4, 3.0, ]
    default_extra_params=dict(Beta12=0.05, Beta23=0.05, Scale=1e3)
    param_names=[ 'Prefactor', 'Index1', 'BreakValue12', 'Index2', 'BreakValue23', 'Index3' ]
    default_mappers=[LogMapper,LinearMapper,LogMapper,LinearMapper,LogMapper,LinearMapper]

    def __call__(self,e):
        prefactor, index1, breakvalue12, index2, breakvalue23, index3 = self.get_all_parameters()

        # convert from pointlike to gtlike definition
        index1, index2, index3 = -index1, -index2, -index3

        # This is exactly the formula from gtlike
        return prefactor*(e/self.Scale)**index1*\
            (1 + (e/breakvalue12)**((index1 - index2)/self.Beta12))**-self.Beta12*\
            (1 + (e/breakvalue23)**((index2 - index3)/self.Beta23))**-self.Beta23

    def full_name(self):
        return '%s, scale=%.0f, beta12=%.3g, beta23=%.3g'% (self.pretty_name,self.Scale,self.Beta12,self.Beta23)


class GaussianSpectrum(Model):
    """ Following the defintion in:
            http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/source_models.html#Gaussian


        >>> g = Gaussian(prefactor=1, mean=0, sigma=1)
        >>> e = np.logspace(1,10,11)
        >>> np.allclose(g(e),norm.pdf(e))
    """
    default_p = [1e-9, 7e4, 1e3]
    default_extra_params = dict()
    param_names = [ "Prefactor", "Mean", "Sigma"]
    default_mappers = [LogMapper, LogMapper, LogMapper]

    def __call__(self,e):
        prefactor, mean, sigma = self.get_all_parameters()
        return prefactor/(sigma*2.0*np.pi)*np.exp(-(e-mean)**2/(2.0*sigma**2))

def convert_exp_cutoff(model):
    """ this function is need for XML parsing. """
    if model.name != 'ExpCutoff':
        raise Exception,'Cannot process %s into PLSuperExpCutoff'%(model.name)
    nm = PLSuperExpCutoff()
    nm._p    = np.append(model._p,0)
    nm.free = np.append(model.free,False)
    nm.internal_cov_matrix[:,:] = 0
    nm.internal_cov_matrix[:-1,:-1] = model.internal_cov_matrix[:,:]
    nm.e0 = model.e0
    return nm
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
