""" A module containing LCPrimitive and its subclasses.  They implement
components of a pulsar light curve.  Includes primitives (Gaussian,
Lorentzian), etc.  as well as more sophisticated holistic templates that
provide single-parameter (location) representations of the light curve.

$Header:
/nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/lcprimitives.py,v 1.3
2010/12/17 18:19:16 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>

"""

#=====================================================================#

import numpy as np
from scipy.special import erf,i0
from scipy.integrate import simps
from scipy.interpolate import interp1d
from scipy.stats import norm
from math import sin,cos
ROOT2PI  = (2*np.pi)**0.5
ROOT2    = 2**0.5
TWOPI    = (2*np.pi)
PI       = np.pi*1
MAXWRAPS = 15

#=====================================================================#

class LCPrimitive(object):
    """ Base class for various components of a light curve.  All "analytic"
        light curve models must inherit and must implement the three
        'virtual' functions below."""

    def __call__(self,phases):
        raise NotImplementedError('Virtual function must be implemented by child class.')

    def integrate(self,x0=0,x1=1):
        raise NotImplementedError('Virtual function must be implemented by child class.')

    def fwhm(self):
        """Return the full-width at half-maximum of the light curve model."""
        raise NotImplementedError('Virtual function must be implemented by child class.')

    def init(self):
        self.p = np.asarray([1])
        self.free = np.asarray([True])
        self.pnames = []
        self.name = 'Default'
         
    def __init__(self,**kwargs):
        self.init()
        self.__dict__.update(kwargs)
        self.p = np.asarray(self.p)
        self.free = np.asarray(self.free)
        self.errors = np.zeros_like(self.p)
        self.shift_mode = False
        self.cache = False

    def enable_cache(self,phases):
        self.free  = np.asarray([False] * len(self.p))
        self.cache_vals = self(phases)
        self.cache = True

    def disable_cache(self):
        self.cache_vals = None
        self.cache      = False

    def set_parameters(self,p): self.p[self.free] = p

    def get_parameters(self): return self.p[self.free]

    def get_location(self,error=False):
        if error: return np.asarray([self.p[-1],self.errors[-1]])
        return self.p[-1]

    def set_location(self,loc):
        self.p[-1] = loc

    def get_norm(self,error=False):
        if error: return np.asarray([self.p[0],self.errors[0]])
        return self.p[0]

    def get_width(self,error=False,fwhm=False):
        scale = self.fwhm()/self.p[1] if fwhm else 1
        if error: return np.asarray([self.p[1],self.errors[1]])*scale
        return self.p[1]*scale
    
    def get_gradient(self,phases):
        g = self.gradient(phases)
        ###N.B. -- the "-1" comes from the normalization constraint!
        g[0] -= 1 
        return g[self.free]

    def gradient(self):
        raise NotImplementedError('No gradient function found for this object.')

    def random(self,n):
        """ Default is accept/reject.""" 
        if n < 1: return 0
        M = self(self.p[-1])
        rvals    = np.empty(n)
        position = 0
        rfunc    = np.random.rand
        while True:
            cand_phases = rfunc(n)
            cand_phases = cand_phases[rfunc(n) < self(cand_phases)/M]
            ncands = len(cand_phases)
            if ncands == 0: continue
            rvals[position:position + ncands] = cand_phases[:n-position]
            position += ncands
            if position >= n: break
        return rvals
      
    def __str__(self):
        m=max([len(n) for n in self.pnames])
        l = []
        errors = self.errors if hasattr(self,'errors') else [0]*len(self.pnames)
        for i in xrange(len(self.pnames)):
            n=self.pnames[i][:m]
            t_n = n+(m-len(n))*' '
            l += [t_n + ': %.4f +\- %.4f'%(self.p[i],errors[i])]
        l = [self.name+'\n------------------'] + l
        return '\n'.join(l)

#=====================================================================#

class LCWrappedFunction(LCPrimitive):
    """ Super-class for profiles derived from wrapped functions."""

    def __call__(self,phases):
        results = self.base_func(phases)
        for i in xrange(1,MAXWRAPS+1):
            results += self.base_func(phases,index= i)
            results += self.base_func(phases,index=-i)
        return results

    def gradient(self,phases):
        """ Return the gradient evaluated at a vector of phases.

            output : a 3xlen(phases) ndarray, the 3-dim gradient at each phase
        """
        results = self.base_grad(phases)
        for i in xrange(1,MAXWRAPS+1):
            for j in [i,-i]:
                results += self.base_grad(phases,index=j)
        return results

    def integrate(self,x1=0,x2=1):
        results = self.base_int(x1,x2,index=0)
        for i in xrange(1,MAXWRAPS+1):
            for j in [i,-i]:
                results += self.base_int(x1,x2,index=j)
        return results

    def base_func(self,phases,index=0):
        raise NotImplementedError('No base_func function found for this object.')
        
    def base_grad(self,phases,index=0):
        raise NotImplementedError('No base_grad function found for this object.')
        
    def base_int(self,phases,index=0):
        raise NotImplementedError('No base_int function found for this object.')

#=====================================================================#

class LCWrappedFunction2(LCWrappedFunction):
    """ Super class for two-sided wrapped functions. """

    def init(self):
        self.p    = np.asarray([1,0.03,0.03,0.5])
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Norm','Width1','Width2','Location']
        self.two_sided = True
        self.lock_widths = False

    def fwhm(self):
        x0,left_p,right_p = self._equalize()
        self.prim.p[:] = left_p
        w_left = self.prim.fwhm()
        self.prim.p[:] = right_p
        return (w_left+self.prim.fwhm())/2
        
    def set_parameters(self,p):
        super(LCWrappedFunction2,self).set_parameters(p)
        if (not self.free[2]) and self.lock_widths:
            self.p[2] = self.p[1]

    def _equalize(self):
        """ Enforce continuity at x0.  This base implementation
            works for gaussians and lorentzians."""
        norm,w1,w2,x0 = self.p
        n1 = 2*norm/(1+w2/w1); n2 = n1*(w2/w1)
        return x0,(n1,w1,x0),(n2,w2,x0)

    def base_func(self,phases,index=0):
        x0,left_p,right_p = self._equalize()
        mask = (phases - x0 + index) < 0
        if ((index < 0) or np.all(mask)):
            self.prim.p[:] = left_p
            return self.prim.base_func(phases,index=index)
        elif ((index > 0) or np.all(~mask)):
            self.prim.p[:] = right_p
            return self.prim.base_func(phases,index=index)
        else:
            rvals = np.empty_like(phases)
            self.prim.p[:] = left_p
            rvals[mask] = self.prim.base_func(phases[mask],index=index)
            self.prim.p[:] = right_p
            rvals[~mask] = self.prim.base_func(phases[~mask],index=index)
            return rvals
            
    def base_grad(self,phases,index=0):
        raise NotImplementedError,'Have not done gradient yet...'

    def base_int(self,x1,x2,index=0):
        x0,left_p,right_p = self._equalize()
        left = x1+index; right = x2+index
        if right <= x0:
            self.prim.p[:] = left_p
            return self.prim.base_int(x1,x2,index=index)
        elif left >= x0:
            self.prim.p[:] = right_p
            return self.prim.base_int(x1,x2,index=index)
        else:
            self.prim.p[:] = left_p
            rval = self.prim.base_int(x1,x0,index=index)
            self.prim.p[:] = right_p
            return rval + self.prim.base_int(x0,x2,index=index)

#=====================================================================#

class LCGaussian(LCWrappedFunction):
    """ Represent a (wrapped) Gaussian peak.

        Parameters
        Norm      :     fraction of photons belonging to peak
        Width     :     the standard deviation parameter of the norm dist.
        Location  :     the mode of the Gaussian distribution
    """

    def init(self):
        self.p    = np.asarray([1,0.03,0.5])
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Norm','Width','Location']
        self.name = 'Gaussian'

    def fwhm(self): return self.p[1]*(8 * np.log(2))**0.5

    def base_func(self,phases,index=0):
        norm,width,x0 = self.p
        z = (phases - x0 + index)/width
        return (norm/(width*ROOT2PI))*np.exp( - 0.5*z**2 )

    def base_grad(self,phases,index=0):
        norm,width,x0 = self.p
        z      = (phases - x0 + index)/width
        fvals  = self.base_func(phases,index=index)
        return np.asarray([fvals/norm,fvals/width*(z**2 - 1.),fvals/width*z])

    def base_int(self,x1,x2,index=0):
        norm,width,x0 = self.p
        z1 = (x1 - x0 + index)/width; z2 = (x2 - x0 + index)/width;
        return 0.5*norm*(erf(z2/ROOT2)-erf(z1/ROOT2))

    def random(self,n):
        return np.mod(norm.rvs(loc=self.p[-1],scale=self.p[1],size=n),1)

#=====================================================================#

class LCLorentzian(LCWrappedFunction):
    """ Represent a Lorentzian peak.
   
        Parameters
        Norm      :     fraction of photons belonging to peak
        Width     :     the standard width paramater of Cauchy distribution, HWHM
        Location  :     the center of the peak in phase
    """

    def init(self):
        self.p = np.asarray([1,0.1,0.1])
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Norm','Width','Location']
        self.name = 'Lorentzian'

    def fwhm(self): return self.p[1]*2

    def base_func(self,phases,index=0):
        norm,gamma,x0 = self.p
        return (gamma/np.pi*norm) / ( (phases - x0 + index)**2 + gamma**2)

    def base_int(self,x1,x2,index=0):
        norm,gamma,x0 = self.p
        z1 = (x1 - x0 + index)/gamma; z2 = (x2 - x0 + index)/gamma
        return norm/np.pi * (np.arctan(z2) - np.arctan(z1))

#=====================================================================#

class LCVonMises(LCPrimitive):
    """ Represent a peak from the von Mises distribution.  This function is
        used in directional statistics and is naturally wrapped.
   
        Parameters:
            Norm      :     fraction of photons belonging to peak
            Width     :     inverse of the 'kappa' parameter in the std. def.
            Location  :     the center of the peak in phase
    """

    def init(self):
        self.p    = np.asarray([1,0.05,0.5])
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Norm','Width','Location']
        self.name = 'VonMises'

    def fwhm(self):
        return np.arccos(self.p[1]*np.log(0.5)+1)/PI

    def __call__(self,phases):
        norm,width,loc = self.p
        return norm*np.exp(np.cos(TWOPI*(phases-loc))/width)/i0(1./width)

    def integrate(self,x0=0,x1=1):
        if (x0==0) and (x1==1): return self.p[0]
                
        dom = np.linspace(0,1,min(500,max(100,5./self.p[1])))
        cod = self(dom)
        return simps(cod,x=dom)

#=====================================================================#

class LCGaussian2(LCWrappedFunction2):
    """ Represent a (wrapped) Gaussian peak with two sides.
        Parameters
        Norm      :  fraction of photons belonging to peak
        Width1    :  the standard deviation parameter of the norm dist.
        Width2    :  the standard deviation parameter of the norm dist.
        Location  :  the mode of the Gaussian distribution
    """

    def init(self):
        super(LCGaussian2,self).init()
        self.p    = np.asarray([1,0.03,0.03,0.5])
        self.name = 'Gaussian2'
        self.prim = LCGaussian()

#=====================================================================#

class LCLorentzian2(LCWrappedFunction2):
    """ Represent a (wrapped) Lorentzian peak with two sides.
        Parameters
        Norm      :  fraction of photons belonging to peak
        Width1    :  the HWHM of Cauchy distribution
        Width2    :  the HWHM of Cauchy distribution
        Location  :  the center of the peak in phase
    """

    def init(self):
        super(LCLorentzian2,self).init()
        self.p = np.asarray([1,0.05,0.05,0.5])
        self.name = 'Lorentzian2'
        self.prim = LCLorentzian()

#=====================================================================#

class LCTopHat(LCPrimitive):
   """Represent a top hat function.
   
      Parameters:
         Norm      :     fraction of photons belonging to peak
         Width     :     right edge minus left edge
         Location  :     center of top hat
   """

   def init(self):
      self.p    = np.asarray([1,0.03,0.5])
      self.free = np.asarray([True]*len(self.p))
      self.pnames = ['Norm','Width','Location']
      self.name = 'TopHat'
      self.fwhm_scale = 1

   def __call__(self,phases,wrap=True):
      norm,width,x0 = self.p
      v = norm/width
      return np.where(np.abs(phases - x0%1) < width/2,v,0)

   def integrate(self,x1=-np.inf,x2=np.inf):
      # achtung -- kluge for now
      norm,width,x0 = self.p
      return norm

#=====================================================================#

class LCHarmonic(LCPrimitive):
    """Represent a sinusoidal shape corresponding to a harmonic in a Fourier expansion.
   
      Parameters:
         Norm      :     coefficient of the term (must be positive)
         Location  :     the phase of maximum

    """

    def init(self):
        self.p    = np.asarray([1.,0.])
        self.free = np.asarray([True]*len(self.p))
        self.order= 1
        self.pnames = ['Norm','Location']
        self.name = 'Harmonic'

    def __call__(self,phases):
        norm,x0 = self.p
        return norm*np.cos( (TWOPI*self.order) * (phases - x0 ) )

    def integrate(self,x1=0,x2=1):
        if x1 == 0 and x2 == 1: return 0
        else:
            # not yet implemented
            raise Exception

#=====================================================================#

class LCEmpiricalFourier(LCPrimitive):
    """ Calculate a Fourier representation of the light curve.
        The only parameter is an overall shift.
        Cannot be used with other LCPrimitive objects!
   
        Parameters:
           Shift     :     overall shift from original template phase
    """

    def init(self):
        self.nharm = 20
        self.p     = np.asarray([0.])
        self.free  = np.asarray([True])
        self.pnames= ['Shift']
        self.name  = 'Empirical Fourier Profile'
        self.shift_mode = True

    def __init__(self,phases=None,input_file=None,**kwargs):
        """Must provide either phases or a template input file!"""
        self.init()
        self.__dict__.update(kwargs)
        if input_file is not None: self.from_file(input_file)
        if phases is not None: self.from_phases(phases)

    def from_phases(self,phases):
        n = float(len(phases))
        harmonics = np.arange(1,self.nharm+1)*(2*np.pi)
        self.alphas = np.asarray([(np.cos(k*phases)).sum() for k in harmonics])
        self.betas  = np.asarray([(np.sin(k*phases)).sum() for k in harmonics])
        self.alphas /= n; self.betas /= n;
        self.harmonics = harmonics

    def from_file(self,input_file):
        if type(input_file) == type(''):
            toks = [line.strip().split() for line in file(input_file) if len(line.strip()) > 0 and '#' not in line]
        else: toks = input_file
        alphas = []
        betas  = []
        for tok in toks:
            if len(tok) != 2: continue
            try:
                a = float(tok[0])
                b = float(tok[1])
                alphas += [a]
                betas  += [b]
            except: pass
        n = len(alphas)
        self.alphas = np.asarray(alphas)
        self.betas  = np.asarray(betas)
        self.nharm = n
        self.harmonics = np.arange(1,n+1)*(2*np.pi)

    def to_file(self,output_file):
        f = file(output_file,'w')
        f.write('# fourier\n')
        for i in xrange(self.nharm):
            f.write('%s\t%s\n'%(self.alphas[i],self.betas[i]))

    def __call__(self,phases):
        shift = self.p[0] ; harm = self.harmonics
        if shift != 0:
            """ shift theorem, for real coefficients
                It's probably a wash whether it is faster to simply 
                subtract from the phases, but it's more fun this way! """
            c = np.cos(harms * shift)
            s = np.sin(harms * shift)
            a = c*self.alphas - s*self.betas
            b = s*self.alphas + c*self.betas
        else: a,b = self.alphas,self.betas

        ak = np.asarray([np.cos(phases*k) for k in harm]).transpose()
        bk = np.asarray([np.sin(phases*k) for k in harm]).transpose()
        return (1 + 2*(a*ak + b*bk).sum(axis=1))

    def integrate(self):
        """ The Fourier expansion by definition includes the entire signal, so
        the norm is always unity."""
        return 1

#=====================================================================#

class LCKernelDensity(LCPrimitive):
    """ Calculate a kernel density estimate of the light curve.
        The bandwidth is empirical, determined from examining several pulsars.
        The only parameter is an overall shift.
        Cannot be used with other LCPrimitive objects!

        Parameters:
            Shift     :     overall shift from original template phase
    """

    def init(self):
        self.bw = None
        self.use_scale    = True
        self.max_contrast = 1
        self.resolution = 0.001 #interpolation sampling resolution
        self.p     = np.asarray([0.])
        self.free  = np.asarray([True])
        self.pnames= ['Shift']
        self.name  = 'Gaussian Kernel Density Estimate'
        self.shift_mode = True

    def __init__(self,phases=None,input_file=None,**kwargs):
        """Must provide either phases or a template input file!"""
        self.init()
        self.__dict__.update(kwargs)
        if input_file is not None: self.from_file(input_file)
        if phases is not None: self.from_phases(phases)

    def from_phases(self,phases):
        n = len(phases)
        # put in "ideal" HE bins after initial calculation of pulsed fraction
        # estimate pulsed fraction
        h  = np.histogram(phases, bins = 100)
        o  = np.sort(h[0])
        p  = float((o[o > o[15]] - o[15]).sum()) / o.sum() # based on ~30% clean offpulse
        b = o[15]
        if self.bw is None:
            self.bw = (0.5 * (p**2 * n)**-0.2)/(2*np.pi)
            print p,self.bw
            local_p = np.maximum(h[0] - b,0).astype(float) / h[0]
            print local_p, b
            bgbw = ((1-p)**2*n)**-0.2/(2*np.pi)
            print bgbw
            self.bw = np.minimum((local_p**2 * h[0])**-0.2/100.,bgbw)

        keys = np.searchsorted(h[1],phases)
        keys[keys==len(h[0])] = len(h[0]) - 1
        bw = self.bw[keys]
        print len(phases),len(bw),type(bw)


        phases = phases.copy()
        self.phases = phases
        self.phases.sort()
        phases = np.asarray(phases)
        self.phases = np.asarray(phases)
        print type(self.phases),type(phases)
        hi_mask = np.asarray(phases > 0.9)
        lo_mask = np.asarray(phases < 0.1)
        self.num = len(phases)
        self.phases = np.concatenate([phases[hi_mask]-1,phases])
        self.phases = np.concatenate([self.phases,1+phases[lo_mask]])

        print len(hi_mask),type(hi_mask),type(bw),len(bw)
        self.bw = np.concatenate([bw[hi_mask],bw])
        self.bw = np.concatenate([self.bw,bw[lo_mask]])

        #if self.bw is None:
        #   self.bw = len(phases)**-0.5

        dom = np.linspace(0,1,int(1./self.resolution))
        vals = self.__all_phases__(dom)
        ip = interp1d(dom,vals)
        mask = (self.phases > 0)&(self.phases < 1)

        """
        # this is a scaling that somehow works very well...
        vals = ip(self.phases[mask])
        scale = vals/(vals.max()-vals.min())*self.max_contrast
        #scale = scale**2
        #scale = (vals/vals.min())**1.5
        if self.use_scale:
         bw = self.bw / scale
        else:
         bw = self.bw * np.ones(len(vals))
        #bw = np.maximum(bw,self.resolution)
        """
        hi_mask = (self.phases[mask] > 0.9)
        lo_mask = (self.phases[mask] < 0.1)
        self.bw = np.concatenate([bw[hi_mask],bw])
        self.bw = np.concatenate([self.bw,bw[lo_mask]])

        vals = self.__all_phases__(dom) #with new bandwidth
        self.interpolator = interp1d(dom,vals)
        self.xvals,self.yvals = dom,vals

    def __all_phases__(self,phases):
        return np.asarray([(np.exp(-0.5 * ((ph - self.phases)/self.bw)**2 )/self.bw).sum() for ph in phases])/((2*np.pi)**0.5*self.num)

    def from_file(self,input_file):
        if type(input_file) == type(''):
            toks = [line.strip().split() for line in file(input_file) if len(line.strip()) > 0 and '#' not in line]
        else: toks = input_file

        xvals,yvals = np.asarray(toks).astype(float).transpose()
        self.xvals,self.yvals = xvals,yvals
        self.interpolator = interp1d(xvals,yvals)   

    def __call__(self,phases):
        shift = self.p[0]
        if shift == 0:  return self.interpolator(phases)
        # think this sign convention consistent with other classes - check.
        phc = np.mod(phases.copy() - shift,1)
        """ MTK changed 25 Jul 2011
        if shift >= 0 : phc[phc<0] += 1
        else: phc[phc > 1] -= 1
        """
        return self.interpolator(phc) 

    def to_file(self,output_file):
        f = file(output_file,'w')
        f.write('# kernel\n')
        for i in xrange(len(self.xvals)):
            f.write('%s\t%s\n'%(self.xvals[i],self.yvals[i]))

    def integrate(self,x0=0,x1=1):
        if (x0==0) and (x1==1): return 1.
        # crude nearest neighbor approximation
        x = self.interpolator.x; y = self.interpolator.y
        mask = (x >= x0) & (x <= x1)
        return simps(y[mask],x=x[mask])
        #return self.interpolator.y[mask].sum()/len(mask)

#=====================================================================#

def convert_primitive(p1,ptype='LCLorentzian'):
    """ Attempt to set the parameters of p2 to give a comparable primitive
        to p1."""
    p2 = eval('%s()'%ptype)
    p2.p[:] = p1.p[:]
    p2.p[1] = p1.fwhm()/p2.fwhm()*p2.p[1]
    return p2
