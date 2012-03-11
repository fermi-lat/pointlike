""" A module containing LCPrimitive and its subclasses.  They implement
components of a pulsar light curve.  Includes primitives (Gaussian,
Lorentzian), etc.  as well as more sophisticated holistic templates that
provide single-parameter (location) representations of the light curve.

$Header:
/nfs/slac/g/glast/ground/cvs/pointlike/python/uw/pulsar/lcprimitives.py,v 1.3
2010/12/17 18:19:16 kerrm Exp $

author: M. Kerr <matthew.kerr@gmail.com>

"""

import numpy as np
from scipy.special import erf,i0
from scipy.integrate import simps,quad
from scipy.interpolate import interp1d
from scipy.stats import norm,cauchy
from math import sin,cos,sinh,cosh,atan,tan
ROOT2PI  = (2*np.pi)**0.5
R2DI     = (2/np.pi)**0.5
ROOT2    = 2**0.5
TWOPI    = (2*np.pi)
PI       = np.pi*1
MAXWRAPS = 15
MINWRAPS = 3
WRAPEPS  = 1e-8

# TODO -- possible "LCBase" class with certain method common to LCPrimitive and LCTemplate

def two_comp_mc(n,w1,w2,loc,func):
    """ Helper function to generate MC photons from a two-peaked distribution.
        n -- total number of photons
        w1 -- scale parameter for func, lefthand peak
        w2 -- scale parameter for func, righthand peak
        loc -- position of peak
        func -- an 'rvs' function from scipy
    """
    frac1 = w1/(w1+w2)
    # number of photons required from left side
    n1 = (np.random.rand(n) < frac1).sum()
    r1 = func(loc=0,scale=w1,size=n1)
    # reflect and relocate photons to right or lef side
    r1 = loc + np.where(r1<=0,r1,-r1)
    r2 = func(loc=0,scale=w2,size=n-n1)
    r2 = loc + np.where(r2>0,r2,-r2)
    return np.mod(np.append(r1,r2),1)

class LCPrimitive(object):
    """ Base class for various components of a light curve.  All "analytic"
        light curve models must inherit and must implement the three
        'virtual' functions below."""

    def __call__(self,phases):
        raise NotImplementedError('Virtual function must be implemented by child class.')

    def integrate(self,x0=0,x1=1):
        """ Base implemention with scipy quad."""
        return quad(self,x0,x1)[0]

    def fwhm(self):
        """Return the full-width at half-maximum of the light curve model."""
        return self.hwhm(0)+self.hwhm(1)

    def hwhm(self,right=False):
        """Return the half-width at half-maximum of the light curve model."""
        raise NotImplementedError('Virtual function must be implemented by child class.')

    def init(self):
        self.p = np.asarray([1])
        self.free = np.asarray([True])
        self.pnames = []
        self.name = 'Default'
         
    def __init__(self,**kwargs):
        self.init()
        if not hasattr(self,'bounds'):
            self.bounds = [[None,None]*len(self.p)] # default
        self.__dict__.update(kwargs)
        self.p = np.asarray(self.p)
        self.free = np.asarray(self.free)
        self.bounds = np.asarray(self.bounds)
        self.errors = np.zeros_like(self.p)
        self.shift_mode = False

    def set_parameters(self,p): self.p[self.free] = p

    def get_parameters(self): return self.p[self.free]

    def get_bounds(self): return self.bounds[self.free]

    def get_location(self,error=False):
        if error: return np.asarray([self.p[-1],self.errors[-1]])
        return self.p[-1]

    def set_location(self,loc):
        self.p[-1] = loc

    def get_norm(self,error=False):
        #if error: return np.asarray([self.p[0],self.errors[0]])
        #return self.p[0]
        return 1

    def get_width(self,error=False,hwhm=False,right=False):
        """ Return the width of the distribution.
            Keyword arguments:
            -----------------
            error   [False] if True, return tuple with value and error
            hwhm    [False] if True, scale width to be HWHM
            right   [False] if True, return "right" component, else "left".
                            There is no distinction for symmetric dists.
        """
        scale = self.hwhm()/self.p[right] if hwhm else 1
        if error: return np.asarray([self.p[right],self.errors[right]])*scale
        return self.p[1+right]*scale
    
    def get_gradient(self,phases):
        g = self.gradient(phases)
        ###N.B. -- the "-1" comes from the normalization constraint!
        #g[0] -= 1 
        return g[self.free]

    def gradient(self):
        raise NotImplementedError('No gradient function found for this object.')

    def random(self,n):
        """ Default is accept/reject.""" 
        if n < 1: return 0
        M = self(np.asarray([self.p[-1]])) # peak amplitude
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
            fstring = '' if self.free[i] else ' [FIXED]'
            n=self.pnames[i][:m]
            t_n = n+(m-len(n))*' '
            l += [t_n + ': %.4f +\- %.4f%s'%(self.p[i],errors[i],fstring)]
        l = [self.name+'\n------------------'] + l
        return '\n'.join(l)

    def approx_gradient(self,phases,eps=1e-5):
        orig_p = self.get_parameters().copy()
        g = np.zeros([len(orig_p),len(phases)])
        weights = np.asarray([-1,8,-8,1])/(12*eps)

        def do_step(which,eps):
            p0 = orig_p.copy()
            p0[which] += eps
            self.set_parameters(p0)
            return self(phases)

        for i in xrange(len(orig_p)):
            # use a 4th-order central difference scheme
            for j,w in zip([2,1,-1,-2],weights):
                g[i,:] += w*do_step(i,j*eps)

        self.set_parameters(orig_p)
        return g

    def check_gradient(self,tol=1e-6,eps=1e-5,quiet=False):
        """ Test gradient function with a set of MC photons."""
        ph = self.random(1000)
        g1 = self.gradient(ph)
        g2 = self.approx_gradient(ph,eps=eps)
        anyfail = False
        for i in xrange(g1.shape[0]):
            d1 = np.abs(g1[i]-g2[i])
            d2 = d1/g1
            fail = np.any((d1>tol) | (d2>tol))
            if not quiet:
                pass_string = 'FAILED' if fail else 'passed'
                print '%d (%s) %.3g (abs) %.3g (frac)'%(i,pass_string,d1.max(),d2.max())
            anyfail = anyfail or fail
        return not anyfail

    def sanity_checks(self,eps=1e-6):
        """ A few checks on normalization, integration, etc. """
        from scipy.integrate import quad
        errfac = 1
        # Normalization test
        y,ye = quad(self,0,1)
        #t1 = abs(self.p[0]-y)<(ye*errfac)
        t1 = abs(1-y)<(ye*errfac)
        # integrate method test
        #t2 = abs(self.p[0]-self.integrate(0,1))<eps
        t2 = abs(1-self.integrate(0,1))<eps
        # FWHM test
        t3 = (self(self.p[-1])*0.5-self(self.p[-1]-self.fwhm()/2))<eps
        # gradient test
        try:
            t4 = self.check_gradient(quiet=True)
        except: t4 = False
        # boundary conditions
        t5 = abs(self(0)-self(1-eps))<eps
        if not t1: print 'Failed Normalization test'
        if not t2: print 'Failed integrate method test'
        if not t3: print 'Failed FWHM test'
        if not t4: print 'Failed gradient test'
        if not t5: print 'Did not pass boundary conditions'
        return np.all([t1,t2,t3,t4,t5])

class LCWrappedFunction(LCPrimitive):
    """ Super-class for profiles derived from wrapped functions.
        
        While some distributions (e.g. the wrapped normal) converge
        quickly, others (e.g. the wrapped Lorentzian) converge very slowly
        and must be truncated before machine precision is reached.

        In order to preserve normalization, the pdf is slightly adjusted:
        f(phi) = sum_(i,-N,N,g(phi+i)) + (1 - int(phi,-N,N,g(phi)) ).

        This introduces an additional parameteric dependence which must
        be accounted for by computation of the gradient.
    """

    def _norm(self,nwraps):
        """ Compute the truncated portion of the template."""
        #return self.p[0]-self.base_int(-nwraps,nwraps+1)
        return 1-self.base_int(-nwraps,nwraps+1)

    def _grad_norm(self,nwraps):
        """ Compute the gradient terms due to truncated portion.
            Default implementation is to ignore these terms, applicable
            for rapidly-converging distributions (e.g. wrapped normal with
            small width parameter)."""
        return None

    def __call__(self,phases):
        """ Return wrapped template + DC component corresponding to truncation."""
        results = self.base_func(phases)
        for i in xrange(1,MAXWRAPS+1):
            t = self.base_func(phases,index= i)
            t += self.base_func(phases,index=-i)
            results += t
            if (i>=MINWRAPS) and (np.all(t < WRAPEPS)): break
        #return results
        #print self._norm(i)
        return results+self._norm(i)

    def gradient(self,phases):
        """ Return the gradient evaluated at a vector of phases.

            output : a num_parameter x len(phases) ndarray, 
                     the num_parameter-dim gradient at each phase
        """
        results = self.base_grad(phases)
        for i in xrange(1,MAXWRAPS+1):
            t = self.base_grad(phases,index=i)
            t += self.base_grad(phases,index=-i)
            results += t
            if (i >= MINWRAPS) and (np.all(t < WRAPEPS)): break
        gn = self._grad_norm(i) 
        if gn is not None:
            for i in xrange(len(gn)):
                results[i,:] += gn[i]
        return results

    def integrate(self,x1=0,x2=1):
        #if(x1==0) and (x2==0): return self.p[0] # this is true by definition, now
        if(x1==0) and (x2==0): return 1. # this is true by definition, now
        # NB -- this method is probably overkill now.
        results = self.base_int(x1,x2,index=0)
        for i in xrange(1,MAXWRAPS+1):
            t = self.base_int(x1,x2,index=i)
            t += self.base_int(x1,x2,index=-i)
            results += t
            if np.all(t < WRAPEPS): 
                break
        return results+(x2-x1)*self._norm(i)

    def base_func(self,phases,index=0):
        raise NotImplementedError('No base_func function found for this object.')
        
    def base_grad(self,phases,index=0):
        raise NotImplementedError('No base_grad function found for this object.')
        
    def base_int(self,phases,index=0):
        raise NotImplementedError('No base_int function found for this object.')

class LCGaussian(LCWrappedFunction):
    """ Represent a (wrapped) Gaussian peak.

        Parameters
        Norm      :     fraction of photons belonging to peak
        Width     :     the standard deviation parameter of the norm dist.
        Location  :     the mode of the Gaussian distribution
    """

    def init(self):
        #self.p    = np.asarray([1,0.03,0.5])
        self.p    = np.asarray([0.03,0.5])
        #self.bounds = [ [0.001,1], [0.005,0.5], [0,1] ]
        self.bounds = [ [0.005,0.5], [0,1] ]
        self.free = np.asarray([True]*len(self.p))
        #self.pnames = ['Norm','Width','Location']
        self.pnames = ['Width','Location']
        self.name = 'Gaussian'

    def hwhm(self,right=False):
        #return self.p[1]*(2 * np.log(2))**0.5
        return self.p[0]*(2 * np.log(2))**0.5

    def base_func(self,phases,index=0):
        #norm,width,x0 = self.p
        width,x0 = self.p
        z = (phases + index - x0)/width
        return (1./(width*ROOT2PI))*np.exp(-0.5*z**2 )

    def base_grad(self,phases,index=0):
        #norm,width,x0 = self.p
        width,x0 = self.p
        z = (phases + index - x0)/width
        f = (1./(width*ROOT2PI))*np.exp(-0.5*z**2 )
        return np.asarray([f/width*(z**2 - 1.),f/width*z])

    def base_int(self,x1,x2,index=0):
        #norm,width,x0 = self.p
        width,x0 = self.p
        z1 = (x1 + index - x0)/width
        z2 = (x2 + index - x0)/width
        return 0.5*1.*(erf(z2/ROOT2)-erf(z1/ROOT2))

    def random(self,n):
        #return np.mod(norm.rvs(loc=self.p[-1],scale=self.p[1],size=n),1)
        return np.mod(norm.rvs(loc=self.p[-1],scale=self.p[0],size=n),1)

class LCGaussian2(LCWrappedFunction):
    """ Represent a (wrapped) two-sided Gaussian peak.
        Parameters
        Norm      :  fraction of photons belonging to peak
        Width1    :  the standard deviation parameter of the norm dist.
        Width2    :  the standard deviation parameter of the norm dist.
        Location  :  the mode of the distribution
    """

    def init(self):
        self.p    = np.asarray([0.03,0.03,0.5])
        self.bounds = [ [0.005,0.5], [0.005,0.5], [0,1] ]
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Width1','Width2','Location']
        self.name = 'Gaussian2'

    def hwhm(self,right=False):
        return (self.p[right])*(2 * np.log(2))**0.5

    def base_func(self,phases,index=0):
        width1,width2,x0 = self.p
        z = (phases + (index - x0))
        z *= np.where(z <= 0, 1./width1, 1./width2)
        return (R2DI/(width1+width2)) * np.exp(-0.5*z**2 )

    def base_grad(self,phases,index=0):
        width1,width2,x0 = self.p
        z = (phases + (index - x0))
        m = (z <= 0)
        w = np.where(m, width1, width2)
        z /= w
        f = (R2DI/(width1+width2)) * np.exp(-0.5*z**2 )
        k = 1./(width1+width2)
        z2w = z**2/w
        t = f*(z2w-k)
        g1 = f*(z2w*( m)-k)
        g2 = f*(z2w*(~m)-k)
        g3 = f*z/w
        return np.asarray([g1,g2,g3])

    def base_int(self,x1,x2,index=0):
        width1,width2,x0 = self.p
        if index==0 and (x1 < x0) and (x2 > x0):
            z1 = (x1 + index - x0)/width1
            z2 = (x2 + index - x0)/width2
            k1 = 2*width1/(width1+width2)
            k2 = 2*width2/(width1+width2)
            return 0.5*(k2*erf(z2/ROOT2)-k1*erf(z1/ROOT2))
        w = width1 if ((x1+index) < x0) else width2
        z1 = (x1 + index - x0)/w
        z2 = (x2 + index - x0)/w
        k = 2*w/(width1+width2)
        return 0.5*k*(erf(z2/ROOT2)-erf(z1/ROOT2))

    def random(self,n):
        """ Use multinomial technique to return random photons from
            both components."""
        return two_comp_mc(n,self.p[0],self.p[1],self.p[-1],norm.rvs)

class LCLorentzian(LCPrimitive):
    """ Represent a (wrapped) Lorentzian peak.
   
        Parameters
        Width     :     the width paramater of the wrapped Cauchy distribution,
                        namely HWHM*2PI for narrow distributions
        Location  :     the center of the peak in phase
    """
    def init(self):
        self.p = np.asarray([0.1,0.5])
        self.bounds = [ [0.005,0.5], [0,1] ]
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Width','Location']
        self.name = 'Lorentzian'

    def hwhm(self,right=False):
        # NB -- bounds on p[1] set such that this is well-defined
        return np.arccos( 2-cosh(self.p[1]) )/TWOPI

    def __call__(self,phases):
        gamma,loc = self.p
        z = TWOPI*(phases-loc)
        return sinh(gamma)/(cosh(gamma)-np.cos(z))

    def gradient(self,phases):
        gamma,loc = self.p
        z = TWOPI*(phases-loc)
        s1 = sinh(gamma); c1 = cosh(gamma)
        c = np.cos(z); s = np.sin(z)
        f = s1/(c1-c)
        f2 = f**2
        g1 = f*(c1/s1) - f2
        g2 = f2*(TWOPI/s1)*s
        return np.asarray([g1,g2])

    def random(self,n):
        return np.mod(cauchy.rvs(loc=self.p[-1],scale=self.p[0]/TWOPI,size=n),1)

    def integrate(self,x0=0,x1=1):
        gamma,loc = self.p
        if (x0==0) and (x1==1): return self.p[0]
        x0 = PI*(x0-loc)
        x1 = PI*(x1-loc)
        t = cosh(gamma/2)/sinh(gamma/2)
        return (atan(t*tan(x1))-atan(t*tan(x0)))/PI

class LCLorentzian2(LCWrappedFunction):
    """ Represent a (wrapped) two-sided Lorentzian peak.
        Parameters
        Width1    :  the HWHM of the distribution (left)
        Width2    :  the HWHM of the distribution (right)
        Location  :  the mode of the distribution
    """

    def init(self):
        self.p    = np.asarray([0.03,0.03,0.5])
        self.bounds = [ [0.005,0.5], [0.005,0.5], [0,1] ]
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Width1','Width2','Location']
        self.name = 'Lorentzian2'

    def hwhm(self,right=False):
        return self.p[right]

    def _grad_norm(self,nwraps):
        gamma1,gamma2,x0 = self.p
        z1 = (-nwraps-x0)/gamma1
        z2 = (nwraps+1-x0)/gamma2
        t = gamma2*atan(z2)-gamma1*atan(z1)
        t1 = 1./(1+z1**2)
        t2 = 1./(1+z2**2)
        k = 2/(gamma1+gamma2)/PI
        f = k*t
        g1 = -1./(gamma1+gamma2)-(atan(z1)-z1*t1)/t
        g2 = -1./(gamma1+gamma2)+(atan(z2)-z2*t2)/t
        g3 = (t1-t2)/t
        return [-f*g1,-f*g2,-f*g3]

    def base_func(self,phases,index=0):
        gamma1,gamma2,x0 = self.p
        z = (phases + (index - x0))
        z *= np.where(z<=0, 1./gamma1, 1./gamma2)
        k = 2/(gamma1+gamma2)/PI
        return k/(1+z**2)

    def base_grad(self,phases,index=0):
        gamma1,gamma2,x0 = self.p
        z = (phases + (index - x0))
        m = z < 0
        g = np.where(m,1./gamma1,1./gamma2)
        t1 = 1+(z*g)**2
        t2 = 2*(z*g)/t1
        g1 = -1/(gamma1+gamma2)+t2*((m*z)/gamma1**2)
        g2 = -1/(gamma1+gamma2)+t2*((~m*z)/gamma2**2)
        g3 = t2*g
        f = (2./(gamma1+gamma2)/PI)/t1
        return np.asarray([f*g1,f*g2,f*g3])

    def base_int(self,x1,x2,index=0):
        gamma1,gamma2,x0 = self.p
        if index==0 and (x1 < x0) and (x2 > x0):
            g1,g2 = gamma1,gamma2
        else:
            g1,g2 = [gamma1]*2 if ((x1+index) < x0) else [gamma2]*2
        z1 = (x1 + index - x0)/g1
        z2 = (x2 + index - x0)/g2
        k = (2./(gamma1+gamma2)/PI)
        return k*(g2*atan(z2)-g1*atan(z1))

    def random(self,n):
        """ Use multinomial technique to return random photons from
            both components."""
        return two_comp_mc(n,self.p[0],self.p[1],self.p[-1],cauchy.rvs)

class LCVonMises(LCPrimitive):
    """ Represent a peak from the von Mises distribution.  This function is
        used in directional statistics and is naturally wrapped.
   
        Parameters:
            Norm      :     fraction of photons belonging to peak
            Width     :     inverse of the 'kappa' parameter in the std. def.
            Location  :     the center of the peak in phase
    """

    def init(self):
        self.p    = np.asarray([0.05,0.5])
        self.free = np.asarray([True]*len(self.p))
        self.pnames = ['Width','Location']
        self.name = 'VonMises'

    def hwhm(self,right=False):
        return 0.5*np.arccos(self.p[0]*np.log(0.5)+1)/TWOPI

    def __call__(self,phases):
        width,loc = self.p
        z = TWOPI*(phases-loc)
        return np.exp(np.cos(z)/width)/i0(1./width)

class LCTopHat(LCPrimitive):
    """ Represent a top hat function.
   
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

def convert_primitive(p1,ptype=LCLorentzian):
    """ Attempt to set the parameters of p2 to give a comparable primitive
        to p1."""
    p2 = ptype()
    p2.p[:2] = p1.p[:2]
    p2.p[-1] = p1.p[-1]
    width_scale = p1.hwhm()/p2.hwhm()
    p2.p[1] = p1.p[1]*width_scale
    if len(p2.p) > 3:
        p = p1.p[1] if (len(p1.p) == 3) else p1.p[2]
        p2.p[2] = p*width_scale
    return p2

