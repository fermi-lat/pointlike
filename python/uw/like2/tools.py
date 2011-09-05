"""
Tools for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/tools.py,v 1.2 2011/09/02 16:30:43 burnett Exp $

"""
import numpy as np
from scipy import optimize, integrate
from skymaps import SkyDir
from uw.like import Models

def ufunc_decorator(f): # this adapts a bound function
    def new_ufunc(self, par):
        return np.array(map(lambda x: f(self,x),par)) if hasattr(par, '__iter__')  else f(self,par)
    return new_ufunc


class SourceFlux(object):
    """ measure the energy-dependent flux for a given source
    An object can be passed to LogLikelihood to measure errors, limit, TS
    It can be set for all bands, just the band(s) at given energy, or a single band:
        see select_band
    But beware: it alters the source's spectral_model if so: be sure to call reset when done!
    """
    def __init__(self, rstat, source_name, quiet=True):
        """ rstat : ROIstat object
            source_name : name of one of the sources in the SourceList
        """
        self.rs = rstat
        self.rs.quiet=quiet
        self.source = rstat.sources.find_source(source_name)
        parname = source_name+'_Norm'
        try:
            self.pindex = list(rstat.parameter_names).index(parname)
        except:
            raise Exception('did not find parameter name, %s, for source flux'%parname)
        self.saved_model=self.source.spectral_model
        self.energies = np.sort(list(set([ sm.band.e for sm in rstat.all_bands])))
        #self.projector = ProjectToLinear(self.rs, [self.pindex])
        self.select_band(None)
        
    def reset(self):
        self.source.spectral_model = self.saved_model
        
    @ufunc_decorator # make this behave like a ufunc
    def __call__(self, eflux):
        """ eflux : double
                energy flux in eV units
        """
        self.source.spectral_model[0] = max(eflux,1e-3)*self.factor
        self.rs.update()
        return self.rs.log_like()
        
    def select_band(self, index, event_class=None):
        """ Select an energy band or bands
        parameters:
            index: None or integer
                an index into the list of energies; if None, select all bands
                and use the current spectral model, otherwise a powerlaw to 
                represent an model-independent flux over the band.
            class : None or integer
                if None, select both front and back, otherwise 0/1 for front/back
        Sets self.factor as conversion factor from flux to eflux in eV    
            
        """
        if index==None: #select all bands, use input model
            self.selected_energy = self.source.spectral_model.e0
            self.rs.select_bands()
            self.factor = self.saved_model[0]/(self.saved_model.eflux*1e6)
            self.source.spectral_model = self.saved_model
        else:
            # selected band(s) at a given energy: use a powerlaw
            self.selected_energy = energy =self.energies[index]
            self.source.spectral_model = Models.PowerLaw(free=[True,False],p=[1e-11,2.1], 
                        e0=self.selected_energy) 
            class_select = lambda x : True if event_class is None else x==event_class
            self.rs.select_bands(lambda b: b.e==energy and class_select(b.ec))
            assert len(self.rs.selected_bands)>0, 'did not find any bands for energy %.1f' % energy
            self.factor = 1.0/(energy**2*1e6) 
        return LogLikelihood(self)
        

class Localization(object):
    """ manage localization of a source
    Implements a minimization interface
    TODO: add call to the elliptical localization
    
    """

    def __init__(self, roistat, source_name):
        """ roistat : an ROIstat object
            source_name : string
                the name of a source
        """
        self.rs = roistat
        self.source_mask = np.array([source.name==source_name for source in roistat.sources])
        if sum(self.source_mask)!=1:
            raise Exception('Localization: source %s not found'%source_name)
        self.rs.initialize(self.source_mask)
        self.rs.update(True)
        self.maxlike=self.rs.log_like()
        self.source = np.array(roistat.sources)[self.source_mask][0]
        self.skydir = self.source.skydir
   
    def log_like(self, skydir):
        """ return log likelihood at the given position"""
        self.source.skydir =skydir
        self.rs.update(True)
        return self.rs.log_like()
   
    def TS(self, skydir):
        """ return the TS at given position, or 
            2x the log(likelihood ratio) from the nominal position
        """
        return 2*(self.log_like(skydir)-self.maxlike)

    def get_parameters(self):
        return np.array([self.source.skydir.ra(), self.source.skydir.dec()])
    
    def set_parameters(self, par):
        self.skydir = SkyDir(par[0],par[1])
        self.source.skydir = self.skydir
        
    def __call__(self, par):
        return -self.TS(SkyDir(par[0],par[1]))
    
    def reset(self):
        """ restore modifications to the ROIstat
        """
        self.source.skydir=self.skydir
        self.rs.update(True)
        self.rs.initialize()
        
        
class LogLikelihood(object):
    """ manage a 1-dimensional likelihood function """

    def __init__(self, loglikelihood,  guess=1.0):
        """ 
        loglikelihood: a function of one parameter: expect to have a maximum
        guess: a guess for setting the scale
        """
        self.loglike = loglikelihood
        
        self.wpeak=0.
        self.maxl = self.maximum(guess)
        self.wpeak = self.logL(self.maxl)
        self.maxflux = self.find_maxflux(guess)
        self.tot=None

    def __str__(self):
        h = tuple('max -sig +sig 95% TS'.split())
        t =  (self.maxl,) + self.errors() +(self.upper_limit(),  self.TS())
        n = len(t)
        return '\n'.join([(n*'  %-9s') % h, ((n-1)*' %10.2e'+'%10.1f') % t])
        
    def logL(self, norm):
        """ evaluate the log likelihood function
        norm: a value, or array of values
        """
        if hasattr(norm, '__iter__'):
            return np.array( map(self.loglike,norm) )
        return self.loglike(norm)
        
    def __call__(self, x):
        """ evalate the likelihood function, normalized to 1 at peak"""
        return np.exp(self.logL(x)-self.wpeak)

    def maximum(self, guess=1.0, disp=0):
        """ find the position of the maximum likelihood
        val : starting value for fmin
        disp: pass to fmin
        
        """
        ret = optimize.fmin( lambda x: -self.logL(x), guess, disp=disp)[0]
        return ret if ret>0 else  0
        
    def find_maxflux(self, guess=1.0, tol=10):
        """ find an appropriate maximum flux for evaluating the likelihood
        """
        wmax = self.logL(self.maxl)
        v = self.maxl*2
        if v>1e-15:
            for i in range(10):
                if self.logL(v) < wmax-tol: break
                v*=np.e
            return v
        # case with just a limit: return flux where log likelihood is -5
            
        z = -self.logL(0)+5
        g = lambda x: self.logL(x)+z
        a = 0.1*guess
        for i in range(10):
            if g(a)>0 : break
            a /= 10
        b =10*guess
        for i in range(10):
            if g(b)<0: break
            b *= 10
        v =  optimize.brentq(g, a, b , xtol=1e-3*guess )
        return v
            
    def errors(self, tol=None):
        """ tuple with lower, upper 1-sigma uncertainty flux values (profile method)"""
        delta = self.wpeak -0.5
        g = lambda x: self.logL(x)-delta # function that should be zero at the 1 sigma points
        xtol = 1e-3*self.maxflux if tol is None else tol
        yl=yu=0
        if self.maxl==0: return (yl,yu)
        try:
            yl= optimize.brentq( g, 0,         self.maxl,    xtol=xtol)
        except: pass
        try:
            assert g(self.maxl)*g(self.maxflux)<0, \
                'bad: peak, max %.2e, %.2e,\n values: %.2e,%.2e' % \
                (self.maxl,self.maxflux, self.logL(self.maxl),self.logL(self.maxflux))
            yu= optimize.brentq( g, self.maxl, self.maxflux, xtol=xtol)
        except: 
            raise
        return (yl,yu) 
            
    def integral(self, x):
        """ the integral of the normalized function"""
        if self.tot is None:
            self.tot=1.0
            self.tot=self.integral(self.maxflux)
        f =  lambda x : integrate.quad(self, 0, x, epsabs=1e-4,epsrel=1e-4)[0]/self.tot
        if hasattr(x, '__iter__'):
            return np.array( map(f,x) )
        return f(x)
        
    
    def upper_limit(self, cl = 0.95):
        """ the flux value at the confidence level"""
        if cl>1: cl /= 100.
        t =0
        try:
            t = optimize.brentq( lambda x: self.integral(x) - cl, 0, self.maxflux, xtol=1e-3*self.maxflux)
        except: pass
        return t
     
    def TS(self):
        return 2*( self.logL(self.maxl) - self.logL(0) )
        
    def plot(self,fignum=10, axes = None, **kwargs):
        """ a simple plot of the likelihood """
        import pylab as plt
        if axes is None:
            plt.close(fignum)
            plt.figure(fignum)
            axes = plt.gca()
        dom = np.linspace(0, self.maxflux, 101)
        axes.plot(dom, self(dom),'-', **kwargs)
        color = kwargs.get('color', 'b')
        a,b = self.errors()
        if b==0: b=self.upper_limit()
        axes.axvspan(a,b, color=color, alpha=0.5)
        axes.grid(True)
        axes.set_ylim((-0.1, 1.1))
    @staticmethod
    def test(mean=1.0, sigma = 0.1, guess=1.0):
        def fun(x):
            return -0.5*((x-mean)/sigma)**2
        t =  LogLikelihood(fun, guess=guess)
        #print t
        return t
        
