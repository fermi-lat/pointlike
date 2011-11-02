"""
Tools for ROI analysis

$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/tools.py,v 1.9 2011/10/20 21:40:26 burnett Exp $

"""
import os
import numpy as np
from scipy import optimize, integrate
from skymaps import SkyDir
from uw.like import quadform
from uw.utilities import keyword_options
from . plotting import tsmap


def ufunc_decorator(f): # this adapts a bound function
    def new_ufunc(self, par):
        return np.array(map(lambda x: f(self,x),par)) if hasattr(par, '__iter__')  else f(self,par)
    return new_ufunc


    
class Localization(object):
    """ manage localization of a source
    Implements a minimization interface,
    see also the localize function, which uses the eliptical fitter
   
    """
    defaults = (
        ('tolerance',1e-3),
        ('verbose',False),
        ('update',False,"Update the source position after localization"),
        ('max_iteration',10,"Number of iterations"),
        #('bandfits',True,"Default use bandfits"),
        ('maxdist',1,"fail if try to move further than this"),
        ('quiet', False, 'set to suppress output'),
    )

    @keyword_options.decorate(defaults)
    def __init__(self, roistat, source_name, **kwargs):
        """ roistat : an ROIstat object
            source_name : string
                the name of a source
        * Warning: sets the roistat to evaluate only likelihood changes for the given source 
        * run its initialize() function after
          *
        """
        keyword_options.process(self, kwargs)
        self.rs = roistat
        self.source = self.rs.select_source(source_name)
        self.rs.select_source(source_name)
        self.rs.update(True)
        self.maxlike=self.rs.log_like()
        self.skydir=self.saved_skydir = self.source.skydir #saved value
        self.name = self.source.name

    def log_like(self, skydir):
        """ return log likelihood at the given position"""
        self.source.skydir =skydir
        self.rs.update(True)
        w = self.rs.log_like()
        self.source.skydir = self.saved_skydir
        return w
   
    def TSmap(self, skydir):
        """ return the TS at given position, or 
            2x the log(likelihood ratio) from the nominal position
        """
        val= 2*(self.log_like(skydir)-self.maxlike)
        return val

    def get_parameters(self):
        return np.array([self.source.skydir.ra(), self.source.skydir.dec()])
    
    def set_parameters(self, par):
        self.skydir = SkyDir(par[0],par[1])
        self.source.skydir = self.skydir
        
    def __call__(self, par):
        return -self.TSmap(SkyDir(par[0],par[1]))
    
    def reset(self):
        """ restore modifications to the ROIstat
        """
        self.source.skydir=self.skydir
        self.rs.select_source(None)
      
    def dir(self):
        return self.skydir

    def errorCircle(self):
        return 0.05 #initial guess

    def spatialLikelihood(self, sd, update=False):
        return -self.log_like(sd)
        
    def localize(self):
        """Localize a source using an elliptic approximation to the likelihood surface.

            return fit position, number of iterations, distance moved, delta TS
        """
        #roi    = self.roi
        #bandfits = self.bandfits
        verbose  = self.verbose
        update    = self.update
        tolerance= self.tolerance
        l   = quadform.Localize(self,verbose = verbose)
        ld  = l.dir

        ll0 = self.spatialLikelihood(self.skydir)

        if not self.quiet:
            fmt ='Localizing source %s, tolerance=%.1e...\n\t'+7*'%10s'
            tup = (self.name, tolerance,)+tuple('moved delta ra     dec    a     b  qual'.split())
            print fmt % tup
            print ('\t'+4*'%10.4f')% (0,0,self.skydir.ra(), self.skydir.dec())
            diff = np.degrees(l.dir.difference(self.skydir))
            print ('\t'+7*'%10.4f')% (diff,diff, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
        
        old_sigma=1.0
        for i in xrange(self.max_iteration):
            try:
                l.fit(update=True)
            except:
                #raise
                l.recenter()
                if not self.quiet: print 'trying a recenter...'
                continue
            diff = np.degrees(l.dir.difference(ld))
            delt = np.degrees(l.dir.difference(self.skydir))
            sigma = l.par[3]
            if not self.quiet: print ('\t'+7*'%10.4f')% (diff, delt, l.par[0],l.par[1],l.par[3],l.par[4], l.par[6])
            if delt>self.maxdist:
                if not self.quiet: print '\t -attempt to move beyond maxdist=%.1f' % self.maxdist
                raise Exception('localize failure: -attempt to move beyond maxdist=%.1f' % self.maxdist)
            if (diff < tolerance) and (abs(sigma-old_sigma) < tolerance):
                break
            ld = l.dir
            old_sigma=sigma

        self.qform    = l
        self.lsigma   = l.sigma
        q = l.par
        self.ellipse = dict(ra=float(q[0]), dec=float(q[1]),
                a=float(q[3]), b=float(q[4]),
                ang=float(q[5]), qual=float(q[6]),
                lsigma = l.sigma)

        ll1 = self.spatialLikelihood(l.dir)
        if not self.quiet: print 'TS change: %.2f'%(2*(ll0 - ll1))

        #roi.delta_loc_logl = (ll0 - ll1)
        # this is necessary in case the fit always fails.
        delt = np.degrees(l.dir.difference(self.skydir))
        if self.update:
            self.source.skydir = l.dir
        self.delta_ts = 2*(ll0-ll1)
        self.delt = delt
        self.niter = i

        
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
        assert g(self.maxl)*g(self.maxflux)<0, \
            'bad: peak, max %.2e, %.2e, values: %.2e,%.2e' % \
            (self.maxl,self.maxflux, self.logL(self.maxl),self.logL(self.maxflux))
        yu= optimize.brentq( g, self.maxl, self.maxflux, xtol=xtol)
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
        if 'lw' not in kwargs or 'linewidth' not in kwargs: kwargs.update(lw=2)
        axes.plot(dom, self(dom),'-', **kwargs)
        color = kwargs.get('color', 'b')
        a,b = self.errors()
        if b==0: b=self.upper_limit()
        axes.axvspan(a,b, color=color, alpha=0.25)
        axes.grid(True)
        axes.set_ylim((-0.1, 1.1))
    @staticmethod
    def test(mean=1.0, sigma = 0.1, guess=1.0):
        def fun(x):
            return -0.5*((x-mean)/sigma)**2
        t =  LogLikelihood(fun, guess=guess)
        #print t
        return t
        

def localize_all(roi, **kwargs):
    """ localize all variable local sources in the roi, make TSmaps if requested
    """
    sources = [s for s in roi.sources if s.skydir is not None and np.any(s.spectral_model.free)]
    tsmap_dir = kwargs.pop('tsmap_dir', None)
    tsfits = kwargs.pop('tsfits', False) #TODO: reimplement this to generate FITS maps
    
    initw = roi.log_like()
    
    for source in sources:
        loc =roi.localize(source.name, quiet=True, update=False, **kwargs)
        source.ellipse = loc.qform.par[0:2]+loc.qform.par[3:7] +[loc.delta_ts] if hasattr(loc,'qform') else None
        if not roi.quiet: 
            print 'Localizating %s: %d iterations, moved %.3f deg, deltaTS: %.1f' % \
                (source.name, loc.niter, loc.delt, loc.delta_ts)
            labels = 'ra dec a b ang qual'.split()
            print (len(labels)*'%10s') % tuple(labels)
            p = loc.qform.par[0:2]+loc.qform.par[3:7]
            print len(p)*'%10.4f' % tuple(p)
        if tsmap_dir is not None:
            loc = source.loc
            tsize = loc.ellipse['a']*15. if hasattr(loc,'ellipse') and loc.ellipse is not None else 1.1
            pixelsize= tsize/15.;
            tsm=tsmap.plot(loc, source.name, center=source.skydir, 
                outdir=tsmap_dir, catsig=0, size=tsize, 
                pixelsize= pixelsize, # was 14: desire to have central pixel
                # todo: fix this
                assoc=source.__dict__.get('adict', None), # either None or a dictionary
                notitle=True, #don't do title
                markersize=10,
                primary_markersize=12,
                )
        loc.reset() # restore source position
    #roi.initialize()
    curw= roi.log_like()
    if (initw-curw)>1.0:
        print 'localize_all: unexpected change in roi state after localization, from %.1f to %.1f' %(initw, curw)
        return False
    else: return True
