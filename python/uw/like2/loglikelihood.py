"""Tools for parameterizing log likelihood curves.

Author(s): Eric Wallace, Matthew Kerr, Toby Burnett
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/loglikelihood.py,v 1.9 2012/11/03 22:42:39 burnett Exp $
"""
__version__ = "$Revision: 1.9 $"

import numpy as np
from scipy import optimize, special, polyfit, stats

class PoissonLogLikelihood(object):
    """A class representing a parametrized log likelihood curve.

    The parametrization used is equivalent to that described in William
    Tompkins' thesis (arxiv: astro-ph/0202141) and Nolan, et al., 2003,
    ApJ 597:615:627. The functional form is that of a Poisson distribution
    with three parameters: logL(s) = e(s_p+b)*ln(e(s+b))-e(s+b). Here s is
    the source flux of interest, s_p is the value of s that maximizes the
    likelihood, and b and e represent an effective background flux and
    exposure, respectively.
    """
    
    def __init__(self,roi,source,profile = False,atol=.5):
        """Parameters:
            roi: an ROI_user instance
            source: specification of the source of interest. This can be
                    anything accepted by roi.get_source()
            profile: if True, evaluate the profile likelihood
            atol:   absolute tolerance for agreement
        """
        self.roi = roi
        self.source = source
        self.model = self.roi.get_source(self.source).spectral_model
        self.profile = profile
        self.atol = atol
        self.cache = dict() 
        self._poiss = None
        self.maximum = self.find_max(False)
        seeds = [1/self.maximum] if self.maximum>0 else [] 
        seeds += [1e13,3e10,1e9]
        for seed in seeds:
            for bg_fac in [10.,1.,0.1,100]:
                try:
                    self._poiss = Poisson(self._do_fit(exp_seed=seed,bg_fac=bg_fac))
                except Exception:
                    self.ok = False
                    continue
                self.ok = self._check_agreement()
                if self.ok: break
            if self.ok: break
        if not self.ok:
            print 'WARNING fit and interpolation are significantly different (%.2f)'%(self.agreement)
    
    def log_like(self,s):
        """Return the log likelihood at s.
        
            Here, s is taken to be the (external, i.e., not log-transformed) 
            normalization parameter of the spectral model for the source of 
            interest.
        """
        try:
            if hasattr(s,'__len__') and len(s)==1:
                s = s[0]
            return self.cache[s]
        except KeyError:
            save_pars = self.roi.get_parameters() #These are logs, for some reason.
            self.model[0] = max(s,1e-20) #No negative fluxes
            self.roi.update()
            if (not self.profile or sum([np.any(source.spectral_model.free)
                                         for source in self.roi.sources])<2):
                ll = self.roi.log_like()
            else:
                self.model.freeze(0)
                self.roi.fit()
                ll = self.roi.log_like()
                self.model.thaw(0)
            self.roi.set_parameters(save_pars)
            self.roi.update()
            self.cache[s]=ll
            return ll

    def __call__(self,x,use_fit=True):
        """Return the log likelihood at the specified point(s).
        
        If use_fit==True, use the Poisson fit, otherwise use the likelihood
        from the roi object.
        """
        if use_fit:
            return self._poiss(x)
        else:
            if hasattr(x,'__len__'):
                return np.array([self.log_like(s) for s in x])
            else:
                return self.log_like(x)

    def find_max(self,use_fit=True):
        """Return the flux value that maximizes the likelihood.

        If use_fit is True, return the appropriate parameter from the Poisson
        fit. Otherwise, use fmin to find the maximum from the "true" likelihood.
        """
        if use_fit:
            if self._poiss is None:
                raise Exception('No fit has been performed yet, try use_fit=False.')
            return self._poiss.p[0]
            return max(0,self._poiss.p[0])
        else:
            return optimize.fmin(lambda s: -self.log_like(s),self.model[0])[0]
            return max(0,optimize.fmin(lambda s: -self.log_like(s),self.model[0])[0])
    
    def _find_domain(self,delta_logl=2,use_fit = False):
        """Find a reasonable domain on which to fit the likelihood.

        We want a region surrounding the maximum (or from zero to some
        positive value, if the maximum is <= 0) which is sufficiently broad
        to produce an accurate representation of the likelihood over a
        reasonably broad range of fluxes, but narrow enough to obviate the
        need for an excessively large number of likelihood calls. In practice,
        this is done by finding the point(s) that give a difference in the log
        likelihood with respect to the maximum of 2 (this value can be adjusted
        with the delta_logl kwarg). The evaluations used to find
        this range are cached and used as a starting set of points for the fit.
        """
        smax = self.maximum
        #ll_max = self.log_like(smax)
        #ll_zero = self.log_like(0)
        #func = lambda s: ll_max-self.log_like(s)-delta_logl
        ll_max = self(smax,use_fit)
        ll_zero = self(0,use_fit)
        func = lambda s: ll_max-self(s,use_fit)-delta_logl
        if ll_max-ll_zero<delta_logl:
            s_low = 0
        else:
            #s_low = optimize.bisect(func,0,smax,xtol=.01*smax)
            s_low = optimize.brentq(func,0,smax,xtol=1e-15)
        if smax>0:
            s_high = smax*10
        else:
            s_high = 1e-15
        while func(s_high)<0: s_high*=2
        #s_high = optimize.bisect(func,smax,s_high,xtol=.01*smax)
        s_high = optimize.brentq(func,smax,s_high,xtol=1e-15)
        if not np.all(np.isreal([s_low,s_high])):
            print('Could not find two roots!')
            return None
        return (s_low,s_high)
    
    def _do_fit(self,exp_seed=3e11,bg_fac=1):
        """Do the fit"""
        #make sure we've cached the max and desired endpoints
        low,high = self._find_domain(10)
        dom = np.asarray(sorted(self.cache.keys()))
        dom = dom[np.logical_and(dom>=low,dom<=high)]
        smax = self.maximum
        cod = self(dom,False)-self(smax,False)
        def fitfunc(p):
            self._poiss = Poisson(p)
            return self._poiss(dom)-cod
        return optimize.leastsq(fitfunc,
                                [smax,exp_seed,smax*bg_fac],
                                maxfev=10000)[0]

    def _check_agreement(self):
        """Check the agreement between the fit and the true likelihood"""
        dom = np.asarray(sorted(self.cache.keys()))
        cod = self(dom,False)-self(self.maximum,False)
        mask = np.where(np.logical_and(cod>=-10,dom>=0))
        diffs = self(dom[mask])-self(dom[mask],use_fit=False)+self(self.maximum,False)
        self.agreement = max_diff = np.abs(diffs).max()
        return max_diff<self.atol

class Poisson(object):
    """log of the three-parameter Poisson-like function use to represent the flux likelihood
    parameters are, in order:
        sp : flux at peak, if positive; if negative, there is only a limit
        e  : normalization factor to convert flux to equivalent counts. must be >0
        b  : background flux: must be >=0
        
    log likelihood expression is for flux s>=0
       w(s; sp,e,b) = e*(sp+b) * log( e*(s+b) ) - e*s + const
    the const is such that w=0 at the peak.
    
    A slightly more elegant expression, used in the cdf function, is to define
       beta = e*b
       mu = e*sp + beta
       x = e*s
    Then the log likelihood is
       w(x) = mu * log( x + beta) - x + const
    
    where the peak is at x=xp=max(0, mu-beta), and the constant is defined so that w(xp)=0 
    
    >>> w = Poisson([10.0, 1.0, 5.0])
    >>> np.array(map(w, np.linspace(0,20,5))).round(2)
    array([-6.48, -1.08,  0.  , -0.68, -2.34])
    >>> w.find_delta()
    (6.452927189667698, 14.213246175302684)
    >>> map(w.cdf, np.linspace(0,20,5))
    [0.0, 0.048674754021320932, 0.43187121882311852, 0.84347606391862051, 0.97770544018425898]
    >>> round(w.percentile(),2)
    18.1
     
    >>> w = Poisson([-5, 1.0, 10.0]) 
    >>> map(w, np.linspace(0, 10, 5))
    [0.0, -1.3842822434289523, -2.9726744594591796, -4.7019210603228885, -6.5342640972002766]
    >>> map(w.cdf, np.linspace(0, 10, 5))
    [0.0, 0.77904655517623111, 0.95837535584402034, 0.9930192704933527, 0.99892810898968387]
    >>> round(w.percentile(),2)
    4.73

    
    """
    def __init__(self,p):
        self.p = p
    
    def __call__(self,dom):
        """Return the value of the fit function for the given domain."""
        sp,e,b = self.p
        b = abs(b) # this is a bit of a kluge to keep log argument positive
        if b==0: b = 1e-20 #another kludge
        e = abs(e) #and another
        r = e*(dom+b)
        r_peak = e*(sp+b)
        if sp > 0:
            const = r_peak*np.log(r_peak) - r_peak
        else:
            #sp=0
            #const = r_peak*np.log(r_peak) - r_peak
            t = e*b
            const = r_peak*np.log(t) - t
        f = r_peak*np.log(r) - r
        return f - const

    def __str__(self):
        e, beta, mu = self.altpars()
        return 'Poisson: mu,beta= %.0f, %.0f' %( mu, beta)
    @property
    def flux(self):
        return max(self.p[0], 0)
    
    @property
    def errors(self):
        return self.find_delta()
        
    @property
    def ts(self):
        return 0 if self.flux<=0 else (self(self.flux)-self(0))*2.0
    
    def altpars(self):
        """ compute alternate parameters """
        e = abs(self.p[1])
        beta = e * abs(self.p[2])
        mu =   e * self.p[0] + beta
        return e,beta,mu
    
    def find_delta(self,delta_logl=.5):
        """Find points where the function decreases by delta from the max"""
        smax = max(0,self.p[0])
        ll_max = self(smax)
        ll_zero = self(0)
        func = lambda s: ll_max-self(s)-delta_logl
        if ll_max-ll_zero<delta_logl:
            s_low = 0
        else:
            #s_low = optimize.bisect(func,0,smax,xtol=.01*smax)
            s_low = optimize.brentq(func,0,smax,xtol=1e-17)
        if smax>0:
            s_high = smax*10
        else:
            s_high = 1e-15
        while func(s_high)<0: s_high*=2
        #s_high = optimize.bisect(func,smax,s_high,xtol=.01*smax)
        s_high = optimize.brentq(func,smax,s_high,xtol=1e-17)
        if not np.all(np.isreal([s_low,s_high])):
            print('Could not find two roots!')
            return None
        return (s_low,s_high)

    def cdf(self, flux ):
        """ cumulative pdf, from flux=0, according to Bayes
        uses incomplete gamma function for integrals. (Note that the scipy function
        is a regularized incomplete gamma function)
        """
        e, beta, mu = self.altpars()
        offset = special.gammainc(mu+1, beta) # Bayes offset if beta>0
        return (special.gammainc(mu+1, beta+flux*e)-offset)/(1-offset)

    def cdfc(self,flux):
        """ complementary cumulative cdf: 1-cdf(flux)"""
        e, beta, mu = self.altpars()
        return special.gammaincc(mu+1, beta+flux*e)/special.gammaincc(mu+1, beta)

    def cdfinv(self, pval):
        """ return the inverse of the cdf 
         pval : float
        """
        e,beta,mu = self.altpars()
        gbar = lambda x : special.gammainc(mu+1, beta+x)
        chatinv = lambda pv : special.gammaincinv(mu+1, pv+gbar(0)*(1-pv))-beta
        return  chatinv(pval)/e
        
    def cdfcinv(self, pvalc): 
        """ return the inverse of cdfc = 1-cdf
        useful for limits with very low probablity
        pvalc : float
            1-pval: zero corresponds to infinite flux
        """
        e,beta,mu = self.altpars()
        if e==0: return np.nan
        gcbar = lambda x : special.gammaincc(mu+1, beta+x)
        #cchat = lambda x : gcbar(x)/gcbar(0)
        cchatinv = lambda pv : special.gammainccinv( mu+1, pv*gcbar(0) )-beta
        return cchatinv(pvalc)/e

    def percentile(self, limit=0.95):
        """Left for compatibility: use cdfinv or cdfcinv
        """
        return self.cdfinv(limit)
        
    def pts(self):
        return 0 if self.flux<=0 else (self(self.flux)-self(0))*2.0
 
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
        """ evalate the likelihood function (*not* log likelihood), normalized to 1 at peak"""
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
            fig=plt.figure(fignum)
            axes = plt.gca()
        else: fig = plt.gcf()
        dom = np.linspace(0, self.maxflux, 101)
        if 'lw' not in kwargs or 'linewidth' not in kwargs: kwargs.update(lw=2)
        axes.plot(dom, self(dom),'-', **kwargs)
        color = kwargs.get('color', 'b')
        a,b = self.errors()
        if b==0: b=self.upper_limit()
        axes.axvspan(a,b, color=color, alpha=0.25)
        axes.grid(True)
        axes.set_ylim((-0.1, 1.1))
        return fig
    @staticmethod
    def test(mean=1.0, sigma = 0.1, guess=1.0):
        def fun(x):
            return -0.5*((x-mean)/sigma)**2
        t =  LogLikelihood(fun, guess=guess)
        return t
         
class PoissonFitter(object):
    """ Helper class to fit a log likelihood function to the Poisson
    
    Try a case with only a limit
    >>> pf = PoissonFitter( Poisson([-1, 1.0, 5.]) )
    >>> np.array(pf.fit()).round(3)
    array([-0.988,  0.996,  4.942])
    
    And one with a peak
    >>> pf2 = PoissonFitter( Poisson([50., 1., 10.]))
    >>> np.array(pf2.fit()).round(3)
    array([ 50.,   1.,  10.])
    
    """
    def __init__(self, func, scale=1., tol=0.02):
        """
        parameters
        ----------
        func : function of one parameter
        scale: float
            estimate for scale to use
        tol : float
            absolute tolerance for fit, within default domain out to delta L of 4
        """
        self.func = func
        self.smax = self.find_max(scale)
        # determine values of the function corresponding to delta L of 0.5, 1, 2, 4
        # depending on how peaked the function is, this will be from 5 to 8 
        # The Poisson will be fit to this set of values
        dom = set()
        for delta in (0.5, 1.0, 2.0, 4.0):
            a,b = self.find_delta(delta, scale, xtol=tol*1e-2)
            dom.add(a); dom.add(b)
        self.dom = np.array(sorted(list(dom)))
        self.fit()
        self.maxdev=self.check(tol)[0]
        
    @property
    def poiss(self):
        return self._poiss
        
    def __call__(self, x):
        if hasattr(x, '__iter__'):
            return map(self.func, x)
        return self.func(x)

    def find_max(self, scale):
        """Return the flux value that maximizes the likelihood.
        """
        t= optimize.fmin(lambda s: -self.func(s), scale,disp=False)[0]
        return t if t>0 else 0
    
    def find_delta(self, delta_logl=.5, scale=1.0, xtol=1e-5):
        """Find positive points where the function decreases by delta from the max
        """
        ll_max = self(self.smax)
        ll_zero = self(0)
        func = lambda s: ll_max-self(s)-delta_logl
        if ll_max-ll_zero<delta_logl:
            s_low = 0
        else:
            s_low = optimize.brentq(func,0, self.smax, xtol=xtol)
        if self.smax>0:
            s_high = self.smax*10
        else:
            s_high = scale
        while func(s_high)<0: s_high*=2
        s_high = optimize.brentq(func,self.smax,s_high, xtol=xtol)
        if not np.all(np.isreal([s_low,s_high])):
            print '%s.find_delta Failed to find two roots!' % self.__class__.__name__
            return None
        if s_high==s_low:
            print '%s.find_delta Failed to find high root with delta=%.1f: %s' % (self.__class__.__name__,delta_logl,s_high)
        return (s_low,s_high)

    def fit(self, mu=30, beta=5):
        """Do the fit, return parameters for a Poisson constructor
        mu, beta: initial parameters for the fit if the peak is positive
        """
        smax = self.smax
        if smax>0:
            # function to fit has positive peak. Fit the drived parameters mu, beta
            cod = self(self.dom)-self(smax)
            #print 'smax=%.2f, w(%s)=%s' % (smax, self.dom, cod)
            def fitfunc(p):
                mu,beta=p
                e=(mu-beta)/smax; b = beta/e
                self._poiss = Poisson([smax, e,b])
                r = self._poiss(self.dom)-cod
                #print'f(%.3f,%.3f): %s' % (mu,beta,r)
                return r
            mu,beta =  optimize.leastsq(fitfunc,[mu,beta], maxfev=10000)[0]
            e = (mu-beta)/smax; b = beta/e
            return [smax, e, b]
        else:
            # maximum is at zero, so only limit.
            x=self.dom; y=self(x)
            # exposure factor estimated from asymptotic behavior
            big= x[-1]*1e3; e = -self(big)/big; 
            # preliminary fit to the quadratic coeficients and estimate parameters from the linear and 2nd order
            pf = polyfit(x,y, 2)
            b,a = pf[:2] 
            beta = -e*(e+a)/(2.*b)
            mu = beta*(1+a/e)
            smax= (mu-beta)/e; b= beta/e 
            # now fit the Poisson with e fixed to the asym. estimate.
            cod = self(self.dom)-self(0)
            pinit = [smax,  b]
            def fitfunc(p):
                self._poiss = Poisson([p[0], e, p[1]])
                return self._poiss(x) -cod
            t =  optimize.leastsq(fitfunc,  pinit,  maxfev=10000)[0]
            return self._poiss.p

    def check(self, tol=0.02):
        deltas = np.array(map( lambda x: self.func(x)-self._poiss(x), self.dom))
        t = np.abs(deltas-deltas.mean()).max()
        if t>tol: raise Exception('PoissonFitter: maximum deviation, %.2f > tolerance, %s'%(t,tol))
        return t, deltas

 
class MultiPoiss(object):
    """
    Manage a source, and moultiple months
    
    """
    def __init__(self, name, logls):
        """ 
           name: text
           logls: dictionary: keys are
                months: list of Poisson objects
                flux: fit flux at reference energy
                e0: reference energy
        """
        self.name = name
        self.months = logls['months']
        self.e0 = logls['e0']
        self.norm = logls['flux']
        self.loglike = LogLikelihood(self)

    def __str__(self):
        s =  '\nLight curve for %s, flux=%.2e at %.2f GeV\n' % (self.name, self.norm, self.e0/1e3)
        s += (' '+5*'%-10s') % tuple('month peak low high 95% '.split()) 
        fluxes = self.fluxes
        errors = self.errors
        months = np.arange(len(fluxes))
        f95    = np.array([m.percentile() for m in self.months])
        t=zip( fluxes/self.norm,  errors[0,:]/self.norm , errors[1,:]/self.norm, f95/self.norm)
        for i, q in enumerate(self.months):
            row = t[i]
            fmt = '\n%5d' + len(row)*'%10.2f'
            s += fmt % ( (i+1,)+row)
        s += '\nlog likelihood difference =%.1f for %d months: prob=%.1e\n'\
                %( self.varindex, len(self.months), self.survival_prob)
        return s+'\n'
        
    def __call__(self, nflux):
        """ combined likelihood for given flux, relative to norm """
        t = np.array([p(nflux*self.norm) for p in self.months]).sum()
        return t
        
    @property
    def fluxes(self):
        return np.array([max(ps.p[0],0) for ps in self.months])
    @property
    def errors(self):
        t=[]
        for ps in self.months:
            t.append(ps.find_delta())
        return np.array(t).T
    @property
    def flux95(self):
        return np.array([ps.percentile() for ps in self.months])
    
    @property
    def varindex(self):
        """ The variability index, the value of the combined likelihood at the peak is
        the log of the likelihood ratio: the doubled value should be distributed like 
        a Chi-squared distribution of len(months)-1 ndf
        """
        return -self(self.loglike.maxl)
    @property
    def survival_prob(self):
        """ calculate the survival probability, for the log likelihood ratio assuming distributed as chi**2 with N-1 deg of f"""
        x = self.varindex
        return stats.chi2.sf(2*x, len(self.months)-1 )
   
    def plot_likelihoods(self, ymin=0, xmax=None):
        import pylab as plt
        if xmax is None: xmax= self.norm*10.
        dom = np.linspace(0, xmax, 100)
        fig , axs = plt.subplots(7,6, figsize=(14,18))
        for i,ax in enumerate(axs.flatten()):
            poiss = self.months[i]
            w = np.exp(poiss(dom))
            sig = poiss.find_delta()
            ax.plot(dom, w, '-')
            ax.plot(sig, [np.exp(poiss(x)) for x in sig],'-o', color='r',lw=2)
            ax.plot(dom, poiss.cdf(dom), 'g')
            plt.setp(ax, xscale='linear', xlim=(0,xmax), ylim=(ymin,1.0))
            ax.axvline(poiss.p[0], color='r')
            f95 = poiss.percentile()
            ax.plot([f95,f95], [ymin, np.exp(poiss(f95))], '-or', lw=2)
            ax.text(0.5*xmax, 0.9, '%d'%(i+1,), fontsize=8)

            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
        return fig

             
        
    def plot_light_curve(self,  fignum=11, axes=None, ymin=0.1, ymax=10,  **kwargs):
        """
        Make a nice light curve plot
        """
        import pylab as plt
        
        def point(axis, x, yvals, xerr=0, TS=None, arrow_size=0.4, **kwargs):
            """
                x, x-value
                yvals - a tuple (y, yl,yu, y95)
            """
            y, yl,yu,y95 = yvals[:4]
            a,b = axis.get_ylim(); 
            fmt =kwargs.pop('fmt','o')
            if yl>2*a or TS>9 and TS is not None: 
                # lower limit is clearly above lower limit: plot the point
                axis.plot([x], [y], fmt, **kwargs)
                axis.plot([x-xerr,x+xerr], [y,y], '-', **kwargs)
                axis.plot([x,x], [yl, yu], '-', **kwargs)
                #axis.errorbar( [x], [y], xerr =xerr, yerr=[[yu-y],[y-yl]] if yu>y else 0, fmt=fmt, capsize=0,**kwargs)
            elif y95>0:
                # will be an upper limit: plot symbol
                axis.plot( [x-xerr,x+xerr], [y95,y95], **kwargs)
                dx = arrow_size
                dy = arrow_size*(b-a) 
                if axis.get_yscale()=='linear':
                    ya,yb,yc = y95, y95-dy, y95-2*dy
                else:
                    ya, yb, yc = y95, y95/1.2, y95/1.2**2
                kwargs.pop('ms', None)
                axis.fill( [x,  x, x+dx, x,  x-dx, x], 
                           [ya, yb,yb,   yc, yb,  yb], **kwargs)
            else:
                axis.plot( [x-xerr,x+xerr], [y95,y95], '*', **kwargs)

        norm = self.norm
        e0 =  self.e0
        months =self.months
        name = self.name
        oldlw = plt.rcParams['axes.linewidth']
        plt.rcParams['axes.linewidth'] = 2
        if axes is None:
            plt.close(fignum)
            fig=plt.figure(fignum, figsize=(8,3)); 
            axes =fig.add_axes((0.22,0.15,0.75,0.72))
        if 'color' not in kwargs: kwargs['color']='k'
        xmax = len(months)+1
        plt.setp(axes, xlim=(0,xmax+1), ylim=(ymin,ymax), yscale='log', 
                    xlabel='months',ylabel='flux /  reference')
        #axes.set_autoscale_on(False)
        for i,month in enumerate(months):
            yvals = [ month.flux] + list(month.errors) + [month.percentile(0.95)] 
            point(axes, i+1, np.array(yvals)/norm, 0.5, TS=None,  **kwargs)
        axes.axhline(1.0, color='k', lw=1)
        axes.set_title('%s' % name, fontsize=12)
        axes.text(0.2*xmax,1.1*ymin, 'reference flux=%.2e at %d (%.1f eV/cm^2/s)' \
            % (norm, e0, norm*e0**2*1e6), fontsize=10)
        yl,yu = self.loglike.errors()
        axes.axhspan(yl,yu, color='r', alpha=0.5)
        axes.grid(True)
        plt.rcParams['axes.linewidth'] = oldlw
if __name__ == "__main__":
    print __doc__
    import doctest
    doctest.testmod()
