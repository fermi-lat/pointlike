"""Tools for parametrizing log likelihood curves.

Author(s): Eric Wallace, Matthew Kerr
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/loglikelihood.py,v 1.2 2012/07/07 14:41:35 burnett Exp $
"""
__version__ = "$Revision: 1.2 $"

import numpy as np
from scipy import optimize, special
from uw.utilities import decorators

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
    
    where the peak is at x=mu, if mu>beta, else at x=0
    
    >>> w = Poisson([10.0, 1.0, 5.0])
    >>> map(w, np.linspace(0,20,5))
    [-6.4791843300216456, -1.0819766216224593, 0.0, -0.68476891322328726, -2.3376156435101407]
    >>> w.find_delta()
    (6.452927189667698, 14.213246175302684)
    >>> map(w.cdf, np.linspace(0,20,5))
    [0.0, 0.048674754021320932, 0.43187121882311852, 0.84347606391862051, 0.97770544018425898]
    >>> w.percentile()
    18.097300356339645
     
    >>> w = Poisson([-5, 1.0, 10.0]) 
    >>> map(w, np.linspace(0, 10, 5))
    [0.0, -1.3842822434289523, -2.9726744594591796, -4.7019210603228885, -6.5342640972002766]
    >>> map(w.cdf, np.linspace(0, 10, 5))
    [0.0, 0.77904655517623111, 0.95837535584402034, 0.9930192704933527, 0.99892810898968387]
    >>> w.percentile()
    4.734595152471901

    
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
        uses incomplete gamma function for integrals
        """
        e = abs(self.p[1])
        beta = e * abs(self.p[2])
        mu =   e * self.p[0] + beta
        norm = 1-special.gammainc(mu+1, beta)
        q = special.gammainc(mu+1, beta + e*flux )
        return (q-1)/norm  +1.

    def percentile(self, limit=0.95):
        """return percentile for given limit, default 95%"""
        e = abs(self.p[1]) 
        f = lambda x : self.cdf(x/e)-limit
        xmax = self.find_delta(10)[1]*e
        ret = optimize.brentq(f, 0, xmax) 
        assert abs(f(ret))<1e-2, 'percentile failed fit: %.1f %.1f %.1e' % (ret, f(ret), xmax)
        return ret/e
        
    
if __name__ == "__main__":
    print __doc__
    import doctest
    doctest.testmod()
