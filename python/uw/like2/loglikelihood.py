"""Tools for parameterizing log likelihood curves.

Author(s): Eric Wallace, Matthew Kerr, Toby Burnett
$Header: /nfs/slac/g/glast/ground/cvs/pointlike/python/uw/like2/loglikelihood.py,v 1.19 2013/12/26 17:32:13 burnett Exp $
"""
__version__ = "$Revision: 1.19 $"

import numpy as np
from scipy import optimize, special, polyfit, stats


class Poisson(object):
    """log of the three-parameter Poisson-like function use to represent the flux likelihood
    parameters are, in order:
        sp : flux at peak, if positive; if negative, there is only a limit
        e  : normalization factor to convert flux to equivalent counts. must be >0
        b  : background flux: must be >=0
        
        This parametrization used is equivalent to that described in William
        Tompkins' thesis (arxiv: astro-ph/0202141) and Nolan, et al., 2003,
        ApJ 597:615:627. The functional form is that of a Poisson distribution

        
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
        return 'Poisson: mu,beta= %.1f, %.1f' %( mu, beta)
    
    def __repr__(self):
        e, beta, mu = self.altpars()
        return '%s.%s: mu,beta=%.1f, %.1f' % (self.__module__, self.__class__.__name__, mu,beta)
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
        """ return alternate parameters: e, beta, mu """
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
        if limit>=0.95: return self.cdfcinv(1-limit)
        return self.cdfinv(limit)
        
    def pts(self):
        return 0 if self.flux<=0 else (self(self.flux)-self(0))*2.0
 
         
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
    def __init__(self, func, scale=1., tol=0.05, delta=0.1, dd=-0.1):
        """
        parameters
        ----------
        func : function of one parameter
        scale: float
            estimate for scale to use
        tol : float
            absolute tolerance in probability amplitude for fit, within default domain out to delta L of 4
        """
        self.func = func
        #first check derivative at zero flux - delta is sensitive
        s = self.wprime = (func(delta)-func(0))/delta
        self.smax = self.find_max(scale) if s>0 else 0.
        # determine values of the function corresponding to delta L of 0.5, 1, 2, 4
        # depending on how peaked the function is, this will be from 5 to 8 
        # The Poisson will be fit to this set of values
        dlist = np.array([0.5, 1.0, 2.0, 4.0])
        if s < dd:
            # large negative derivative: this will be just an exponential
            if s < -100: s=-100. #cut off for now
            self.dom = - dlist/s
            self._poiss= Poisson([-1, -s, 1])
            self.maxdev=0
            return #no test in this case
        else:
            dom = set()
            for delta in dlist:
                a,b = self.find_delta(delta, scale, xtol=tol*1e-2)
                dom.add(a); dom.add(b)
            self.dom = np.array(sorted(list(dom)))
            self.fit()
        self.maxdev=self.check(tol)[0]
        
    def __repr__(self):
        return '%s.%s : wprime=%.3e maxdev=%.2f, %s' % (self.__module__,self.__class__.__name__,
            self.wprime, self.maxdev,  str(self._poiss))
    
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
        if self.func(0) > self.func(scale/10.):
            return 0
        r= optimize.fmin(lambda s: -self.func(s), scale, ftol=0.01, xtol=0.01, 
                disp=False, full_output=True, retall=True)
        t = r[0][0]
        #if t==scale:
        #    raise Exception('Failure to find max value: %s' % list(r))
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
        while func(s_high)<0 and s_high<1e6: 
            s_high*=2
        s_high = optimize.brentq(func,self.smax,s_high, xtol=xtol)
        if not np.all(np.isreal([s_low,s_high])):
            msg = '%s.find_delta Failed to find two roots!' % self.__class__.__name__
            print msg
            raise Exception( msg)
        if s_high==s_low:
            msg= '%s.find_delta Failed to find high root with delta=%.1f: %s' % (self.__class__.__name__,delta_logl,s_high)
            print msg
            print 'wprime: %.3e' % self.wprime
            raise Exception(msg)
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
            mu,beta =  optimize.leastsq(fitfunc,[mu,beta], ftol=1e-6,xtol=1e-6, maxfev=10000)[0]
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
            t =  optimize.leastsq(fitfunc,  pinit,  xtol=1e-6,ftol=1e-6, maxfev=10000)[0]
            return self._poiss.p

    def check(self, tol=0.05):
        offset = self(self.smax)
        deltas = np.array(map( lambda x: np.exp(self.func(x)-offset)-np.exp(self._poiss(x)), self.dom))
        t = np.abs(deltas).max()
        if t>tol: raise Exception('PoissonFitter: maximum deviation, %.2f > tolerance, %s'%(t,tol))
        return t, deltas
    
    def plot(self, ax=None):
        """Return a figure showing the fit"""
        import matplotlib.pyplot as plt
        xp = self.dom
        x = np.linspace(0, xp[-1])
        if ax is None:
            fig, ax = plt.subplots(figsize=(3,3))
        else: fig = ax.figure
        pfmax = self(self.smax)
        ax.plot(x, np.exp(self(x)-pfmax), '-', label='Input')
        ax.plot(xp, np.exp(self._poiss(xp)), 'o', label='approx')
        ax.legend(prop =dict(size=8) )        ax.set_xticks([0, xp[-1]])
        ax.grid()
        return fig
        

 
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
